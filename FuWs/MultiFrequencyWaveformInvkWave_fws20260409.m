% ===================== Waveform Inversion + Virtual Elements =========GIF=================
%  版本: FreqVE
%  新增特性: 按频率动态切换虚拟阵元（VE_freqMode: always / threshold / manual）
%  修改说明:
%    1. 顶部新增 VE_freqMode / VE_freqThresh / VE_manualFlags 控制参数
%    2. VE代码块改为预计算 geo_orig / geo_ve 两套几何，不再覆盖全局变量
%    3. 主循环内按 VE_flags(f_idx) 动态选择当前套装，跨套装切换时强制重置CG
%    4. GIF 预创建窗口，复用图像句柄，大幅减少每帧开销
%    5. [修改] 已删除 PP 模块及其相关代码
%    6. [新增] 方向A2：VE通道置信度加权（beta_ve），降低虚拟Rx插值偏差对梯度的污染
%    7. [新增] 方向D ：照明补偿预条件（enableIllumComp），均衡稀疏阵角向不均匀照明
%    8. [新增] 失配函数：misfitType 支持 'L2'（标准）和 'HV_polar'（极坐标分解）
%             HV_polar 将复数 phasor 分解为幅度残差（→衰减）和相位残差（→声速），
%             分别用 alpha_hv / alpha_hv_atten 加权，缓解 L2 的周期跳跃敏感性。
%             参考：Neumann & Yang, arXiv:2505.01817 (2025)，频域FWI + USCT验证。
%    9. [修改] 方向E ：低频先验正则化升级为三阶段版本（enableLFPrior），
%             对应师兄时域FWI"均匀分三段低频复用，效果较好"的频域等价实现：
%               阶段1（低频，< f_stage1_cutoff）：无先验约束，自由重建背景声速
%               阶段2（中频，f_stage1~f_stage2）：约束拉向阶段1末尾参考模型
%               阶段3（高频，> f_stage2_cutoff）：约束拉向阶段2末尾参考模型
%             各阶段内 lambda_k = lambda_stage·(f_stageN_cutoff/f)²，频率越高约束越弱。
%             相对归一化消除量纲不匹配（见主循环内注释）。
%   10. [新增] 方向F ：各向同性TV正则化（enableTV），Huber平滑近似，
%             抑制128稀疏阵梯度相干增强产生的同心环状伪影，保留真实声速边界。
%             [修改] 相对归一化消除量纲不匹配；lambda_tv增至2e-3增强环状抑制效果。
%   11. [新增] 方向G ：极坐标径向TV正则化（enablePolarReg），仅对径向梯度分量
%             施加TV惩罚，比各向同性TV更定向地抑制同心环，对真实边界干扰更小。
%             [修改] 相对归一化消除量纲不匹配；lambda_polar增至2e-3与方向F协同。
%   12. [修改] 方向E/F/G 均改用相对幅值归一化，lambda_* 直接表示占数据梯度比例。
%   13. 其余算法逻辑（Helmholtz、CG）与原版一致
% ===========================================================================================
clear; clc;
scriptTimer = tic;   % 整个脚本总时长（含预计算、主反演、保存）

% Add Functions to Path
addpath(genpath('Functions'));

% ------------------ 用户：核心开关 ------------------
enableVirtualElements = true;          % true=启用虚拟阵元；false=原始128阵元基线
useVirtualTx          = false;         % true=虚拟阵元也参与发射；false=仅真实Tx
virtualElementsMode   = 'unilateral';  % 'unilateral' = 单侧交替（128→256）
                                       % 'bilateral'  = 两侧对称（128→384）
verTag                = 'V3_5';        % 保存版本号
% ----------------------------------------------------------------------

% ------------------ fws: Near-end exclusion ------------------
exclMode  = 'frac';   % 'angle' | 'frac'
frac_excl = 1/4;
% -------------------------------------------------------------

% ------------------ fws: 失配函数选择 ------------------
% 'L2'       : 标准最小二乘，逐通道复数 phasor 残差的 L2 范数（原版）
% 'HV_polar' : 极坐标分解失配（HV度量简化形式），将复数 phasor 分解为：
%              ① 幅度残差 ΔA = |p_sim| - |p_obs|  → 主要约束衰减
%              ② 相位残差 Δφ = angle(p_sim·conj(p_obs)) ∈ (-π,π] → 主要约束声速
%              伴随源 = (1-α)·ΔA·exp(iφ_sim) + α·Δφ·i·exp(iφ_sim)
%              其中 exp(iφ_sim) 为幅度方向，i·exp(iφ_sim) 为正交相位方向。
%              α=alpha_hv（SoS阶段）或 alpha_hv_atten（衰减阶段）
% 物理依据：L2 对相位误差与幅度误差的梯度方向不加区分；
%           HV_polar 用相位分量提供指向真实声速的切向修正，
%           增大凸性盆地，抑制周期跳跃（参考 Neumann & Yang 2025）。
misfitType = 'HV_polar';   % 'L2' | 'HV_polar'

% alpha_hv：SoS 更新阶段的相位权重
%   0.0 = 纯幅度（退化为 L2 幅度版）
%   1.0 = 纯相位（对周期跳跃最鲁棒，但完全忽略幅度匹配）
%   建议起点：0.7（强调相位约束声速；消融可测试 [0.3, 0.5, 0.7, 0.9]）
alpha_hv = 0.5;

% alpha_hv_atten：衰减更新阶段的相位权重
%   衰减主要由幅度约束，应降低相位权重；建议起点：0.3
alpha_hv_atten = 0.1;
% -------------------------------------------------------

% ------------------ fws: VE通道置信度加权（方向A2）------------------
% 物理依据：虚拟Rx通道由插值生成，不是独立测量，其残差包含插值误差引入的偏差。
%           对真实Rx保持单位权重，对虚拟Rx施加 beta_ve < 1 的降权，
%           使伴随源注入中虚假相关成分减少，梯度更准确地反映真实物理约束。
% 取值范围：1.0 = 等权（与原版等价）；0.0 = 仅真实Rx参与梯度（最保守）
% 建议起点：0.5，消融实验可测试 [0.0, 0.3, 0.5, 0.7, 1.0]
beta_ve = 0.9;
% -------------------------------------------------------------------

% ------------------ fws: 照明补偿预条件（方向D）------------------
% 物理依据：128稀疏阵的角向照明极不均匀，靠近发射源的侧面照明强、
%           远离发射源的区域照明弱，导致梯度图像出现角向不均匀环状偏差。
%           用各发射源正向波场强度之和 illum = Σ|u_i|² 估计伪Hessian对角，
%           对梯度做归一化预条件，等效改善目标函数曲率。
% eps_illum：归一化照明的截断参数（相对值），防止零照明区除零。
%            值越大→正则化越强→抑制越多高空间频率→图像越平滑；
%            值越小→更接近真实照明补偿→可能放大低照明区的噪声。
%            建议起点：0.05，可测试 [0.01, 0.05, 0.10, 0.20]
enableIllumComp = true;    % true=开启照明补偿；false=关闭（退回标准梯度）
eps_illum       = 0.10;    % 照明补偿正则化参数（相对归一化照明的截断值）
% ----------------------------------------------------------------

% [修改] ------------------ fws: 三阶段低频先验正则化（方向E）------------------
% 物理依据：对应师兄时域FWI"均匀分三个阶段进行低频复用，效果较好"的频域等价。
%           时域三段方案：Low(0.3~0.6 MHz)→Mid(0.625~0.95 MHz)→High(0.975~1.25 MHz)
%           频域实现：在阶段1→2和阶段2→3的边界处各保存一次参考模型，
%           后续阶段施加Tikhonov约束把当前模型拉向前一阶段的末尾结果：
%             阶段1（低频）：无先验约束，自由重建背景声速
%             阶段2（中频）：约束拉向 VEL_stage1_ref（阶段1末尾）
%             阶段3（高频）：约束拉向 VEL_stage2_ref（阶段2末尾）
%           约束强度：lambda_k = lambda_stage · (f_stageN_cutoff / fDATA(f_idx))²
%           → 各阶段内频率越高约束越弱，阶段切换时自动切换参考锚点。
%           [修改] 相对归一化：同之前单阶段版本，lambda_stage 含义为占数据梯度的比例。
%                   三阶段版本因阶段2/3分别有独立参考，约束更精准，副作用更小。
% f_stage1_cutoff：阶段1/2分界 [Hz]，对应频率表第一段末（≈0.475 MHz）
% f_stage2_cutoff：阶段2/3分界 [Hz]，对应频率表第二段末（≈0.725 MHz）
% lambda_stage：各阶段约束强度基准（消融可测试 [5e-4, 1e-3, 5e-3]）
enableLFPrior    = true;     % true=开启三阶段LF先验；false=关闭
lambda_stage     = 1e-3;     % 各阶段约束强度（归一化后：0.1%数据梯度；消融可测试 [5e-4, 1e-3, 5e-3]）
f_stage1_cutoff  = 0.475e6;  % 阶段1/2分界 [Hz]（≈Low段末，与 frac_shift_f_split 对齐）
f_stage2_cutoff  = 0.850e6;  % 阶段2/3分界 [Hz]（≈Mid段末，对应 fDATA_SoS 第二组末尾）
% -----------------------------------------------------------------------

% [修改] ------------------ fws: TV正则化（方向F：各向同性Huber-TV）-----------
% 物理依据：128稀疏阵梯度相干增强在特定环形轨迹上产生同心环状伪影，
%           其空间特征为高频振荡（径向方向）。
%           各向同性TV对所有方向的梯度振荡施加稀疏性约束：
%             R_TV = lambda_tv · ∫ |∇v|_Huber dΩ
%           梯度方向（作为伴随项添加到 gradient_img）：
%             reg_grad = lambda_tv · (-div(∇v / |∇v|_Huber))
%           Huber平滑（eps_tv 参数）避免 |∇v|→0 时梯度奇异。
%           真实声速边界处 |∇v| 大，TV惩罚自动降权，边界保留性优于Laplacian。
%           [修改] 相对归一化消除量纲不匹配；lambda_tv 提高至 2e-3 以增强环状伪影抑制。
%           与方向G（极坐标径向TV）同时开启，两者协同压制同心环。
% lambda_tv：TV强度（归一化后含义；消融可测试 [1e-3, 2e-3, 5e-3]）
% eps_tv：Huber平滑参数（建议 1e-3~1e-2；越小越接近纯TV，越大越接近L2平滑）
enableTV  = true;     % true=开启各向同性TV；false=关闭
lambda_tv = 2e-3;     % TV强度（归一化后：2.0%数据梯度；消融可测试 [1e-3, 2e-3, 5e-3]）
eps_tv    = 1e-3;     % Huber平滑参数（建议 1e-3~1e-2）
% -----------------------------------------------------------------------

% [修改] ------------------ fws: 极坐标径向TV（方向G：同心环定向抑制）---------
% 物理依据：同心环特征为径向梯度剧烈、切向梯度光滑，与真实声速边界的
%           各向同性梯度特征不同。仅对径向梯度分量施加TV惩罚：
%             reg_grad = lambda_polar · (-div_r( (∂v/∂r) / |∂v/∂r|_Huber ))
%           比方向F更定向地抑制同心环，对真实边界的干扰更小。
%           [修改] 相对归一化消除量纲不匹配；lambda_polar 提高至 2e-3，
%           与方向F（各向同性TV）同时开启，协同抑制环状伪影。
%           两者分工：方向F从各向同性角度抑制整体振荡，
%                     方向G专门针对径向周期结构（同心环的定向特征）强化惩罚。
% lambda_polar：径向TV强度（消融可测试 [1e-3, 2e-3, 5e-3]）
% eps_polar_reg：Huber平滑参数（同 eps_tv 数量级）
enablePolarReg = true;     % true=开启极坐标径向TV；false=关闭
lambda_polar   = 2e-3;     % 径向TV强度（归一化后：2.0%数据梯度；消融可测试 [1e-3, 2e-3, 5e-3]）
eps_polar_reg  = 1e-3;     % Huber平滑参数（建议 1e-3~1e-2）
% -----------------------------------------------------------------------

% ------------------ fws: 高频更新可信度监控与门控 ------------------
% 目标：避免“高频把低波数背景误差投影到高波数环状伪影”。
% 核心监控量：
%   1) ΔΦ_k = Φ_{k-1} - Φ_k：数据失配是否真实下降
%   2) τ_k  ：平均时移残差（由相位残差换算）是否继续下降
%   3) R_k  ：背景环带的径向伪影比（越大越像同心环）
% 门控策略：
%   - 若 ΔΦ_k <= 0 或 τ_k 上升：缩小步长（alpha_gate）
%   - 若 ΔΦ_k 很小且 R_k 上升：进一步缩步长 + 增强径向TV（lambda_polar_gate）
enableTrustGate        = true;
phi_small_drop_ratio   = 1e-3;  % “失配下降很小”的相对阈值
tau_rise_ratio         = 1.01;  % τ_k > 1.01*τ_{k-1} 视为上升
ring_rise_ratio        = 1.02;  % R_k > 1.02*R_{k-1} 视为环伪影增强
alpha_gate_badfit      = 0.50;  % 拟合退化时步长缩放
alpha_gate_ringgrowth  = 0.70;  % 环伪影增长时额外缩放
lambda_polar_gate_boost = 1.25; % 环伪影增长时径向TV增强倍率
lambda_polar_gate_max   = 3.00; % 径向TV最大增强上限（相对 lambda_polar）
% -----------------------------------------------------------------------

% ------------------ fws: 高频边界保护（抑制“中心污染外圈”） ------------------
% 现象：高频阶段若低波数背景仍有偏差，中心高对比分量可能通过有限孔径/建模误差
%      把更新能量“扩散”到外围背景，表现为外圈泛红/环纹。
% 处理：
%   1) 外环锚定正则（edge anchor）：把外环拉向低频参考，防止外圈被中心拖偏；
%   2) 外环步长阻尼（edge step damping）：高频阶段对外环更新做软抑制。
enableEdgeGuard    = true;
f_edge_guard_start = 0.85e6;  % 从高频段起生效（可与 f_stage2_cutoff 对齐）
edge_r_inner_ratio = 0.70;    % 外环保护起点（相对半径）
edge_r_outer_ratio = 0.92;    % 外环保护终点（相对半径）
lambda_edge_anchor = 2e-3;    % 外环锚定强度（相对归一化）
edge_step_damp_max = 0.60;    % 外环最大步长衰减比例（0~1）
edge_blend_base    = 0.08;    % 后处理外环回拉比例（每步）
edge_blend_max     = 0.25;    % 后处理外环回拉上限
edge_contam_rise_ratio = 1.03;% 外环污染判据（较上一步上升）
% -----------------------------------------------------------------------

% ------------------ GIF 保存 ------------------
% saveGIF=true 时，每次迭代后截取 figure(99) 帧并追加到 GIF
% 用于观察环状伪影从何频率段开始出现、VE 如何影响其演化
% 优化：主循环前预创建 figure(99) 和图像句柄，每帧只更新 CData，避免重建渲染链
saveGIF        = true;   % true=保存迭代 GIF；false=跳过
gif_delay      = 0.20;    % 每帧停留时间 [s]（建议 0.15~0.5）
gif_save_every = 1;       % 每隔多少次迭代保存一帧（1=每次，2=隔一次，减少文件大小）
% GIF 路径将在 result_dir 定义后自动生成（见主循环前）
% -----------------------------------------------

% Red-Blue colormap (blue=low, red=high)
% Soft Blue-White-Red diverging colormap
n_cmap = 256; n2 = n_cmap/2;
cmap_rb = [
    linspace(0.17, 1.00, n2)', linspace(0.51, 1.00, n2)', linspace(0.83, 1.00, n2)';
    linspace(1.00, 0.84, n2)', linspace(1.00, 0.18, n2)', linspace(1.00, 0.15, n2)'
];

% Load Extracted Dataset
filename = 'kWave_BreastCT_half'; % 'kWave_BreastCT', 'kWave_BreastMRI'
load(['D:\Document_ING_fws\WaveformInversionUST\Simulations\datasets\', filename, '.mat'], ...
    'xi_orig', 'yi_orig', 'C', 'atten', ...
    'time', 'transducerPositionsXY', 'full_dataset');

% === fws: 提前保存原始阵元数（避免后面 VE 扩展后丢失）===
orig_numElements = size(transducerPositionsXY, 2);
% ======================================================

% Assemble Full Dataset
numElements = size(transducerPositionsXY, 2); % 128
fsa_rf_data = single(full_dataset);
clearvars full_dataset;

% % Extract the Desired Frequency for Waveform Inversion
% fDATA_SoS      = (0.3:0.05:1.25)*(1e6);         % SoS-only
% fDATA_SoSAtten = (0.325:0.05:1.275)*(1e6);      % SoS+Atten
% fDATA          = [fDATA_SoS, fDATA_SoSAtten];   % All

% fws-Modified: Multi-stage step size with same overall frequency range
% Low freq (0.30–0.60 MHz): coarse 50 kHz steps
% Mid freq (0.63–0.90 MHz): medium 40 kHz steps
% High freq (0.93–1.25 MHz): fine 30 kHz steps
% Total range remains the same as original (0.3–1.25 MHz)
% fDATA_SoS = [ (0.30:0.025:0.475), (0.525:0.05:0.725), (0.80:0.075:1.25) ] * 1e6;
fDATA_SoS = [ (0.30:0.025:0.45), (0.50:0.05:0.70), (0.775:0.075:1.225), 1.250 ] * 1e6;
% fDATA_SoS = [ ...
%     (0.30:0.075:0.45), ...      % 3个：低频粗步长
%     (0.50:0.05:0.70),  ...      % 5个：中频正常步长
%     (0.725:0.025:0.825), ...    % 5个：ring onset 过渡区加密
%     (0.875:0.05:1.125), ...     % 6个：后段恢复中等步长
%     1.250 ] * 1e6;              % 1个：终点频率
% fDATA_SoS = [ (0.30:0.025:0.45), 0.350, (0.50:0.05:0.65), 0.350, (0.85:0.075:1.225), 0.350 ] * 1e6;
fDATA_SoSAtten = fDATA_SoS + 25e3;  % same +25 kHz shift as original
fDATA = [fDATA_SoS, fDATA_SoSAtten]; % all frequencies [Hz]

% [调试] 频率表校验
fprintf('[CHECK] numel(fDATA_SoS)=%d | SoS总迭代=%d | Niter预估=%d\n', ...
    numel(fDATA_SoS), numel(fDATA_SoS)*3, numel(fDATA_SoS)*3 + numel(fDATA_SoSAtten)*6);

% ------------------ 频率相关 frac_shift 调度 ------------------
% 物理依据：低频波场空间变化慢，小偏移插值误差小且已足够提供角向补偿；
%           高频波场空间变化快，需保持原有 frac_shift 以提供足够角向支撑，
%           但不应继续增大（插值相位误差 ∝ f × frac_shift × dtheta）。
% 因此采用"低频保守、高频标准"的两档调度，而非全程统一值。
frac_shift_lo      = 0.10;    % 低频段保守偏移（f < frac_shift_f_split）
frac_shift_hi      = 0.20;    % 高频段标准偏移（f >= frac_shift_f_split）
frac_shift_f_split = 0.45e6;  % 低/高频分界 [Hz]
% 注：两档各预计算一套 geo_ve，主循环内按频率自动选择
% ---------------------------------------------------------------

% ------------------ VE 频率开关 ------------------
% 'always'    : 所有频率均开VE
% 'threshold' : fDATA >= VE_freqThresh 的频率开VE
% 'manual'    : 由 VE_manualFlags 逐频指定
VE_freqThresh  = NaN;    % 仅 threshold 模式使用
% -------------------------------------------------

% ------------------ VE_manualFlags 在此处填写 ------------------
% 每次切换频率时交替开关 VE
VE_freqMode    = 'always';
VE_manualFlags = mod(0:numel(fDATA)-1, 2) == 1;  % 偶数索引关，奇数索引开（交替）
% 其他示例：
%   全程开VE：VE_freqMode = 'always';
%   阈值开VE：VE_freqMode = 'threshold'; VE_freqThresh = 0.7e6;
%   前半关后半开：VE_manualFlags = [false(1,numel(fDATA_SoS)), true(1,numel(fDATA_SoSAtten))];
% 注意：manual模式下 VE_manualFlags 长度必须等于 numel(fDATA)
% ---------------------------------------------------------------

% Iteration schedule
niterSoSPerFreq   = [3*ones(size(fDATA_SoS)),   3*ones(size(fDATA_SoSAtten))];
niterAttenPerFreq = [0*ones(size(fDATA_SoS)),   3*ones(size(fDATA_SoSAtten))];

% DTFT (not FFT)
DTFT = exp(-1i*2*pi*fDATA'*time)*mean(diff(time));

% Geometric TOF Based on Sound Speed
c_geom = 1540; % % Sound Speed [m/s] (原脚本注释写mm/us，但这里按原仓库逻辑使用即可)
transducerPositionsX = transducerPositionsXY(1,:);
transducerPositionsY = transducerPositionsXY(2,:);
geomTOFs = zeros(size(transducerPositionsXY,2));
for col = 1:size(transducerPositionsXY,2)
    geomTOFs(:,col) = sqrt((transducerPositionsX-transducerPositionsX(col)).^2 + ...
                           (transducerPositionsY-transducerPositionsY(col)).^2)/c_geom;
end
clearvars col transducerPositionsXY;

% Window Time Traces - Extract the Frequencies for Waveform Inversion
REC_DATA = zeros(numElements, numElements, numel(fDATA), 'single');
sos_perc_change_pre  = 0.05;
sos_perc_change_post = Inf;
twinpre  = sos_perc_change_pre  * max(geomTOFs(:));
twinpost = sos_perc_change_post * max(geomTOFs(:));

for tx_element = 1:numElements
    times_tx = geomTOFs(tx_element,:);
    [TIMES_TX, TIME] = meshgrid(times_tx, time);
    window = exp(-(1/2)*(subplus(TIME-TIMES_TX)/twinpost + ...
                         subplus(TIMES_TX-TIME)/twinpre).^2);
    REC_DATA(tx_element,:,:) = permute(DTFT*(window.*fsa_rf_data(:,:,tx_element)),[2,1]);
end
clearvars twin window fsa_rf_data time times_tx TIME TIMES_TX DTFT;

%% ===================== Waveform Inversion Grid =====================
% Create Sound Speed Map and Transducer Ring
dxi  = 0.3e-3; xmax = 120e-3;
xi   = -xmax:dxi:xmax; yi = xi;
Nxi  = numel(xi); Nyi = numel(yi);
[Xi, Yi] = meshgrid(xi, yi);
R_grid = sqrt(Xi.^2 + Yi.^2);
% Discretize transducer ring onto grid (coarse grid)
x_circ = transducerPositionsX;
y_circ = transducerPositionsY;
x_idx  = dsearchn(xi(:), x_circ(:));
y_idx  = dsearchn(yi(:), y_circ(:));
ind    = sub2ind([Nyi, Nxi], y_idx, x_idx);
msk    = zeros(Nyi, Nxi); msk(ind) = 1;

% -------- 备份原始粗网格（ Helmholtz / 索引全程使用这个）--------
xi_original = xi; yi_original = yi; dxi_original = dxi;


%% ===================== Parameters =====================
h = dxi; % [m]
g = 1; % (grid spacing in Y)/(grid spacing in X)
alphaDB = 0.0; % Attenuation [dB/(MHz*mm)]
alphaNp = (log(10)/20)*alphaDB*((1e3)/(1e6)); % Attenuation [Np/(Hz*m)]
ATTEN = alphaNp*ones(Nyi,Nxi); % Make Spatially Varying Attenuation [Np/(Hz m)]
% Solve Options for Helmholtz Equation
sign_conv = -1; % Sign Convention
a0 = 10; % PML Constant
L_PML = 9.0e-3; % Thickness of PML  

% Conversion of Units for Attenuation Map
Np2dB = 20/log(10);
slow2atten = (1e6)/(1e2); % Hz to MHz; m to cm

crange = [1350, 1600]; % For reconstruction display [m/s]
attenrange = 10*[-1,1]; % For reconstruction display [dB/(cm MHz)]

%% ===================== Phase Screen (for discretization error) =====================
% Compute Phase Screen to Correct for Discretization of Element Positions
% Times for Discretized Positions Assuming Constant Sound Speed
x_circ_disc = xi(x_idx); y_circ_disc = yi(y_idx);
[x_circ_disc_tx, x_circ_disc_rx] = meshgrid(x_circ_disc, x_circ_disc);
[y_circ_disc_tx, y_circ_disc_rx] = meshgrid(y_circ_disc, y_circ_disc);
geomTOFs_disc  = sqrt((x_circ_disc_tx-x_circ_disc_rx).^2 + ...
                      (y_circ_disc_tx-y_circ_disc_rx).^2)/c_geom;
% Time-of-Flight Error for Discretized Positions
geomTOFs_error = geomTOFs_disc - geomTOFs;

% Phase Screen and Gain Correction
PS = zeros(numElements, numElements, numel(fDATA));
for f_idx = 1:numel(fDATA)
    PS(:,:,f_idx) = exp(1i*sign_conv*2*pi*fDATA(f_idx)*geomTOFs_error);
end
REC_DATA = REC_DATA .* PS;% Phase Screen and Gain Corrections
clearvars PS geomTOFs_error geomTOFs_disc geomTOFs;

REC_DATA(isnan(REC_DATA)) = 0;

%% ===================== dwnsmp（放在预计算之前，两套共用）=====================
% Which Subset of Transmits to Use
dwnsmp = 1; % can be 1, 2, or 4 (faster with more downsampling)
            % NOTE: dwnsmp = 1 to get the results in the paper

%% ===================== 预计算两套几何 =====================
%  geo_orig      : 原始128阵元套装
%  geo_ve_lo     : VE扩展，frac_shift_lo（低频保守插值）
%  geo_ve_hi     : VE扩展，frac_shift_hi（高频标准插值）
%
%  每套结构体字段：
%    .numElements          阵元总数
%    .transducerPositionsX/Y  阵元位置
%    .x_idx / y_idx        网格离散索引
%    .ind                  线性索引
%    .elemInclude          近端排除矩阵 (numElements x numElements)
%    .isRealRx             真实接收标志 (1 x numElements)
%    .tx_include           发射阵元索引
%    .REC_DATA             切片后数据 (numel(tx_include) x numElements x numel(fDATA))
%    .frac_shift           本套装使用的插值偏移量（仅VE套装有此字段）

% ---- angle 模式共用的排除角 ----
% % Extract Subset of Times within Acceptance Angle
% numElemLeftRightExcl = 63;
% elemLeftRightExcl = -numElemLeftRightExcl:numElemLeftRightExcl;
% elemInclude = true(numElements, numElements);
% for tx_element = 1:numElements 
%     elemLeftRightExclCurrent = elemLeftRightExcl + tx_element;
%     elemLeftRightExclCurrent(elemLeftRightExclCurrent<1) = numElements + ...
%          elemLeftRightExclCurrent(elemLeftRightExclCurrent<1);
%     elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements) = ...
%         elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements) - numElements;
%     elemInclude(tx_element,elemLeftRightExclCurrent) = false;
% end
numElem_Ali     = 512;
numExclSide_Ali = 63;
theta_excl      = numExclSide_Ali * (2*pi/numElem_Ali); % ≈43.5 deg

% ============================================================
%  套装A：原始128阵元
% ============================================================
geo_orig.numElements          = orig_numElements;
geo_orig.transducerPositionsX = transducerPositionsX;
geo_orig.transducerPositionsY = transducerPositionsY;
geo_orig.x_idx                = x_idx;
geo_orig.y_idx                = y_idx;
geo_orig.ind                  = ind;
geo_orig.isRealRx             = true(1, orig_numElements);

numELRE_orig     = round(orig_numElements * frac_excl / 2);
elemExcl_orig    = -numELRE_orig:numELRE_orig;
elemInclude_orig = true(orig_numElements, orig_numElements);
theta_all_orig   = atan2(transducerPositionsY(:), transducerPositionsX(:));

switch lower(exclMode)
    case 'angle'
        for tx_element = 1:orig_numElements
            dtheta = theta_all_orig - theta_all_orig(tx_element);
            dtheta = atan2(sin(dtheta), cos(dtheta));
            elemInclude_orig(tx_element, abs(dtheta) <= theta_excl) = false;
        end
    case 'frac'
        for tx_element = 1:orig_numElements
            ec = elemExcl_orig + tx_element;
            ec(ec < 1) = orig_numElements + ec(ec < 1);
            ec(ec > orig_numElements) = ec(ec > orig_numElements) - orig_numElements;
            elemInclude_orig(tx_element, ec) = false;
        end
end
geo_orig.elemInclude = elemInclude_orig;
geo_orig.tx_include  = 1:dwnsmp:orig_numElements;

% Remove Outliers from Observed Signals Prior to Waveform Inversion
perc_outliers = 0.99; % Confidence Interval Cutoff
REC_DATA_orig_full = REC_DATA;
for f_idx = 1:numel(fDATA)
    slice = REC_DATA_orig_full(geo_orig.tx_include,:,f_idx);
    mag   = elemInclude_orig(geo_orig.tx_include,:) .* abs(slice);
    n_out = ceil((1-perc_outliers)*numel(mag));
    [~,idx_out] = maxk(mag(:), n_out);
    slice(idx_out) = 0;
    REC_DATA_orig_full(geo_orig.tx_include,:,f_idx) = slice;
end
geo_orig.REC_DATA = REC_DATA_orig_full(geo_orig.tx_include,:,:);

% ============================================================
%  套装B：Virtual Elements Scheme（按 frac_shift 各建一套）
%  单侧模式（unilateral）：每个真实阵元单侧交替增加1个虚拟阵元，128→256
%    排列：[真实1, 虚拟1', 真实2, 虚拟2', ...]，左右交替 jitter 破坏周期性相干
%  双侧模式（bilateral）：每个真实阵元两侧各增加1个虚拟阵元，128→384
%    排列：[左虚1, 真实1, 右虚1, ...]，对称偏移
%  核心：用"角向卷积重采样 H_theta"生成更密通道，
%        再用"离散几何TOF差"做相位搬运（分数延时）
%  frac_shift_lo/hi 分别对应低频保守和高频标准两套几何
% ============================================================
if enableVirtualElements

    virtualInterpKernel = 'zoh';   % 'zoh' | 'linear' | 'cubic'
    keys_a = -0.5;                 % Keys cubic parameter（仅 'cubic' 时使用）

    frac_shift_build_list = [frac_shift_lo, frac_shift_hi];
    geo_ve_built = cell(1, 2);

    for fsi = 1:2
        fs_cur = frac_shift_build_list(fsi);

        N0         = orig_numElements;
        x_base     = transducerPositionsX(:).';
        y_base     = transducerPositionsY(:).';
        x_idx_base = x_idx(:);
        y_idx_base = y_idx(:);

        theta_base  = unwrap(atan2(y_base, x_base));
        dtheta_next = circshift(theta_base,-1) - theta_base;
        dtheta_next(dtheta_next<=0) = dtheta_next(dtheta_next<=0) + 2*pi;
        dtheta_prev = theta_base - circshift(theta_base,1);
        dtheta_prev(dtheta_prev<=0) = dtheta_prev(dtheta_prev<=0) + 2*pi;
        r_base = sqrt(x_base.^2 + y_base.^2);

        switch lower(virtualElementsMode)
            case 'unilateral'
                % ---- 单侧交替：128→256，[真实, 虚拟] 交替排列，sgn 左右交替 ----
                upsample_factor = 2;
                sgn = ones(size(theta_base)); sgn(2:2:end) = -1;
                dtheta_use = dtheta_next; dtheta_use(sgn<0) = dtheta_prev(sgn<0);
                theta_virtual = theta_base + sgn.*fs_cur.*dtheta_use;
                x_virtual = r_base.*cos(theta_virtual);
                y_virtual = r_base.*sin(theta_virtual);
                % 交错插入：[真实1, 虚拟1', 真实2, 虚拟2', ...]
                tx_ve = reshape([x_base; x_virtual],1,[]);
                ty_ve = reshape([y_base; y_virtual],1,[]);

            case 'bilateral'
                % ---- 双侧对称：128→384，[左虚, 真实, 右虚] 排列 ----
                upsample_factor = 3;
                theta_left  = theta_base - fs_cur.*dtheta_prev;
                theta_right = theta_base + fs_cur.*dtheta_next;
                x_left  = r_base.*cos(theta_left);  y_left  = r_base.*sin(theta_left);
                x_right = r_base.*cos(theta_right); y_right = r_base.*sin(theta_right);
                % 交错插入：[左虚1, 真实1, 右虚1, 左虚2, 真实2, 右虚2, ...]
                tx_ve = reshape([x_left; x_base; x_right],1,[]);
                ty_ve = reshape([y_left; y_base; y_right],1,[]);

            otherwise
                error('Unknown virtualElementsMode. Use ''unilateral'' or ''bilateral''.');
        end

        L_ve = upsample_factor;
        N_ve = N0;
        numElements_ve = L_ve * N_ve;
        base_map = repelem(1:N_ve, L_ve);   % [1 1 2 2 ...] 或 [1 1 1 2 2 2 ...]

        x_idx_ve = dsearchn(xi(:), tx_ve(:));
        y_idx_ve = dsearchn(yi(:), ty_ve(:));
        ind_ve   = sub2ind([Nyi, Nxi], y_idx_ve, x_idx_ve);

        % 构造 H_theta（L_ve*N × N）
        % 行排列（unilateral）：2n-1=真实，2n=虚拟
        % 行排列（bilateral）：3n-2=左虚，3n-1=真实，3n=右虚
        H_theta = zeros(L_ve*N_ve, N_ve, 'single');
        for n = 1:N_ve
            switch lower(virtualElementsMode)
                case 'unilateral'
                    H_theta(2*n-1, n) = 1;  % real sample
                    switch lower(virtualInterpKernel)
                        case 'zoh'
                            H_theta(2*n, n) = 1;
                        case 'linear'
                            if sgn(n) > 0
                                nb = n+1; if nb>N_ve, nb=1; end
                            else
                                nb = n-1; if nb<1, nb=N_ve; end
                            end
                            H_theta(2*n, n)  = 1-fs_cur;
                            H_theta(2*n, nb) = fs_cur;
                        case 'cubic'
                            if sgn(n) > 0
                                idx_v = mod([n-1,n,n+1,n+2]-1,N_ve)+1;
                                w_v   = keysCubicWeights(fs_cur, keys_a);
                            else
                                idx_v = mod([n-2,n-1,n,n+1]-1,N_ve)+1;
                                w_v   = keysCubicWeights(1-fs_cur, keys_a);
                            end
                            w_v = w_v/(sum(w_v)+eps);
                            for kk=1:4, H_theta(2*n, idx_v(kk)) = w_v(kk); end
                        otherwise
                            error('Unknown virtualInterpKernel.');
                    end

                case 'bilateral'
                    H_theta(3*n-1, n) = 1;  % real sample
                    nb_prev = n-1; if nb_prev<1, nb_prev=N_ve; end
                    nb_next = n+1; if nb_next>N_ve, nb_next=1; end
                    switch lower(virtualInterpKernel)
                        case 'zoh'
                            H_theta(3*n-2, n) = 1;
                            H_theta(3*n,   n) = 1;
                        case 'linear'
                            H_theta(3*n-2, n)       = 1-fs_cur;
                            H_theta(3*n-2, nb_prev) = fs_cur;
                            H_theta(3*n,   n)       = 1-fs_cur;
                            H_theta(3*n,   nb_next) = fs_cur;
                        case 'cubic'
                            idx_L = mod([n-2,n-1,n,n+1]-1,N_ve)+1;
                            w_L   = keysCubicWeights(1-fs_cur,keys_a);
                            w_L   = w_L/(sum(w_L)+eps);
                            for kk=1:4, H_theta(3*n-2, idx_L(kk)) = w_L(kk); end
                            idx_R = mod([n-1,n,n+1,n+2]-1,N_ve)+1;
                            w_R   = keysCubicWeights(fs_cur,keys_a);
                            w_R   = w_R/(sum(w_R)+eps);
                            for kk=1:4, H_theta(3*n,   idx_R(kk)) = w_R(kk); end
                        otherwise
                            error('Unknown virtualInterpKernel.');
                    end
            end
        end

        % 用 H_theta 生成扩展后的 REC_DATA_ve_full（角向卷积重采样）
        REC_DATA_ve_full = zeros(L_ve*N_ve, L_ve*N_ve, numel(fDATA), 'like', REC_DATA);
        for f_idx = 1:numel(fDATA)
            D0 = REC_DATA(:,:,f_idx);
            REC_DATA_ve_full(:,:,f_idx) = H_theta * D0 * H_theta.';
        end

        % 离散几何TOF差相位搬运
        % 注意：是在"插值后整体"做相位搬运（锚点近似）
        % 更严格的做法是对每个贡献项分别搬运再加权，但复杂度会明显上升
        enableCopyPhaseTransport = true;
        if enableCopyPhaseTransport
            x_base_disc = xi(x_idx_base); y_base_disc = yi(y_idx_base);
            [Xb_tx, Xb_rx] = meshgrid(x_base_disc, x_base_disc);
            [Yb_tx, Yb_rx] = meshgrid(y_base_disc, y_base_disc);
            TOF_base_disc        = sqrt((Xb_tx-Xb_rx).^2+(Yb_tx-Yb_rx).^2)/c_geom;
            TOF_base_disc_mapped = TOF_base_disc(base_map, base_map);

            x_new_disc = xi(x_idx_ve); y_new_disc = yi(y_idx_ve);
            [Xn_tx, Xn_rx] = meshgrid(x_new_disc, x_new_disc);
            [Yn_tx, Yn_rx] = meshgrid(y_new_disc, y_new_disc);
            TOF_new_disc = sqrt((Xn_tx-Xn_rx).^2+(Yn_tx-Yn_rx).^2)/c_geom;
            dTOF = single(TOF_new_disc - TOF_base_disc_mapped);

            for f_idx = 1:numel(fDATA)
                REC_DATA_ve_full(:,:,f_idx) = REC_DATA_ve_full(:,:,f_idx) .* ...
                    exp(1i*sign_conv*2*pi*fDATA(f_idx)*dTOF);
            end
        end

        % 真实Rx掩码：unilateral→奇数位（1,3,...）；bilateral→2:3:end（2,5,8,...）
        isRealRx_ve = false(1, numElements_ve);
        switch lower(virtualElementsMode)
            case 'unilateral'; isRealRx_ve(1:L_ve:end) = true;
            case 'bilateral';  isRealRx_ve(2:L_ve:end) = true;
        end

        % tx_include（仅真实Tx 或 全部Tx，取决于 useVirtualTx）
        if useVirtualTx
            tx_include_ve = 1:dwnsmp:numElements_ve;                   %#ok<UNRCH>
        else
            switch lower(virtualElementsMode)
                case 'unilateral'; tx_include_ve = 1:L_ve:numElements_ve; % 1,3,5,...,255
                case 'bilateral';  tx_include_ve = 2:L_ve:numElements_ve; % 2,5,8,...,383
            end
            tx_include_ve = tx_include_ve(1:dwnsmp:end);
        end

        % Near-end exclusion（VE版本）
        theta_all_ve   = atan2(ty_ve(:), tx_ve(:));
        elemInclude_ve = true(numElements_ve, numElements_ve);
        numELRE_ve     = round(numElements_ve * frac_excl / 2);
        switch lower(exclMode)
            case 'angle'
                for tx_element = 1:numElements_ve
                    dtheta = theta_all_ve - theta_all_ve(tx_element);
                    dtheta = atan2(sin(dtheta), cos(dtheta));
                    elemInclude_ve(tx_element, abs(dtheta)<=theta_excl) = false;
                end
            case 'frac'
                elemExcl_ve = -numELRE_ve:numELRE_ve;
                for tx_element = 1:numElements_ve
                    ec = elemExcl_ve + tx_element;
                    ec(ec<1) = numElements_ve + ec(ec<1);
                    ec(ec>numElements_ve) = ec(ec>numElements_ve) - numElements_ve;
                    elemInclude_ve(tx_element, ec) = false;
                end
        end

        % Remove Outliers（VE套装）
        REC_DATA_ve_full(isnan(REC_DATA_ve_full)) = 0;
        for f_idx = 1:numel(fDATA)
            slice = REC_DATA_ve_full(tx_include_ve,:,f_idx);
            mag   = elemInclude_ve(tx_include_ve,:) .* abs(slice);
            n_out = ceil((1-perc_outliers)*numel(mag));
            [~,idx_out] = maxk(mag(:), n_out);
            slice(idx_out) = 0;
            REC_DATA_ve_full(tx_include_ve,:,f_idx) = slice;
        end

        % 打包本套装
        g_tmp.numElements          = numElements_ve;
        g_tmp.transducerPositionsX = tx_ve;
        g_tmp.transducerPositionsY = ty_ve;
        g_tmp.x_idx                = x_idx_ve;
        g_tmp.y_idx                = y_idx_ve;
        g_tmp.ind                  = ind_ve;
        g_tmp.elemInclude          = elemInclude_ve;
        g_tmp.isRealRx             = isRealRx_ve;
        g_tmp.tx_include           = tx_include_ve;
        g_tmp.REC_DATA             = REC_DATA_ve_full(tx_include_ve,:,:);
        g_tmp.frac_shift           = fs_cur;  % 记录本套装使用的偏移量

        geo_ve_built{fsi} = g_tmp;
        fprintf('  [VE预计算] frac_shift=%.2f 完成（套装 %d/2）\n', fs_cur, fsi);
    end

    % 命名两套（低频保守 / 高频标准）
    geo_ve_lo = geo_ve_built{1};   % frac_shift = frac_shift_lo
    geo_ve_hi = geo_ve_built{2};   % frac_shift = frac_shift_hi

end

%% ===================== 构造 VE_flags 向量 =====================
%  VE_flags(f_idx) = true  → 该频率使用虚拟阵元（低频→geo_ve_lo，高频→geo_ve_hi）
%  VE_flags(f_idx) = false → 该频率使用 geo_orig
VE_flags = false(1, numel(fDATA));
if enableVirtualElements
    switch lower(VE_freqMode)
        case 'always'
            VE_flags(:) = true;

        case 'threshold'
            VE_flags = (fDATA >= VE_freqThresh);

        case 'manual'
            if ~exist('VE_manualFlags','var') || isempty(VE_manualFlags)
                error('VE_freqMode=manual 时必须定义 VE_manualFlags，长度=%d', numel(fDATA));
            end
            if numel(VE_manualFlags) ~= numel(fDATA)
                error('VE_manualFlags 长度(%d)与 fDATA 长度(%d)不匹配', ...
                      numel(VE_manualFlags), numel(fDATA));
            end
            VE_flags = logical(VE_manualFlags(:).');
        otherwise
            error('未知 VE_freqMode');
    end
end

% 打印频率-VE分配预览（含 frac_shift 信息）
fprintf('\n=== VE 频率分配预览 ===\n');
fprintf('模式: %s | frac_shift: %.2f(<%.2fMHz) / %.2f(>=%.2fMHz)\n', ...
    VE_freqMode, frac_shift_lo, frac_shift_f_split/1e6, ...
    frac_shift_hi, frac_shift_f_split/1e6);
fprintf('%-10s  %-10s  %-6s  %-10s\n','freq(MHz)','type','VE','frac_shift');
for fi = 1:numel(fDATA)
    if niterAttenPerFreq(fi) == 0, it = 'SoS'; else, it = 'SoS+Att'; end
    if VE_flags(fi) && enableVirtualElements
        if fDATA(fi) < frac_shift_f_split, fs_disp = frac_shift_lo;
        else,                              fs_disp = frac_shift_hi; end
    else
        fs_disp = 0;
    end
    fprintf('  %-8.3f  %-10s  %d     %.2f\n', fDATA(fi)/1e6, it, VE_flags(fi), fs_disp);
end
fprintf('======================\n\n');

%% ===================== Initialize Model =====================
% Initial Constant Sound Speed Map [m/s]
c_init = 1480; % Initial Homogeneous Sound Speed [m/s] Guess
VEL_INIT   = c_init*ones(Nyi,Nxi);

% Initial Constant Attenuation [Np/(Hz m)]
ATTEN_INIT = 0*alphaNp*ones(Nyi,Nxi);

% (Nonlinear) Conjugate Gradient
search_dir = zeros(Nyi,Nxi); % Conjugate Gradient Direction
gradient_img_prev = zeros(Nyi,Nxi); % Previous Gradient Image;

VEL_ESTIM   = VEL_INIT;
ATTEN_ESTIM = ATTEN_INIT;
% Initial Slowness Image [s/m]
SLOW_ESTIM  = 1./VEL_ESTIM + 1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi);

c0 = mean(VEL_ESTIM(:)); cutoff = 0.75; ord = Inf; % Parameters for Ringing Removal Filter
% Values to Save at Each Iteration
Niter = sum(niterSoSPerFreq) + sum(niterAttenPerFreq);

% 背景环带（用于环伪影比 R_k 统计）：避开中心结构与最外层边界/PML
Rmax = max(R_grid(:));
ring_bg_inner_ratio = 0.65;
ring_bg_outer_ratio = 0.92;
ring_bg_mask = (R_grid >= ring_bg_inner_ratio*Rmax) & (R_grid <= ring_bg_outer_ratio*Rmax);
% 外环保护软掩膜（0~1）：中心0，外圈逐步增加到1
edge_guard_mask = (R_grid - edge_r_inner_ratio*Rmax) ./ ...
                  max((edge_r_outer_ratio - edge_r_inner_ratio)*Rmax, eps);
edge_guard_mask = min(max(edge_guard_mask, 0), 1);
edge_outer_mask = edge_guard_mask > 0.95;  % 仅统计最外环污染

% 粗网格迭代结果
VEL_ESTIM_ITER   = zeros(Nyi, Nxi, Niter);
ATTEN_ESTIM_ITER = zeros(Nyi, Nxi, Niter);
GRAD_IMG_ITER    = zeros(Nyi, Nxi, Niter);
SEARCH_DIR_ITER  = zeros(Nyi, Nxi, Niter);
PHI_ITER         = nan(1, Niter);
DELTA_PHI_ITER   = nan(1, Niter);
TAU_SHIFT_ITER   = nan(1, Niter);
R_RING_ITER      = nan(1, Niter);
GAMMA_TRUST_ITER = nan(1, Niter);
ALPHA_GATE_ITER  = ones(1, Niter);
LAMBDA_POLAR_EFF_ITER = nan(1, Niter);
LAMBDA_EDGE_EFF_ITER  = nan(1, Niter);
EDGE_DAMP_EFF_ITER    = nan(1, Niter);
EDGE_CONTAM_ITER      = nan(1, Niter);
EDGE_BLEND_EFF_ITER   = nan(1, Niter);

%% ---------- 用户控制：迭代执行选项 ----------
runAllIterations    = false;
requestedIterations = 60;
if requestedIterations >= Niter, runAllIterations = true; end
savedIters   = [];
stopNow      = false;
prev_VE_flag = false;   % 用于检测跨频率套装切换
% [修改] 三阶段LF先验参考模型保存状态标志（方向E）
% stage1_ref_saved：阶段1→2边界是否已保存参考（阶段2使用此参考约束）
% stage2_ref_saved：阶段2→3边界是否已保存参考（阶段3使用此参考约束）
stage1_ref_saved = false;
stage2_ref_saved = false;
prev_phi_k       = nan;
prev_tau_k       = nan;
prev_ring_k      = nan;
prev_pred_drop   = nan;
prev_edge_contam = nan;
% -----------------------------------------------

%% ---------- 路径 & GIF 初始化（主循环前统一定义）----------
result_dir = 'D:\Document_ING_fws\WaveformInversionUST\Results\start20260303\';
if ~exist(result_dir,'dir'), mkdir(result_dir); end

if saveGIF
    gif_filepath = [result_dir, filename, '_', verTag, '_NOPP_VEL_anim.gif'];
    fprintf('[GIF] 将保存至: %s\n', gif_filepath);
    gif_initialized = false;   % 标记：首帧是否已写入
    gif_iter_count  = 0;       % 帧计数器

    hgif = figure(99);
    set(hgif, 'Visible','on', 'Position',[100 100 640 560]);
    axgif = axes('Parent', hgif);
    hgif_im = imagesc(axgif, xi_original, yi_original, VEL_ESTIM, crange);
    axis(axgif, 'image');
    colorbar(axgif);
    colormap(axgif, cmap_rb);
    xlabel(axgif, 'Lateral [m]');
    ylabel(axgif, 'Axial [m]');
end

mainLoopTimer = tic;   % 主反演循环总时长（只统计迭代过程）
% -------------------------------------------------------

%% ===================== Main Inversion Loop =====================
for f_idx = 1:numel(fDATA)
    % Iterations at Each Frequency
    for iter_f_idx = 1:(niterSoSPerFreq(f_idx)+niterAttenPerFreq(f_idx))
        tic;
        % Step 1: Accumulate Backprojection Over Each Element
        iter = iter_f_idx + sum(niterSoSPerFreq(1:f_idx-1)) + sum(niterAttenPerFreq(1:f_idx-1));
 
        if ~runAllIterations && (numel(savedIters) >= requestedIterations)
            iter    = numel(savedIters);
            stopNow = true;
            disp(['Stopped early after ', num2str(numel(savedIters)), ' iters.']);
            break;
        end

        % ---- 按频率动态选择几何套装 ----
        % 低频 → geo_ve_lo（frac_shift_lo，保守插值，降低虚假相位）
        % 高频 → geo_ve_hi（frac_shift_hi，标准插值，提供角向支撑）
        if VE_flags(f_idx) && enableVirtualElements
            if fDATA(f_idx) < frac_shift_f_split
                geo_cur    = geo_ve_lo;
                cur_fs_tag = sprintf('VE_lo(%.2f)', frac_shift_lo);
            else
                geo_cur    = geo_ve_hi;
                cur_fs_tag = sprintf('VE_hi(%.2f)', frac_shift_hi);
            end
        else
            geo_cur    = geo_orig;
            cur_fs_tag = 'orig';
        end
        numElements_cur = geo_cur.numElements;
        tx_include_cur  = geo_cur.tx_include;
        ind_cur         = geo_cur.ind;
        elemInclude_cur = geo_cur.elemInclude;
        isRealRx_cur    = geo_cur.isRealRx;
        REC_DATA_cur    = geo_cur.REC_DATA;   % (numel(tx_include) x numElements x numel(fDATA))
        % ----------------------------------

        % Reset CG at Each Frequency (SoS and Attenuation)
        % 情况1：每个频率 SoS 阶段第一次迭代
        % 情况2：同频率从 SoS 切换到 SoS+Atten 阶段
        % 情况3：跨频率发生 VE 套装切换（几何/数据变化，梯度不可延续）
        do_reset_CG = false;
        if iter_f_idx == 1
            do_reset_CG = true;
            if (f_idx > 1) && (VE_flags(f_idx) ~= prev_VE_flag)
                fprintf('>>> 频率 %.3f MHz: VE状态切换 (%d→%d)，强制重置CG\n', ...
                    fDATA(f_idx)/1e6, prev_VE_flag, VE_flags(f_idx));
            end
        elseif iter_f_idx == 1 + niterSoSPerFreq(f_idx)
            do_reset_CG = true;
        end

        if do_reset_CG
            search_dir        = zeros(Nyi, Nxi);
            gradient_img_prev = zeros(Nyi, Nxi);
            prev_phi_k        = nan;
            prev_tau_k        = nan;
            prev_ring_k       = nan;
            prev_pred_drop    = nan;
            prev_edge_contam  = nan;
        end

        if iter_f_idx == 1
            prev_VE_flag = VE_flags(f_idx);
        end

        % [修改] ---- 三阶段LF参考模型保存（方向E）----
        % 阶段1→2边界（首次越过 f_stage1_cutoff）：保存 VEL_stage1_ref
        %   时机：iter_f_idx==1 时检测，VEL_ESTIM 为上一频率迭代完成后的状态，
        %         即阶段1末尾的完整重建结果，对应师兄时域方案的"Low段末"。
        % 阶段2→3边界（首次越过 f_stage2_cutoff）：保存 VEL_stage2_ref
        %   时机：同上，对应师兄时域方案的"Mid段末"。
        % 两个标志位确保各只保存一次，不因阶段内多个频率重复覆盖。
        if enableLFPrior && ~stage1_ref_saved && fDATA(f_idx) > f_stage1_cutoff
            VEL_stage1_ref   = VEL_ESTIM;
            stage1_ref_saved = true;
            fprintf('[LF 3-stage] 阶段1→2：已保存参考模型（f_stage1=%.3f MHz，当前 f=%.3f MHz）\n', ...
                f_stage1_cutoff/1e6, fDATA(f_idx)/1e6);
        end
        if enableLFPrior && ~stage2_ref_saved && fDATA(f_idx) > f_stage2_cutoff
            VEL_stage2_ref   = VEL_ESTIM;
            stage2_ref_saved = true;
            fprintf('[LF 3-stage] 阶段2→3：已保存参考模型（f_stage2=%.3f MHz，当前 f=%.3f MHz）\n', ...
                f_stage2_cutoff/1e6, fDATA(f_idx)/1e6);
        end
        % -----------------------------------------------------------------------

        updateAttenuation = (iter_f_idx > niterSoSPerFreq(f_idx));
        gradient_img = zeros(Nyi, Nxi);

        % Sources (only for selected tx_include_cur)
        SRC = zeros(Nyi, Nxi, numel(tx_include_cur));
        for elmt_idx = 1:numel(tx_include_cur)
            x_idx_src = geo_cur.x_idx(tx_include_cur(elmt_idx));
            y_idx_src = geo_cur.y_idx(tx_include_cur(elmt_idx));
            SRC(y_idx_src, x_idx_src, elmt_idx) = 1;
        end

        % IMPORTANT: Helmholtz uses coarse grid
        HS = HelmholtzSolver(xi_original, yi_original, VEL_ESTIM, ATTEN_ESTIM, ...
            fDATA(f_idx), sign_conv, a0, L_PML);
        [WVFIELD, VIRT_SRC] = HS.solve(SRC, false);

        if updateAttenuation
            VIRT_SRC = 1i*sign(sign_conv)*VIRT_SRC;
        end

        % fws: HV_polar 当前 pass 使用的相位权重
        % SoS 阶段（pass 1）：声速主要由相位约束 → 高 alpha_hv，增强相位分量
        % Atten 阶段（pass 2）：衰减主要由幅度约束 → 低 alpha_hv_atten，增强幅度分量
        % misfitType='L2' 时此变量不参与计算，不影响原有流程
        if updateAttenuation
            alpha_hv_cur = alpha_hv_atten;
        else
            alpha_hv_cur = alpha_hv;
        end

        % Build Adjoint Sources
        scaling = zeros(numel(tx_include_cur), 1);
        ADJ_SRC = zeros(Nyi, Nxi, numel(tx_include_cur));
        phi_k_accum = 0;
        tau_k_accum = 0;
        tau_k_count = 0;

        for elmt_idx = 1:numel(tx_include_cur)
            WVFIELD_elmt = WVFIELD(:,:,elmt_idx);
            tx_id = tx_include_cur(elmt_idx);

            rx_mask_all = elemInclude_cur(tx_id, :);
            idx_all     = ind_cur(rx_mask_all);
            REC_SIM_all = WVFIELD_elmt(idx_all);
            REC_all     = REC_DATA_cur(elmt_idx, rx_mask_all, f_idx);
            REC_all     = REC_all(:);

            % fws: scaling computed only on REAL receivers to avoid virtual-Rx bias
            rx_mask_scale = rx_mask_all & isRealRx_cur;
            idx_scale     = ind_cur(rx_mask_scale);
            REC_SIM_scale = WVFIELD_elmt(idx_scale);
            REC_scale     = REC_DATA_cur(elmt_idx, rx_mask_scale, f_idx);
            REC_scale     = REC_scale(:);

            scaling(elmt_idx) = (REC_SIM_scale(:)'*REC_scale(:)) / ...
                                (REC_SIM_scale(:)'*REC_SIM_scale(:));

            % fws: 计算失配残差（支持 L2 与 HV_polar 两种模式）
            % L2：标准复数残差，逐样本比较实部/虚部
            % HV_polar：极坐标分解
            %   p_sim_sc = γ·p_sim：经过源标定的模拟 phasor（复数列向量）
            %   phi_dir  = p_sim_sc / |p_sim_sc|：单位相量（幅度方向基矢）
            %   ΔA       = |p_sim_sc| - |p_obs|：幅度残差（实数）
            %              → 对应 L2 的径向分量，约束衰减
            %   Δφ       = angle(p_sim_sc · conj(p_obs))：包裹相位差 ∈ (-π,π]
            %              → 对应 L2 的切向分量，约束声速，凸性优于 L2
            %   伴随源  = (1-α)·ΔA·phi_dir + α·Δφ·(i·phi_dir)
            %              其中 i·phi_dir 是 phi_dir 的正交切向方向
            %   当 α=0 且幅度误差主导时退化为 L2；
            %   当 α=1 时为纯相位失配，对周期跳跃最鲁棒。
            p_sim_sc = scaling(elmt_idx) * REC_SIM_all;
            switch lower(misfitType)
                case 'l2'
                    residual = p_sim_sc - REC_all;
                case 'hv_polar'
                    amp_sim  = abs(p_sim_sc) + eps;           % |p_sim|，加 eps 防零除
                    phi_dir  = p_sim_sc ./ amp_sim;           % exp(i·φ_sim)，单位相量
                    amp_res  = amp_sim - abs(REC_all);         % ΔA [实数]
                    dphi     = angle(p_sim_sc .* conj(REC_all)); % Δφ ∈ (-π,π] [实数]
                    residual = (1 - alpha_hv_cur) .* amp_res .* phi_dir + ...
                                alpha_hv_cur      .* dphi    .* (1i * phi_dir);
                otherwise
                    error('Unknown misfitType: %s. Use ''L2'' or ''HV_polar''.', misfitType);
            end
            phi_k_accum = phi_k_accum + sum(abs(residual).^2);
            % 走时残差近似：单频相位差换算时移 Δt = |Δφ|/(2πf)
            dphi_tau = angle(p_sim_sc .* conj(REC_all));
            tau_k_accum = tau_k_accum + sum(abs(dphi_tau) / (2*pi*fDATA(f_idx)));
            tau_k_count = tau_k_count + numel(dphi_tau);

            ADJ_SRC_elmt = zeros(Nyi, Nxi);

            % fws: VE通道置信度加权（方向A2）
            % 真实Rx权重=1，虚拟Rx权重=beta_ve，降低插值偏差对伴随源注入的污染
            % geo_orig套装（isRealRx全1）下此操作等价于乘以全1，无额外开销
            w_channel = ones(sum(rx_mask_all), 1);
            w_channel(~isRealRx_cur(rx_mask_all)) = beta_ve;
            ADJ_SRC_elmt(idx_all) = w_channel .* residual;

            ADJ_SRC(:,:,elmt_idx) = ADJ_SRC_elmt;
        end

        % ---------- 可信度监控量（数据一致性 + 环伪影趋势） ----------
        phi_k = real(phi_k_accum) / max(tau_k_count,1);
        tau_k = tau_k_accum / max(tau_k_count,1);

        if ~isnan(prev_phi_k)
            delta_phi_k = prev_phi_k - phi_k;
        else
            delta_phi_k = nan;
        end

        if ~isnan(delta_phi_k) && ~isnan(prev_pred_drop) && abs(prev_pred_drop) > eps
            gamma_k = delta_phi_k / prev_pred_drop;  % 实际下降 / 预测下降
        else
            gamma_k = nan;
        end

        % Backproject error
        ADJ_WVFIELD = HS.solve(ADJ_SRC, true);
        SCALING     = repmat(reshape(scaling,[1,1,numel(scaling)]), [Nyi,Nxi,1]);
        BACKPROJ    = -real(conj(SCALING.*VIRT_SRC).*ADJ_WVFIELD);

        for elmt_idx = 1:numel(tx_include_cur)
            gradient_img = gradient_img + BACKPROJ(:,:,elmt_idx);
        end

        % fws: 照明补偿预条件（方向D）
        % 用各发射源正向波场强度之和估计伪Hessian对角，均衡稀疏阵角向不均匀照明。
        % 对梯度做归一化预条件后再做Ringing去除，使环状照明偏差在进入CG前已被修正。
        % enableIllumComp=false 时完全跳过，退回标准梯度，方便消融对比。
        if enableIllumComp
            illum = zeros(Nyi, Nxi);
            for elmt_idx = 1:numel(tx_include_cur)
                illum = illum + abs(WVFIELD(:,:,elmt_idx)).^2;
            end
            % 归一化到 [0,1]，用 eps_illum 截断防止低照明区除零放大噪声
            illum_norm   = illum / (max(illum(:)) + eps);
            gradient_img = gradient_img ./ (illum_norm + eps_illum);
        end

        % [修改] ---- 三阶段低频先验正则化梯度项（方向E）----
        % 根据当前频率所处阶段，选择对应的参考模型施加 Tikhonov 约束：
        %   阶段2（f_stage1 < f ≤ f_stage2）：约束拉向 VEL_stage1_ref
        %   阶段3（f > f_stage2）           ：约束拉向 VEL_stage2_ref
        % 各阶段内：lambda_k = lambda_stage · (f_stageN_cutoff / fDATA(f_idx))²
        %   → 阶段内频率越高，约束越弱（符合"低频定背景、高频补细节"的物理直觉）
        % 相对归一化：消除速度偏差[m/s]与gradient_img量纲不匹配的问题。
        %   归一化公式：lf_prior_reg = -lambda_k · max(|grad|) · (lf_diff / max(|lf_diff|))
        %   lambda_stage 含义：LF约束贡献占当前数据梯度最大幅值的比例（无量纲）
        % 注：仅在 SoS 更新阶段施加（衰减阶段不对声速先验施加额外约束）
        if enableLFPrior && ~updateAttenuation
            lf_ref_cur   = [];         % 当前阶段使用的参考模型
            f_cutoff_cur = [];         % 当前阶段的分界频率（用于 lambda_k 计算）
            if stage2_ref_saved && fDATA(f_idx) > f_stage2_cutoff
                % 阶段3：使用阶段2末尾参考
                lf_ref_cur   = VEL_stage2_ref;
                f_cutoff_cur = f_stage2_cutoff;
            elseif stage1_ref_saved && fDATA(f_idx) > f_stage1_cutoff
                % 阶段2：使用阶段1末尾参考
                lf_ref_cur   = VEL_stage1_ref;
                f_cutoff_cur = f_stage1_cutoff;
            end
            % 阶段1（f <= f_stage1_cutoff）：lf_ref_cur 为空，不施加约束

            if ~isempty(lf_ref_cur)
                lambda_k      = lambda_stage * (f_cutoff_cur / fDATA(f_idx))^2;
                lf_diff       = VEL_ESTIM - lf_ref_cur;               % 速度偏差 [m/s]
                grad_scale_lf = max(abs(gradient_img(:))) + eps;       % 数据梯度幅值锚点
                lf_diff_scale = max(abs(lf_diff(:))) + eps;            % 速度差幅值锚点
                lf_prior_reg  = -lambda_k * grad_scale_lf * (lf_diff / lf_diff_scale);
                gradient_img  = gradient_img + lf_prior_reg;
            end
        end
        % -----------------------------------------------------------------------

        % [新增] ---- 各向同性TV正则化梯度项（方向F）----
        % 对当前声速图 VEL_ESTIM 计算 Huber-TV 的负散度，叠加到 gradient_img：
        %   ∂TV/∂v = -div(∇v / |∇v|_Huber)
        % 相对归一化：消除 div(∇v/|∇v|) 量纲 [1/m] 与 gradient_img 不匹配的问题。
        %   归一化公式：lambda_tv · max(|grad|) · (-div_term / max(|div_term|))
        %   lambda_tv 含义：TV贡献占数据梯度最大幅值的比例（此版本设为 2e-3，即 2%）
        % enableTV=false 时完全跳过，退回标准梯度，方便消融对比。
        if enableTV
            [gx_tv, gy_tv] = gradient(VEL_ESTIM, dxi_original, dxi_original);
            norm_tv        = sqrt(gx_tv.^2 + gy_tv.^2 + eps_tv^2);   % Huber 平滑范数
            px_tv          = gx_tv ./ norm_tv;                        % 归一化梯度 x 分量
            py_tv          = gy_tv ./ norm_tv;                        % 归一化梯度 y 分量
            [dpx_x, ~]     = gradient(px_tv, dxi_original, dxi_original);
            [~, dpy_y]     = gradient(py_tv, dxi_original, dxi_original);
            div_term_tv    = dpx_x + dpy_y;                           % div(∇v/|∇v|_Huber) [1/m]
            grad_scale_tv  = max(abs(gradient_img(:))) + eps;         % 数据梯度幅值锚点（含LF prior后）
            div_tv_scale   = max(abs(div_term_tv(:))) + eps;          % TV散度幅值锚点
            gradient_img   = gradient_img + lambda_tv * grad_scale_tv * (-div_term_tv / div_tv_scale);
        end
        % -----------------------------------------------------------------------

        % [新增] ---- 极坐标径向TV正则化梯度项（方向G）----
        % 仅对径向梯度分量 ∂v/∂r = ∇v · er 施加 Huber-TV 约束，
        % 专门抑制同心环状伪影，对切向方向（真实边界）影响更小。
        % 相对归一化：同方向F，消除 div_r 量纲 [1/m] 与 gradient_img 不匹配的问题。
        %   lambda_polar 此版本设为 2e-3（与方向F同幅度，协同抑制同心环）。
        % enablePolarReg=false 时完全跳过，方便消融对比。
        grad_radial = [];
        ring_k      = nan;
        if enablePolarReg
            % 计算极坐标径向单位向量（Xi, Yi 来自成像网格 meshgrid，全程有效）
            R_pol  = sqrt(Xi.^2 + Yi.^2) + eps;   % 到中心距离，加 eps 防零除
            er_x   = Xi ./ R_pol;                 % 径向单位向量 x 分量
            er_y   = Yi ./ R_pol;                 % 径向单位向量 y 分量
            % 径向梯度分量（用方向F已算好的 gx_tv/gy_tv 复用，避免重复 gradient 调用）
            if ~enableTV
                % 若方向F未开启，需在此补算空间梯度
                [gx_tv, gy_tv] = gradient(VEL_ESTIM, dxi_original, dxi_original);
            end
            grad_radial    = gx_tv .* er_x + gy_tv .* er_y;          % ∂v/∂r [标量]
            norm_rad       = abs(grad_radial) + eps_polar_reg;        % Huber 平滑
            % 构造归一化径向梯度矢量场
            px_pol         = (grad_radial ./ norm_rad) .* er_x;
            py_pol         = (grad_radial ./ norm_rad) .* er_y;
            [dpx_pol_x, ~] = gradient(px_pol, dxi_original, dxi_original);
            [~, dpy_pol_y] = gradient(py_pol, dxi_original, dxi_original);
            div_r          = dpx_pol_x + dpy_pol_y;                   % 径向TV散度 [1/m]
            grad_mag       = sqrt(gx_tv.^2 + gy_tv.^2);
            num_ring       = mean(abs(grad_radial(ring_bg_mask)), 'omitnan');
            den_ring       = mean(grad_mag(ring_bg_mask), 'omitnan') + eps;
            ring_k         = num_ring / den_ring;
            grad_scale_pol = max(abs(gradient_img(:))) + eps;         % 数据梯度幅值锚点（含前序正则化）
            div_r_scale    = max(abs(div_r(:))) + eps;                 % 径向TV散度幅值锚点
            lambda_polar_eff = lambda_polar;
            if enableTrustGate && ~isnan(delta_phi_k) && ~isnan(prev_ring_k)
                small_drop = (delta_phi_k < phi_small_drop_ratio * max(prev_phi_k, eps));
                ring_rise  = ~isnan(ring_k) && (ring_k > ring_rise_ratio * prev_ring_k);
                if small_drop && ring_rise
                    lambda_polar_eff = min(lambda_polar * lambda_polar_gate_max, ...
                        lambda_polar * lambda_polar_gate_boost);
                end
            end
            gradient_img   = gradient_img + lambda_polar_eff * grad_scale_pol * (-div_r / div_r_scale);
        else
            lambda_polar_eff = 0;
        end
        % -----------------------------------------------------------------------

        % [新增] ---- 高频外环锚定（边界保护）----
        % 仅在高频SoS阶段启用：抑制中心结构误差向背景外圈扩散
        lambda_edge_eff = 0;
        edge_contam_k   = nan;
        edge_ref_cur    = [];
        if enableEdgeGuard && ~updateAttenuation && (fDATA(f_idx) >= f_edge_guard_start)
            if stage2_ref_saved
                edge_ref = VEL_stage2_ref;   % 优先使用中频末参考
            elseif stage1_ref_saved
                edge_ref = VEL_stage1_ref;
            else
                edge_ref = VEL_INIT;
            end
            edge_ref_cur = edge_ref;
            edge_diff = (VEL_ESTIM - edge_ref) .* edge_guard_mask;
            edge_scale = max(abs(edge_diff(:))) + eps;
            grad_scale_edge = max(abs(gradient_img(:))) + eps;
            lambda_edge_eff = lambda_edge_anchor * (f_edge_guard_start / fDATA(f_idx))^2;
            edge_anchor_reg = -lambda_edge_eff * grad_scale_edge * (edge_diff / edge_scale);
            gradient_img = gradient_img + edge_anchor_reg;
            edge_contam_k = mean(abs(edge_diff(edge_outer_mask)), 'omitnan');
        end
        % -----------------------------------------------------------------------

        % Remove Ringing from Gradient Image
        gradient_img = ringingRemovalFilt(xi_original, yi_original, ...
            gradient_img, c0, fDATA(f_idx), cutoff, ord);

        % % Step 2: Compute New Conjugate Gradient Search Direction from Gradient
        % % Conjugate Gradient Direction Scaling Factor for Updates
        % if do_reset_CG
        %     beta = 0;
        % else
        %     betaPR = (gradient_img(:)'*(gradient_img(:)-gradient_img_prev(:))) / ...
        %              (gradient_img_prev(:)'*gradient_img_prev(:));
        %     betaFR = (gradient_img(:)'*gradient_img(:)) / ...
        %              (gradient_img_prev(:)'*gradient_img_prev(:));
        %     beta = min(max(betaPR,0), betaFR);
        % end
        % search_dir        = beta*search_dir - gradient_img;
        % gradient_img_prev = gradient_img;

        % Step 2: Compute New Conjugate Gradient Search Direction from Gradient
        % fws: HV_polar 调试阶段临时禁用 CG 动量，避免"旧L2动量 + 新HV梯度"叠加失控
        if strcmpi(misfitType, 'HV_polar')
            beta = 0;
        elseif do_reset_CG
            beta = 0;
        else
            betaPR = (gradient_img(:)'*(gradient_img(:)-gradient_img_prev(:))) / ...
                     (gradient_img_prev(:)'*gradient_img_prev(:));
            betaFR = (gradient_img(:)'*gradient_img(:)) / ...
                     (gradient_img_prev(:)'*gradient_img_prev(:));
            beta = min(max(betaPR,0), betaFR);
        end
        search_dir        = beta*search_dir - gradient_img;
        gradient_img_prev = gradient_img;

        % Step 3: Compute Forward Projection of Current Search Direction
        % Forward projection of search direction
        PERTURBED_WVFIELD = HS.solve(VIRT_SRC.*search_dir, false);
        dREC_SIM = zeros(numel(tx_include_cur), numElements_cur);

        for elmt_idx = 1:numel(tx_include_cur)
            PERTURBED_WVFIELD_elmt = PERTURBED_WVFIELD(:,:,elmt_idx);
            mask_elmt = elemInclude_cur(tx_include_cur(elmt_idx),:);
            dREC_SIM(elmt_idx, mask_elmt) = -permute( ...
                scaling(elmt_idx) * PERTURBED_WVFIELD_elmt(ind_cur(mask_elmt)), [2,1]);
        end

        % % Step 4: Perform a Linear Approximation of Exact Line Search
        % perc_step_size = 1; % (<1/2) Introduced to Improve Compliance with Strong Wolfe Conditions 
        % % Line search
        % % den = real(dREC_SIM(:)'*dREC_SIM(:));
        % % alpha = -(gradient_img(:)'*search_dir(:)) / (den + eps(den));
        % alpha = -(gradient_img(:)'*search_dir(:)) / (dREC_SIM(:)'*dREC_SIM(:));

        % Step 4: Perform a Linear Approximation of Exact Line Search
        perc_step_size = 1; % (<1/2) Introduced to Improve Compliance with Strong Wolfe Conditions 
        % Line search
        % den = real(dREC_SIM(:)'*dREC_SIM(:));
        % alpha = -(gradient_img(:)'*search_dir(:)) / (den + eps(den));

        % fws: HV_polar 梯度量纲（弧度）与 L2 phasor 扰动量纲不一致，
        %      原 L2 线搜索会给出错误步长（通常高估 2~3 数量级）。
        %      临时改用归一化固定步长调试：将 search_dir 按最大绝对值归一化，
        %      再乘以 alpha_fix，使每步 max|Δv| ≈ 1~5 m/s。
        %      待 HV_polar 收敛行为稳定后，可换回正确的 Gauss-Newton 线搜索。
        %      alpha_fix_sos 经验起点：5e-5（保守），可向上调至 2e-4。
        if strcmpi(misfitType, 'HV_polar')
            alpha_fix_sos   = 8e-6;   % SoS 阶段固定步长（调整目标：max|Δv| ≈ 1~3 m/s）
            alpha_fix_atten = 1e-5;   % Atten 阶段固定步长
            sd_max = max(abs(search_dir(:))) + eps;
            search_dir = search_dir / sd_max;   % 归一化：只让 alpha 控制更新尺度
            if updateAttenuation
                alpha = alpha_fix_atten;
            else
                alpha = alpha_fix_sos;
            end
        else
            alpha = -(gradient_img(:)'*search_dir(:)) / (dREC_SIM(:)'*dREC_SIM(:));
        end

        % 可信度门控：根据 ΔΦ/τ/R 动态缩步与增强径向TV
        alpha_gate = 1.0;
        if enableTrustGate && ~isnan(delta_phi_k)
            bad_fit = (delta_phi_k <= 0);
            tau_rise = ~isnan(prev_tau_k) && (tau_k > tau_rise_ratio * prev_tau_k);
            small_drop = (delta_phi_k < phi_small_drop_ratio * max(prev_phi_k, eps));
            ring_rise = ~isnan(prev_ring_k) && ~isnan(ring_k) && (ring_k > ring_rise_ratio * prev_ring_k);
            if bad_fit || tau_rise
                alpha_gate = min(alpha_gate, alpha_gate_badfit);
            end
            if small_drop && ring_rise
                alpha_gate = min(alpha_gate, alpha_gate_ringgrowth);
            end
        end
        alpha = alpha * alpha_gate;

        % 高频外环步长阻尼：减少外围被中心结构“带跑”
        edge_damp_eff = 0;
        if enableEdgeGuard && ~updateAttenuation && (fDATA(f_idx) >= f_edge_guard_start)
            edge_damp_eff = edge_step_damp_max;
            if ~isnan(edge_contam_k) && ~isnan(prev_edge_contam) && ...
                    (edge_contam_k > edge_contam_rise_ratio * prev_edge_contam)
                % 外环污染继续上升时，额外增强阻尼
                edge_damp_eff = min(0.90, edge_damp_eff * 1.20);
            end
            search_dir = search_dir .* (1 - edge_damp_eff * edge_guard_mask);
        end

        % fws: 调试 —— 打印本轮最大速度更新量，确认步长合理
        % Δs = α·search_dir → Δv ≈ v² · Δs（一阶近似）
        delta_s_max  = perc_step_size * alpha * max(abs(search_dir(:)));
        delta_v_max  = mean(VEL_ESTIM(:))^2 * delta_s_max;
        fprintf(['  Φ=%.3e | ΔΦ=% .3e | τ=%.3eus | R=%.3f | γ=% .3f' ...
                 ' | EdgeContam=%.3f λpolar=%.2e λedge=%.2e damp=%.2f' ...
                 ' | alpha=%.2e(gate=%.2f) | max|Δv|≈%.3f m/s\n'], ...
            phi_k, delta_phi_k, tau_k*1e6, ring_k, gamma_k, edge_contam_k, ...
            lambda_polar_eff, lambda_edge_eff, edge_damp_eff, alpha, alpha_gate, delta_v_max);

        % 本步的“预测下降”供下一步计算 γ_{k+1}
        pred_drop_cur = -perc_step_size * alpha * real(gradient_img(:)'*search_dir(:));
        
        % Update slowness
        if updateAttenuation
            SI = sign(sign_conv) * imag(SLOW_ESTIM) + perc_step_size * alpha * search_dir;
            SLOW_ESTIM = real(SLOW_ESTIM) + 1i * sign(sign_conv) * SI;
        else
            SLOW_ESTIM = SLOW_ESTIM + perc_step_size * alpha * search_dir;
        end
        VEL_ESTIM   = 1./real(SLOW_ESTIM);
        % 外环后处理回拉：在更新后直接抑制边界漂移（高频SoS阶段）
        edge_blend_eff = 0;
        if enableEdgeGuard && ~updateAttenuation && (fDATA(f_idx) >= f_edge_guard_start) && ~isempty(edge_ref_cur)
            edge_blend_eff = edge_blend_base * (f_edge_guard_start / fDATA(f_idx))^2;
            if ~isnan(edge_contam_k) && ~isnan(prev_edge_contam) && ...
                    (edge_contam_k > edge_contam_rise_ratio * prev_edge_contam)
                edge_blend_eff = min(edge_blend_max, edge_blend_eff * 1.5);
            end
            VEL_ESTIM = VEL_ESTIM .* (1 - edge_blend_eff * edge_guard_mask) + ...
                        edge_ref_cur .* (edge_blend_eff * edge_guard_mask);
        end
        ATTEN_ESTIM = 2*pi*imag(SLOW_ESTIM)*sign(sign_conv);
        SLOW_ESTIM  = 1./VEL_ESTIM + 1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi);

        % Save coarse iteration
        VEL_ESTIM_ITER(:,:,iter)   = VEL_ESTIM;
        ATTEN_ESTIM_ITER(:,:,iter) = ATTEN_ESTIM;
        GRAD_IMG_ITER(:,:,iter)    = gradient_img;
        SEARCH_DIR_ITER(:,:,iter)  = search_dir;
        PHI_ITER(iter)             = phi_k;
        DELTA_PHI_ITER(iter)       = delta_phi_k;
        TAU_SHIFT_ITER(iter)       = tau_k;
        R_RING_ITER(iter)          = ring_k;
        GAMMA_TRUST_ITER(iter)     = gamma_k;
        ALPHA_GATE_ITER(iter)      = alpha_gate;
        LAMBDA_POLAR_EFF_ITER(iter)= lambda_polar_eff;
        LAMBDA_EDGE_EFF_ITER(iter) = lambda_edge_eff;
        EDGE_DAMP_EFF_ITER(iter)   = edge_damp_eff;
        EDGE_CONTAM_ITER(iter)     = edge_contam_k;
        EDGE_BLEND_EFF_ITER(iter)  = edge_blend_eff;
        savedIters(end+1)          = iter; %#ok<SAGROW>

        % 迭代状态滚动到下一次
        prev_phi_k     = phi_k;
        prev_tau_k     = tau_k;
        prev_ring_k    = ring_k;
        prev_pred_drop = pred_drop_cur;
        prev_edge_contam = edge_contam_k;

        % 迭代状态滚动到下一次
        prev_phi_k     = phi_k;
        prev_tau_k     = tau_k;
        prev_ring_k    = ring_k;
        prev_pred_drop = pred_drop_cur;

        % Visualize coarse 
        % 频率/VE信息显示在 figure 窗口标题栏
        hfig = figure(1);
        set(hfig, 'Name', sprintf('f=%.3fMHz | %s | iter=%d/%d | pass=%s | misfit=%s(α=%.2f)', ...
            fDATA(f_idx)/1e6, cur_fs_tag, iter, Niter, ...
            ternary_str(updateAttenuation, 'SoS+Atten', 'SoS'), misfitType, alpha_hv_cur));

        subplot(2,3,1); imagesc(xi_original, yi_original, VEL_ESTIM, crange); axis image; colorbar; colormap(cmap_rb);
        title(['Estimated Wave Velocity ', num2str(iter)]); xlabel('Lateral [m]'); ylabel('Axial [m]');
        subplot(2,3,4); imagesc(xi_original, yi_original, Np2dB*slow2atten*ATTEN_ESTIM, attenrange); axis image; colorbar; colormap(cmap_rb);
        title(['Estimated Attenuation ', num2str(iter)]); xlabel('Lateral [m]'); ylabel('Axial [m]');
        subplot(2,3,2); imagesc(xi_orig, yi_orig, C, crange); axis image; colorbar; colormap(cmap_rb);
        title('True Wave Velocity'); xlabel('Lateral [m]'); ylabel('Axial [m]');
        subplot(2,3,5); imagesc(xi_orig, yi_orig, atten, attenrange); axis image; colorbar; colormap(cmap_rb);
        title('Attenuation'); xlabel('Lateral [m]'); ylabel('Axial [m]');
        subplot(2,3,3); imagesc(xi_original, yi_original, search_dir); axis image; colorbar; colormap(cmap_rb);
        title(['Search Direction ', num2str(iter)]); xlabel('Lateral [m]'); ylabel('Axial [m]');
        subplot(2,3,6); imagesc(xi_original, yi_original, -gradient_img); axis image; colorbar; colormap(cmap_rb);
        title(['Gradient ', num2str(iter)]); xlabel('Lateral [m]'); ylabel('Axial [m]');
        drawnow;

        % ---- GIF 帧保存 ----
        % 复用预创建的 hgif/axgif/hgif_im，只更新 CData，不重建渲染链
        % gif_save_every 控制抽帧频率，减少文件大小和写盘开销
        if saveGIF
            gif_iter_count = gif_iter_count + 1;
            if mod(gif_iter_count - 1, gif_save_every) == 0
                set(hgif_im, 'CData', VEL_ESTIM);
                title(axgif, sprintf('iter=%d  f=%.3fMHz  %s  pass=%s', ...
                    iter, fDATA(f_idx)/1e6, cur_fs_tag, ...
                    ternary_str(updateAttenuation,'SoS+Atten','SoS')), 'FontSize', 10, 'Interpreter', 'none');
                drawnow limitrate nocallbacks;

                frame = getframe(hgif);
                im    = frame2im(frame);
                [A, map] = rgb2ind(im, 128);
                if ~gif_initialized
                    imwrite(A, map, gif_filepath, 'gif', 'LoopCount', Inf, 'DelayTime', gif_delay);
                    gif_initialized = true;
                else
                    imwrite(A, map, gif_filepath, 'gif', 'WriteMode', 'append', 'DelayTime', gif_delay);
                end
            end
        end
        % -------------------

        disp(['Iteration ', num2str(iter), ...
              ' | f=', num2str(fDATA(f_idx)/1e6,'%.3f'), 'MHz', ...
              ' | ', cur_fs_tag, ...
              ' | Atten=', num2str(updateAttenuation), ...
              ' | misfit=', misfitType, sprintf('(α=%.2f)', alpha_hv_cur)]);
        toc;
    end

    if stopNow, break; end
end

if saveGIF && gif_initialized
    fprintf('[GIF] 保存完成：%s\n', gif_filepath);
end

mainLoopElapsedSec = toc(mainLoopTimer);

fprintf('\n================ 反演计时汇总 ================\n');
fprintf('已完成迭代次数: %d\n', numel(savedIters));
fprintf('主反演循环总时长: %.2f 秒 (%.2f 分钟)\n', mainLoopElapsedSec, mainLoopElapsedSec/60);
fprintf('平均每次迭代: %.2f 秒/iter\n', mainLoopElapsedSec/max(numel(savedIters),1));
fprintf('=============================================\n\n');

%% ===================== Final Coarse Display =====================
figure(4);
subplot(2,2,1); imagesc(xi_original, yi_original, VEL_ESTIM, crange); axis image; colorbar; colormap(cmap_rb);
title(['Estimated Wave Velocity ', num2str(iter)]); xlabel('Lateral [m]'); ylabel('Axial [m]');
subplot(2,2,2); imagesc(xi_original, yi_original, Np2dB*slow2atten*ATTEN_ESTIM, attenrange); axis image; colorbar; colormap(cmap_rb);
title(['Estimated Attenuation ', num2str(iter)]); xlabel('Lateral [m]'); ylabel('Axial [m]');
subplot(2,2,3); imagesc(xi_orig, yi_orig, C, crange); axis image; colorbar; colormap(cmap_rb);
title('True Wave Velocity'); xlabel('Lateral [m]'); ylabel('Axial [m]');
subplot(2,2,4); imagesc(xi_orig, yi_orig, atten, attenrange); axis image; colorbar; colormap(cmap_rb);
title('Attenuation'); xlabel('Lateral [m]'); ylabel('Axial [m]');


%% ===================== Axial Sound Speed Profile (cutting line) =====================
% 参考 Wu et al. (2023) Fig.17(c): 沿垂直切线对比重建声速与真实声速
figure(5); clf;

cut_x_m = 0;   % 切线的 x 坐标 [m]，默认取 x=0 中心线

[~, ix_coarse] = min(abs(xi_original - cut_x_m));
[~, ix_orig]   = min(abs(xi_orig     - cut_x_m));

vel_cut_true = C(:, ix_orig);
vel_cut_nopp = VEL_ESTIM(:, ix_coarse);

plot(vel_cut_true, yi_orig,     'k-',  'LineWidth', 2.0, 'DisplayName', 'True'); hold on;
plot(vel_cut_nopp, yi_original, 'b--', 'LineWidth', 1.5, 'DisplayName', 'NoPP');

set(gca, 'YDir', 'reverse');
xlabel('Sound Speed [m/s]'); ylabel('Axial position [m]');
title(sprintf('Cutting line x=%.1fmm (iter=%d)', cut_x_m*1e3, iter));
legend('Location','best'); grid on;
xlim(crange); ylim([yi_original(1), yi_original(end)]);

scriptElapsedSec = toc(scriptTimer);
fprintf('整个脚本总时长: %.2f 秒 (%.2f 分钟)\n', scriptElapsedSec, scriptElapsedSec/60);

%% ===================== Horizontal Sound Speed Profile (contour line) =====================
% 参考 Long et al. (2023) Fig.6: 在SOS图上给出一条水平切线，
% 对比该切线上的重建声速曲线与真实声速曲线。
% 左图：重建图像 + 水平切线
% 右图：沿该切线的 SOS profile

figure(6); clf;

% ------------------ 用户可调参数 ------------------
cut_y_m = 0.01;          % 水平切线的 y 坐标 [m]，默认取 y = 0 中心线
showStageRefs = true; % 若存在阶段参考模型，则一并画出 profile
lineColor = [1 1 0];  % 切线颜色（黄色）
% -------------------------------------------------

% ---- 在粗网格/真值网格上分别找最接近 cut_y_m 的行 ----
[~, iy_coarse] = min(abs(yi_original - cut_y_m));
[~, iy_true]   = min(abs(yi_orig     - cut_y_m));

% ---- 提取 profile ----
vel_prof_est  = VEL_ESTIM(iy_coarse, :);
vel_prof_true = C(iy_true, :);

% 若存在阶段参考模型，可选叠加
has_stage1 = exist('VEL_stage1_ref', 'var') == 1;
has_stage2 = exist('VEL_stage2_ref', 'var') == 1;

if has_stage1
    vel_prof_stage1 = VEL_stage1_ref(iy_coarse, :);
end
if has_stage2
    vel_prof_stage2 = VEL_stage2_ref(iy_coarse, :);
end

% ------------------ 左图：图像 + 水平切线 ------------------
subplot(1,2,1);
imagesc(xi_original, yi_original, VEL_ESTIM, crange);
axis image; colorbar; colormap(cmap_rb);
hold on;
plot([xi_original(1), xi_original(end)], [yi_original(iy_coarse), yi_original(iy_coarse)], ...
    '-', 'Color', lineColor, 'LineWidth', 1.5);
hold off;
xlabel('Lateral [m]');
ylabel('Axial [m]');
title(sprintf('Estimated SOS with horizontal line (y = %.1f mm)', yi_original(iy_coarse)*1e3));

% ------------------ 右图：轮廓曲线 ------------------
subplot(1,2,2);
plot(xi_orig,     vel_prof_true, 'k-',  'LineWidth', 2.0, 'DisplayName', 'True'); hold on;
plot(xi_original, vel_prof_est,  'r--', 'LineWidth', 1.8, 'DisplayName', 'Estimated');

if showStageRefs
    if has_stage1
        plot(xi_original, vel_prof_stage1, 'b-.', 'LineWidth', 1.3, 'DisplayName', 'Stage-1 ref');
    end
    if has_stage2
        plot(xi_original, vel_prof_stage2, 'm-.', 'LineWidth', 1.3, 'DisplayName', 'Stage-2 ref');
    end
end

grid on;
xlabel('Lateral position [m]');
ylabel('Sound Speed [m/s]');
title(sprintf('Horizontal SOS profile at y = %.1f mm', yi_original(iy_coarse)*1e3));
legend('Location','best');
xlim([xi_original(1), xi_original(end)]);
ylim(crange);

%% ===================== Save =====================
suffix = 'WaveformInversionResults';
filename_results = [result_dir, filename, '_', verTag, '_NOPP_', suffix, '.mat'];

if ~isempty(savedIters)
    savedIters = sort(savedIters);
    VEL_ESTIM_ITER   = VEL_ESTIM_ITER(:,:,savedIters);
    ATTEN_ESTIM_ITER = ATTEN_ESTIM_ITER(:,:,savedIters);
    GRAD_IMG_ITER    = GRAD_IMG_ITER(:,:,savedIters);
    SEARCH_DIR_ITER  = SEARCH_DIR_ITER(:,:,savedIters);
    PHI_ITER         = PHI_ITER(savedIters);
    DELTA_PHI_ITER   = DELTA_PHI_ITER(savedIters);
    TAU_SHIFT_ITER   = TAU_SHIFT_ITER(savedIters);
    R_RING_ITER      = R_RING_ITER(savedIters);
    GAMMA_TRUST_ITER = GAMMA_TRUST_ITER(savedIters);
    ALPHA_GATE_ITER  = ALPHA_GATE_ITER(savedIters);
    LAMBDA_POLAR_EFF_ITER = LAMBDA_POLAR_EFF_ITER(savedIters);
    LAMBDA_EDGE_EFF_ITER  = LAMBDA_EDGE_EFF_ITER(savedIters);
    EDGE_DAMP_EFF_ITER    = EDGE_DAMP_EFF_ITER(savedIters);
    EDGE_CONTAM_ITER      = EDGE_CONTAM_ITER(savedIters);
    EDGE_BLEND_EFF_ITER   = EDGE_BLEND_EFF_ITER(savedIters);
end

% Save coarse result
xi = xi_original; yi = yi_original;
VEL_ESTIM_FINAL_NOPP   = VEL_ESTIM;
ATTEN_ESTIM_FINAL_NOPP = ATTEN_ESTIM;

save(filename_results, '-v7.3', ...
    'xi','yi','xi_original','yi_original', ...
    'VEL_ESTIM_FINAL_NOPP','ATTEN_ESTIM_FINAL_NOPP', ...
    'fDATA','niterAttenPerFreq','niterSoSPerFreq', ...
    'VEL_ESTIM_ITER','ATTEN_ESTIM_ITER','GRAD_IMG_ITER','SEARCH_DIR_ITER', ...
    'PHI_ITER','DELTA_PHI_ITER','TAU_SHIFT_ITER','R_RING_ITER', ...
    'GAMMA_TRUST_ITER','ALPHA_GATE_ITER','LAMBDA_POLAR_EFF_ITER', ...
    'LAMBDA_EDGE_EFF_ITER','EDGE_DAMP_EFF_ITER', ...
    'EDGE_CONTAM_ITER','EDGE_BLEND_EFF_ITER', ...
    'savedIters', 'enableVirtualElements', 'useVirtualTx', 'orig_numElements', ...
    'virtualElementsMode', 'exclMode','frac_excl', ...
    'VE_freqMode','VE_freqThresh','VE_manualFlags','VE_flags', ...
    'frac_shift_lo','frac_shift_hi','frac_shift_f_split', ...
    'beta_ve', 'enableIllumComp', 'eps_illum', ...
    'misfitType', 'alpha_hv', 'alpha_hv_atten', ...
    'enableLFPrior', 'lambda_stage', 'f_stage1_cutoff', 'f_stage2_cutoff', ...
    'stage1_ref_saved', 'stage2_ref_saved', ...
    'enableTV', 'lambda_tv', 'eps_tv', ...
    'enablePolarReg', 'lambda_polar', 'eps_polar_reg', ...
    'enableTrustGate','phi_small_drop_ratio','tau_rise_ratio','ring_rise_ratio', ...
    'alpha_gate_badfit','alpha_gate_ringgrowth', ...
    'lambda_polar_gate_boost','lambda_polar_gate_max', ...
    'enableEdgeGuard','f_edge_guard_start', ...
    'edge_r_inner_ratio','edge_r_outer_ratio', ...
    'lambda_edge_anchor','edge_step_damp_max', ...
    'edge_blend_base','edge_blend_max','edge_contam_rise_ratio', ...
    'ring_bg_inner_ratio','ring_bg_outer_ratio', ...
    'mainLoopElapsedSec', 'scriptElapsedSec');

% [修改] 三阶段参考模型追加保存（变量存在时才追加，避免未进入对应阶段时报错）
if exist('VEL_stage1_ref', 'var')
    save(filename_results, '-append', 'VEL_stage1_ref');
end
if exist('VEL_stage2_ref', 'var')
    save(filename_results, '-append', 'VEL_stage2_ref');
end

disp(['Save NOPP completed: ', filename_results]);


%% ============================ Local Functions ============================
function s = ternary_str(cond, s_true, s_false)
    % 三元字符串选择（用于 figure Name 拼接）
    if cond, s = s_true; else, s = s_false; end
end

function w = keysCubicWeights(alpha, a)
% alpha in [0,1], Keys cubic convolution kernel parameter a (often -0.5)
% For position n+alpha, weights correspond to samples [n-1, n, n+1, n+2]
    p = [-1, 0, 1, 2];
    d = abs(alpha - p);
    w = zeros(1,4);
    for i = 1:4
        x = d(i);
        if x < 1
            w(i) = (a+2)*x^3 - (a+3)*x^2 + 1;
        elseif x < 2
            w(i) = a*x^3 - 5*a*x^2 + 8*a*x - 4*a;
        else
            w(i) = 0;
        end
    end
end
