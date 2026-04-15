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
%   12. [修改] 保留三阶段LF先验，其余附加正则/门控模块已移除。
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
verTag                = 'V3_5_noReg_20260414';        % 保存版本号
% ----------------------------------------------------------------------

% ------------------ fws: 邻频高斯平滑衔接（f_k -> f_{k+1}） ------------------
% 目的：仅在三阶段边界（阶段1→2、阶段2→3）做一次高斯低通，
%       作为下一阶段起始频率的初值，避免每个频率都平滑导致累积过平滑。
% 注意：仅在频率切换点生效；不在单频迭代内平滑梯度。
enableInterFreqGaussian = true;         % true=在相邻频率之间做高斯平滑衔接
gaussSigmaScale        = 0.15;          % sigma倍率（阶段边界平滑建议0.10~0.20）
saveInterFreqInitMat   = true;          % true=导出每次邻频衔接模型到mat文件
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
%   PolarPhase: 极坐标相位加权失配（Polar-decomposed Phase-weighted misfit）
%   将复数 phasor 分解为幅度残差 ΔA（→衰减）和相位残差 Δφ（→声速），
%   伴随源 = (1-α)·ΔA·φ_dir + α·Δφ·i·φ_dir
%   参考思路来自 HV geometry（Neumann & Yang 2025），但实现为轻量极坐标近似。
misfitType = 'PolarPhase';   % 'L2' | 'PolarPhase'

% ------------------ fws: 优化器设置（新增） ------------------
% 用户反馈“更新太小、速度图几乎不变”时，可用CG而非纯梯度下降。
% 这里统一启用非线性共轭梯度（PR+ / FR 截断），并保留步长安全夹紧：
%   1) alpha_ls  : 线性化线搜索步长（与原L2流程一致）
%   2) alpha_cap : 按期望 max|Δv| 估计的上限，避免一步过大导致发散
% 最终 alpha = min(alpha_ls, alpha_cap)，兼顾“能动起来”与稳定性。
optimizerType      = 'CG_PR_FR';  % 目前支持: 'CG_PR_FR'
target_dv_per_iter = [];          % []或<=0: 自动估计；>0: 手动指定每步max|Δv| [m/s]
auto_dv_ratio      = 3e-3;        % 自动模式: 目标Δv = auto_dv_ratio*mean(VEL_ESTIM)（建议2e-3~8e-3）
auto_dv_min        = 1.0;         % 自动模式下限 [m/s]
auto_dv_max        = 50.0;        % 自动模式上限 [m/s]
alpha_floor        = 1e-8;        % 步长下限，避免数值退化为零步
% --------------------------------------------------------------

% alpha_hv：SoS 更新阶段的相位权重
%   0.0 = 纯幅度（退化为 L2 幅度版）
%   1.0 = 纯相位（对周期跳跃最鲁棒，但完全忽略幅度匹配）
%   建议起点：0.5（折中相位/幅度约束；消融可测试 [0.3, 0.5, 0.7, 0.9]）
alpha_hv = 0.6;

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
beta_ve = 1.0;
% -------------------------------------------------------------------



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
lambda_stage     = 5e-3;     % 各阶段约束强度（归一化后：0.1%数据梯度；消融可测试 [5e-4, 1e-3, 5e-3]）
f_stage1_cutoff  = 0.475e6;  % 阶段1/2分界 [Hz]（≈Low段末，与 frac_shift_f_split 对齐）
f_stage2_cutoff  = 0.850e6;  % 阶段2/3分界 [Hz]（≈Mid段末，对应 fDATA_SoS 第二组末尾）
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


% 粗网格迭代结果
VEL_ESTIM_ITER   = zeros(Nyi, Nxi, Niter);
ATTEN_ESTIM_ITER = zeros(Nyi, Nxi, Niter);
GRAD_IMG_ITER    = zeros(Nyi, Nxi, Niter);
SEARCH_DIR_ITER  = zeros(Nyi, Nxi, Niter);
PHI_ITER         = nan(1, Niter);

%% ---------- 用户控制：迭代执行选项 ----------
runAllIterations    = true;
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
            grad_stage1_ref  = gradient_img_prev;  % 记录低频阶段末梯度方向参考
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
                case 'polarphase'    
                   amp_sim  = abs(p_sim_sc) + eps;           % |p_sim|，加 eps 防零除
                   phi_dir  = p_sim_sc ./ amp_sim;           % exp(i·φ_sim)，单位相量
                   amp_res  = amp_sim - abs(REC_all);         % ΔA [实数]
                   dphi     = angle(p_sim_sc .* conj(REC_all)); % Δφ ∈ (-π,π] [实数]
                   residual = (1 - alpha_hv_cur) .* amp_res .* phi_dir + ...
                                alpha_hv_cur      .* dphi    .* (1i * phi_dir);
                otherwise
                    error('Unknown misfitType: %s. Use ''L2'' or ''PolarPhase''.', misfitType);
            end

            phi_k_accum = phi_k_accum + sum(abs(residual).^2);

            ADJ_SRC_elmt = zeros(Nyi, Nxi);

            % fws: VE通道置信度加权（方向A2）
            % 真实Rx权重=1，虚拟Rx权重=beta_ve，降低插值偏差对伴随源注入的污染
            % geo_orig套装（isRealRx全1）下此操作等价于乘以全1，无额外开销
            w_channel = ones(sum(rx_mask_all), 1);
            w_channel(~isRealRx_cur(rx_mask_all)) = beta_ve;
            ADJ_SRC_elmt(idx_all) = w_channel .* residual;

            ADJ_SRC(:,:,elmt_idx) = ADJ_SRC_elmt;
        end

        % Backproject error
        ADJ_WVFIELD = HS.solve(ADJ_SRC, true);
        SCALING     = repmat(reshape(scaling,[1,1,numel(scaling)]), [Nyi,Nxi,1]);
        BACKPROJ    = -real(conj(SCALING.*VIRT_SRC).*ADJ_WVFIELD);

        for elmt_idx = 1:numel(tx_include_cur)
            gradient_img = gradient_img + BACKPROJ(:,:,elmt_idx);
        end

        % Remove Ringing from Data Gradient (先滤数据梯度，再叠加正则，避免正则被低通抹掉)
        gradient_img = ringingRemovalFilt(xi_original, yi_original, ...
            gradient_img, c0, fDATA(f_idx), cutoff, ord);
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

        % -----------------------------------------------------------------------

        % -----------------------------------------------------------------------

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
        % 使用 PR+ / FR 混合截断的非线性共轭梯度（NLCG）：
        %   betaPR = g_k^T(g_k-g_{k-1}) / ||g_{k-1}||^2
        %   betaFR = ||g_k||^2 / ||g_{k-1}||^2
        %   beta   = min(max(betaPR,0), betaFR)
        % 与纯最速下降(beta=0)相比，NLCG可在谷底方向累积“惯性”，提升有效更新幅度。
        if do_reset_CG || ~strcmpi(optimizerType, 'CG_PR_FR')
            beta = 0;
        else
            den_cg = real(gradient_img_prev(:)'*gradient_img_prev(:));
            if den_cg <= eps
                beta = 0;
            else
                betaPR = real(gradient_img(:)'*(gradient_img(:)-gradient_img_prev(:))) / den_cg;
                betaFR = real(gradient_img(:)'*gradient_img(:)) / den_cg;
                beta   = min(max(betaPR,0), betaFR);
            end
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

        % 统一线搜索：先做线性化步长，再按目标 max|Δv| 做上限夹紧。
        den_ls   = real(dREC_SIM(:)'*dREC_SIM(:));
        num_ls   = -real(gradient_img(:)'*search_dir(:));
        alpha_ls = num_ls / (den_ls + eps);

        % PolarPhase 线搜索分母修正：
        % 分子来自 PolarPhase 梯度，若直接复用 L2 的 dREC_SIM 二范数作为分母，
        % 会导致相位/幅度混合量纲不匹配。这里在不改变现有流程的前提下，
        % 通过 alpha_hv_cur 的加权平方和对分母做一次近似缩放修正。
        if strcmpi(misfitType, 'PolarPhase')
            dPhi_sq_accum = 0;
            for elmt_idx = 1:numel(tx_include_cur)
                dREC_elmt = dREC_SIM(elmt_idx, :).';
                dPhi_sq_accum = dPhi_sq_accum + sum(abs(dREC_elmt).^2);
            end
            scale_pp  = (1 - alpha_hv_cur)^2 + alpha_hv_cur^2;
            den_ls_pp = dPhi_sq_accum * scale_pp;
            alpha_ls  = num_ls / (den_ls_pp + eps);
        end

        % 步长上限: 将 max|Δv| 控制在 target_dv_eff 内（手动或自动）
        vmean = mean(VEL_ESTIM(:));
        if isempty(target_dv_per_iter) || (target_dv_per_iter <= 0)
            target_dv_eff = min(max(auto_dv_ratio * vmean, auto_dv_min), auto_dv_max);
        else
            target_dv_eff = target_dv_per_iter;
        end
        sd_max      = max(abs(search_dir(:))) + eps;
        alpha_cap   = target_dv_eff / ((vmean^2) * sd_max + eps);
        alpha       = min(alpha_ls, alpha_cap);
        alpha       = max(alpha, alpha_floor);

        % fws: 调试 —— 打印本轮最大速度更新量，确认步长合理
        % Δs = α·search_dir → Δv ≈ v² · Δs（一阶近似）
        delta_s_max  = perc_step_size * alpha * max(abs(search_dir(:)));
        delta_v_max  = mean(VEL_ESTIM(:))^2 * delta_s_max;
        fprintf('  Φ=%.3e | alpha=%.2e | max|Δv|≈%.3f m/s\n', phi_k, alpha, delta_v_max);
        
        % Update slowness
        if updateAttenuation
            SI = sign(sign_conv) * imag(SLOW_ESTIM) + perc_step_size * alpha * search_dir;
            SLOW_ESTIM = real(SLOW_ESTIM) + 1i * sign(sign_conv) * SI;
        else
            SLOW_ESTIM = SLOW_ESTIM + perc_step_size * alpha * search_dir;
        end
        VEL_ESTIM   = 1./real(SLOW_ESTIM);
        ATTEN_ESTIM = 2*pi*imag(SLOW_ESTIM)*sign(sign_conv);
        SLOW_ESTIM  = 1./VEL_ESTIM + 1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi);

        % Save coarse iteration
        VEL_ESTIM_ITER(:,:,iter)   = VEL_ESTIM;
        ATTEN_ESTIM_ITER(:,:,iter) = ATTEN_ESTIM;
        GRAD_IMG_ITER(:,:,iter)    = gradient_img;
        SEARCH_DIR_ITER(:,:,iter)  = search_dir;
        PHI_ITER(iter)             = phi_k;
        savedIters(end+1)          = iter; %#ok<SAGROW>


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

    % 邻频高斯平滑衔接：仅在三阶段边界触发（阶段1->2、阶段2->3）
    is_stage_boundary = enableInterFreqGaussian && (f_idx < numel(fDATA)) && ...
        ((fDATA(f_idx) <= f_stage1_cutoff && fDATA(f_idx + 1) > f_stage1_cutoff) || ...
         (fDATA(f_idx) <= f_stage2_cutoff && fDATA(f_idx + 1) > f_stage2_cutoff));

    if is_stage_boundary
        f_cur = fDATA(f_idx);
        f_next = fDATA(f_idx + 1);
        lambda_cur = c_init / f_cur;
        sigma_gauss = gaussSigmaScale * lambda_cur / (2 * dxi_original);
        r_g = max(1, ceil(3 * sigma_gauss));
        [xx_g, yy_g] = meshgrid(-r_g:r_g, -r_g:r_g);
        ker_g = exp(-(xx_g.^2 + yy_g.^2) / (2 * sigma_gauss^2));
        ker_g = ker_g / (sum(ker_g(:)) + eps);

        % 仅在阵元环内做平滑，避免PML边界伪影向成像域扩散
        R_trans = max(sqrt(transducerPositionsX.^2 + transducerPositionsY.^2));
        smooth_mask = (R_grid <= R_trans * 0.98);  % 留少量边缘余量

        VEL_smooth = conv2(VEL_ESTIM, ker_g, 'same');
        VEL_ESTIM = VEL_ESTIM .* (~smooth_mask) + VEL_smooth .* smooth_mask;
        SLOW_ESTIM = 1./VEL_ESTIM + 1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi);

        if saveInterFreqInitMat
            inter_freq_init_path = sprintf('%s%s_%s_VEL_init_f%.3f_to_f%.3fMHz.mat', ...
                result_dir, filename, verTag, f_cur/1e6, f_next/1e6);
            save(inter_freq_init_path, 'VEL_ESTIM', 'sigma_gauss', 'lambda_cur', 'f_cur', 'f_next', 'dxi_original');
        end

        fprintf(['[InterFreq] 已生成邻频高斯初值: f=%.3f -> %.3f MHz | ' ...
                 'lambda=%.3f mm, sigma=%.3f grid | mask内平滑\n'], ...
                f_cur/1e6, f_next/1e6, lambda_cur*1e3, sigma_gauss);
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
plot(xi_original, vel_prof_est,  'r-',  'LineWidth', 1.8, 'DisplayName', 'Estimated');

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

%% ===================== Quantitative Metrics (PSNR / RMSE / SSIM / RD) =====================
% 说明：
% 1) 为保证与重建网格一致，先将真值 C 插值到 (xi_original, yi_original)
% 2) 目标/背景区域由真值相对背景声速偏差自动生成（分位数阈值）
% 3) RD 采用图中定义的 ROI 向量范数形式：
%    RD = ||c_recon - c_true|| / ||c_bkgnd - c_true|| * 100%

% --- 可调参数（用于目标/背景掩膜自动生成） ---
targetDevPct = 85;   % 偏差绝对值 >= P85 视为目标区域
bgDevPct     = 40;   % 偏差绝对值 <= P40 视为背景区域
bg_ref_tol   = 5.0;  % [m/s] 背景参考容差（用于兜底）

% --- 将真值插值到重建网格 ---
[X_true, Y_true]     = meshgrid(xi_orig, yi_orig);
[X_recon, Y_recon]   = meshgrid(xi_original, yi_original);
C_true_on_reconGrid  = interp2(X_true, Y_true, C, X_recon, Y_recon, 'linear', nan);

% 仅在有效区域计算指标
valid_mask = isfinite(C_true_on_reconGrid) & isfinite(VEL_ESTIM);
v_true  = C_true_on_reconGrid(valid_mask);
v_recon = VEL_ESTIM(valid_mask);

% --- PSNR ---
mse_vel = mean((v_recon - v_true).^2);
peakVel = max(v_true);
if mse_vel > 0
    PSNR_dB = 10 * log10((peakVel^2) / mse_vel);
else
    PSNR_dB = inf;
end

% --- SSIM（全局统计形式，避免依赖工具箱） ---
L_dyn = max(v_true) - min(v_true);
if L_dyn <= 0
    SSIM_val = 1.0;
else
    C1_ssim = (0.01 * L_dyn)^2;
    C2_ssim = (0.03 * L_dyn)^2;
    mu_t = mean(v_true);
    mu_r = mean(v_recon);
    var_t = var(v_true, 1);  % 使用总体方差
    var_r = var(v_recon, 1);
    cov_tr = mean((v_true - mu_t) .* (v_recon - mu_r));
    SSIM_val = ((2*mu_t*mu_r + C1_ssim) * (2*cov_tr + C2_ssim)) / ...
               ((mu_t^2 + mu_r^2 + C1_ssim) * (var_t + var_r + C2_ssim));
end

% --- RD 需要目标与背景区域（CNR已替换为RMSE） ---
bg_ref = median(C_true_on_reconGrid(valid_mask), 'omitnan');
dev_true = abs(C_true_on_reconGrid - bg_ref);
dev_valid = dev_true(valid_mask);
dev_target_th = prctile(dev_valid, targetDevPct);
dev_bg_th     = prctile(dev_valid, bgDevPct);

target_mask = valid_mask & (dev_true >= dev_target_th);
bg_mask     = valid_mask & (dev_true <= dev_bg_th);

% 兜底：若分位数导致空集，退化为固定阈值法
if nnz(target_mask) < 10
    target_mask = valid_mask & (abs(C_true_on_reconGrid - bg_ref) > bg_ref_tol);
end
if nnz(bg_mask) < 10
    bg_mask = valid_mask & (abs(C_true_on_reconGrid - bg_ref) <= bg_ref_tol);
end

RMSE_val = sqrt(mse_vel);

if nnz(target_mask) > 1 && nnz(bg_mask) > 1
    target_recon = VEL_ESTIM(target_mask);
    target_true  = C_true_on_reconGrid(target_mask);
    c_bkgnd_scalar = mean(VEL_ESTIM(bg_mask));

    numerator_RD   = norm(target_recon - target_true);
    denominator_RD = norm(c_bkgnd_scalar - target_true);
    RD_percent = numerator_RD / max(denominator_RD, eps) * 100;
else
    RD_percent = nan;
end

fprintf('\n================ 定量指标 (SoS) ================\n');
fprintf('PSNR : %.4f dB\n', PSNR_dB);
fprintf('RMSE : %.6f m/s\n', RMSE_val);
fprintf('SSIM : %.6f\n', SSIM_val);
fprintf('RD   : %.4f %%\n', RD_percent);
fprintf('目标像素数: %d, 背景像素数: %d\n', nnz(target_mask), nnz(bg_mask));
fprintf('===============================================\n\n');

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
    'PHI_ITER', ...
    'savedIters', 'enableVirtualElements', 'useVirtualTx', 'orig_numElements', ...
    'virtualElementsMode', 'exclMode','frac_excl', ...
    'VE_freqMode','VE_freqThresh','VE_manualFlags','VE_flags', ...
    'frac_shift_lo','frac_shift_hi','frac_shift_f_split', ...
    'misfitType', 'alpha_hv', 'alpha_hv_atten', ...
    'enableLFPrior', 'lambda_stage', 'f_stage1_cutoff', 'f_stage2_cutoff', ...
    'stage1_ref_saved', 'stage2_ref_saved', ...
    'mainLoopElapsedSec', 'scriptElapsedSec', ...
    'PSNR_dB','RMSE_val','SSIM_val','RD_percent', ...
    'targetDevPct','bgDevPct','bg_ref_tol');

% [修改] 三阶段参考模型追加保存（变量存在时才追加，避免未进入对应阶段时报错）
if exist('VEL_stage1_ref', 'var')
    save(filename_results, '-append', 'VEL_stage1_ref');
end
if exist('grad_stage1_ref', 'var')
    save(filename_results, '-append', 'grad_stage1_ref');
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
