% ===================== Waveform Inversion + Virtual Elements =========GIF=================
%  版本: FreqVE
%  新增特性: 按频率动态切换虚拟阵元（VE_freqMode: always / threshold / manual）
%  修改说明:
%    1. 顶部新增 VE_freqMode / VE_freqThresh / VE_manualFlags 控制参数
%    2. VE代码块改为预计算 geo_orig / geo_ve 两套几何，不再覆盖全局变量
%    3. 主循环内按 VE_flags(f_idx) 动态选择当前套装，跨套装切换时强制重置CG
%    4. GIF 预创建窗口，复用图像句柄，大幅减少每帧开销
%    5. [修改] 已删除 PP 模块及其相关代码
%    6. [新增] 方向A2：VE通道置信度加权（beta_ve_sos_schedule / beta_ve_atten_schedule），降低虚拟Rx插值偏差对梯度的污染
%    7. [新增] 失配函数：misfitStage_* 三变量按大阶段调度（始终启用），
%             Stage1 SoS-only / Stage2 SoS更新 / Stage2 Atten更新分别独立指定；
%             全程统一 misfit 时三变量设为相同值即可，无需额外开关。
%             PolarPhase 将复数 phasor 分解为幅度残差（→衰减）和相位残差（→声速），
%             分别用 alpha_hv_src / alpha_hv_ls 控制，缓解 L2 的周期跳跃敏感性。
%             参考：Neumann & Yang, arXiv:2505.01817 (2025)，频域FWI + USCT验证。
%    8. [修改] 方向E ：低频先验正则化升级为三阶段版本（enableLFPrior），
%             对应师兄时域FWI"均匀分三段低频复用，效果较好"的频域等价实现：
%               阶段1（低频，< f_stage1_cutoff）：无先验约束，自由重建背景声速
%               阶段2（中频，f_stage1~f_stage2）：约束拉向阶段1末尾参考模型
%               阶段3（高频，> f_stage2_cutoff）：约束拉向阶段2末尾参考模型
%             各阶段内 lambda_k = lambda_stage·(f_stageN_cutoff/f)²，频率越高约束越弱。
%             相对归一化消除量纲不匹配（见主循环内注释）。
%    9. [修改] 已删除方向D/F/G及正则化诊断打印相关代码
%   10. 其余算法逻辑（Helmholtz、CG）与原版一致
%   11. [修改] kWave 步长改成 KCI 风格：PolarPhase 使用 ΔA/Δφ 自洽线搜索
%              （保留文本1的平方权重形式，对应“残差层面加权后取平方范数”）
%   12. [修改] Near-end exclusion 改为 KCI manual 版本，支持 manual / angle / frac 三种模式
%   13. [新增] 定量指标同时计算 SoS 与 Attenuation，并分开显示/保存
%   14. [新增] GIF 改为 SoS + Attenuation 同图保存，便于同步观察声速/衰减迭代演化
%   15. [修改] 新增按 misfit 类型选择完整迭代链：PolarPhase走当前链，L2回到VE+L2原版链
%   16. [修改] 将 alpha_hv 拆分为两个独立参数，语义更清晰：
%              alpha_hv_src ：SoS阶段伴随源的相位权重（控制梯度方向的相位/幅度比例）
%              alpha_hv_ls  ：PolarPhase 自洽线搜索分母权重（控制步长大小，与梯度方向解耦）
%   17. [修改] VE通道降权与misfit链解耦：beta_ve_cur 现在对 L2 链和 PolarPhase 链均生效
%   18. [冗余修复] 删除顶层 misfitType 变量和 enableStageMisfitSchedule 开关：
%              两者在原版中同时存在时 misfitType 被完全覆盖，属于逻辑冗余。
%              现在始终使用三个 misfitStage_* 变量按大阶段调度；
%              若想全程统一 misfit，将三变量设为相同值即可。
%   19. [冗余修复] 删除从未调用的 ternary_str 本地函数
%   20. [冗余修复] 删除单独固定 beta_ve_sos / beta_ve_atten 变量（SoS/Atten统一使用schedule）
%   21. [冗余修复] perc_step_size 移至主循环外（循环内从不改变）
%   22. [新增] SoS/Atten阶段 beta_ve 均按频率分段调度：固定beta时只需把三段填成同一值；
%              Atten默认低频保留VE稳定化，中高频逐渐减弱虚拟Rx幅度误差固化。
%   23. [删除] 原 frac_shift_atten_* / enableFracShiftSchedule_Att / geo_ve_atten_lo/hi
%              以及 frac_excl 的 SoS/Atten 分离调度逻辑——
%              物理原因：SoS 与 Atten 共享同一目标函数 J(c,α)，
%              若两者从不同阵元集合 A_sos ≠ A_atten 分别求偏导，
%              实际是在优化两个不自洽的子目标，把各自更新叠加到同一介质上。
%              SoS路径一变，Atten幅度拟合的参考态也随之改变，
%              导致 Δd_amp 不再反映真实衰减误差，而是混入 c 更新不一致带来的幅度偏差。
%              保留统一几何/数据集；差异化需求改由 Atten梯度角度软加权（方向A3）实现。
%   24. [新增] 方向A3：Atten 梯度角度软加权（enableAttenAngleWeight）
%              物理依据：SoS 与 Atten 共用同一套数据 A 和同一套 frac_excl hard aperture，
%              但 Atten 主要由幅度残差驱动，对近端直达波/界面散射/虚拟通道幅度误差更敏感。
%              为避免再引入 A_sos ≠ A_atten 的硬数据集分离，本文改为在 Atten 伴随源注入层面
%              对每个 Tx-Rx 对施加基于角向间隔 dtheta 的软权重：
%                dtheta = |wrap(theta_rx - theta_tx)|,  dtheta ∈ [0, π]
%                w_atten_sr = w_min + (1-w_min)·sin(dtheta/2)^attenAngleWeightExp
%              近端接收器 dtheta 小 → 权重接近 w_min；远端透射路径 dtheta≈π → 权重接近 1。
%              这样既保持 SoS/Atten 使用同一套阵元剔除方案，又能软性降低 Atten 对近端幅度污染的敏感性。
%              attenAngleWeightExp=0 → 等权；=2 → 平滑二次递增（推荐起点）；=4 → 更强调远端透射路径。
%   25. [新增] 分阶段任意失配函数频率调度（enableMisfitFreqSchedule_*）：
%              默认 Stage1 SoS-only / Stage2 SoS / Stage2 Atten 仍由 misfitStage_* 决定，
%              同时允许各阶段按"频率段"和"特定频率点"覆盖为指定方案。
%              每条规则均可指定 'L2' 或 'PolarPhase'，因此不再限定只能局部切到 L2。
%              支持 PolarPhase→L2→PolarPhase、L2→PolarPhase→L2 等任意夹心式消融调度。
%   26. [新增] Near-end exclusion 按频率三阶段动态调度：
%              frac_excl_schedule = [0.50, 0.25, 0.12]，对应低频更保守、中频逐步开放、
%              高频保留少量近端排除，避免真实数据近端反射/散射直接污染梯度。
%              SoS 与 Atten 共用同一套动态 frac_excl（保持目标函数一致性）。
%   27. [修改] Atten阶段PolarPhase改为衰减特异的 log-amplitude 残差：
%              SoS阶段仍使用线性幅度残差 ΔA = |p_sim|-|p_obs|；
%              Atten阶段使用 ΔlogA = log|p_sim|-log|p_obs|，更贴近衰减的指数幅度损失。
%              同步修改PolarPhase自洽线搜索中的幅度项一阶导数 a_amp = a/|p_sim|，
%              避免"残差是log幅度、线搜索仍按线性幅度"的不自洽。
%   28. [冗余修复] 删除旧变量名兼容字段：
%              不再保留 alpha_hv_atten / L2FreqWindow / LowFreqL2Warmup / 旧指标名 / 旧frac_shift字段；
%              仅保存当前主逻辑真正使用的变量。
%   29. [新增] 定量指标加入阵元内部 ROI 掩膜：
%              主指标只在阵元环内、向内收缩后的有效成像区域计算，避免阵元外/PML/边界水域
%              稀释RMSE或污染SSIM/RD；ROI参数和掩膜同步保存，便于后处理复现。
%   30. [新增] X/Y轮廓线示意图叠加主定量ROI虚线圆：
%              profile取样线与PSNR/RMSE/SSIM/RD实际统计区域在同一图中显示，避免指标区域歧义。
%   31. [修改] profile绘图拆分：
%              figure(5) 作为轮廓示意图，仅显示True与最终Estimated两条曲线；
%              figure(6) 单独显示两阶段参考profile，包含True / Estimated / Stage1 ref / Stage2 ref四条曲线。
%   32. [修改] L2链线搜索改为加权自洽分母：
%              由于 ADJ_SRC 中已使用 w_channel .* residual，L2分母同步使用
%              REC_WEIGHT_CUR .* |Jd|²，使梯度方向与线搜索曲率对应同一个加权L2目标函数。
%   33. [修改] PolarPhase线搜索通道权重统一接入 REC_WEIGHT_CUR：
%              SoS阶段体现VE置信度权重，Atten阶段体现VE置信度 + A3角度软权重。
%              L2/PolarPhase 两条链的伴随源与线搜索均使用同一通道权重。
%   34. [修改] 按文本2逻辑删除 StepCap / alpha_floor：
%              不再用 max|Δv| 人工夹紧步长，max|Δv| 仅保留为日志诊断量；
%              若 alpha_ls 非法或非正，则本轮 alpha=0，直接跳过更新。
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
verTag                = 'V_ours_011_noCap_wLS';  % 保存版本号（wLS=加权自洽线搜索）
% ----------------------------------------------------------------------

% ------------------ fws: 邻频高斯平滑衔接（f_k -> f_{k+1}） ------------------
% 目的：仅在三阶段边界（阶段1→2、阶段2→3）做一次高斯低通，
%       作为下一阶段起始频率的初值，避免每个频率都平滑导致累积过平滑。
% 注意：仅在频率切换点生效；不在单频迭代内平滑梯度。
enableInterFreqGaussian = true;         % true=在相邻频率之间做高斯平滑衔接
gaussSigmaScale        = 0.15;          % sigma倍率（阶段边界平滑建议0.10~0.20）
saveInterFreqInitMat   = false;          % true=导出每次邻频衔接模型到mat文件
% ----------------------------------------------------------------------

% ------------------ Near-end exclusion controls ------------------
% 'manual' : exclude a fixed number of receivers on each side of the transmitter
% 'angle'  : use Ali-equivalent acceptance angle（保留原方案，当前注释不用）
% 'frac'   : exclude a fixed fraction of receivers nearest the transmitter（保留原方案，当前注释不用）
exclMode = 'frac';   % 'manual' | 'angle' | 'frac'
numElemLeftRightExcl = 34; % Manual exclusion count on each side
frac_excl = 1/4;       % only used when exclMode='frac' and enableFracExclStageSchedule=false

% ------------------ fws: 逐阶开放接收阵元（动态 frac_excl 调度） ------------------
% 目的：低频阶段先保守使用远端透射/散射信息，避免近端强直达/界面反射污染背景更新；
%       随频率升高逐步开放更多接收阵元，提高角向约束和局部细节恢复能力。
%
% 逻辑：仅在 exclMode='frac' 时生效；manual / angle 模式保持原有定义。
%   Stage 1: f <= stage_freq_boundary(1)  → frac_excl_schedule(1)
%   Stage 2: f <= stage_freq_boundary(2)  → frac_excl_schedule(2)
%   Stage 3: f >  stage_freq_boundary(2)  → frac_excl_schedule(3)
%
% [修改] SoS 与 Atten 共用同一套动态 frac_excl：
%        保持目标函数 J(c,α) 对同一套数据 A 求偏导的一致性。
%        外围阵元对 Atten 幅度梯度的差异化处理改由 enableAttenAngleWeight 软加权实现。
enableFracExclStageSchedule = true;
frac_excl_schedule          = [0.30, 0.25, 0.18];       % Stage1 / Stage2 / Stage3
stage_freq_boundary         = [0.650e6, 1.000e6];      % [Stage1上界, Stage2上界]，单位Hz

if numel(frac_excl_schedule) ~= 3 || numel(stage_freq_boundary) ~= 2
    error('frac_excl_schedule 必须为3个值，stage_freq_boundary 必须为2个频率边界。');
end
if any(frac_excl_schedule < 0) || any(frac_excl_schedule >= 1)
    error('frac_excl_schedule 取值应满足 0 <= frac_excl < 1。');
end
if stage_freq_boundary(1) > stage_freq_boundary(2)
    error('stage_freq_boundary 必须按从低到高排列。');
end
% -------------------------------------------------------------

% ------------------ fws: 失配函数阶段调度（始终启用） ------------------
% [冗余修复] 原版中 misfitType + enableStageMisfitSchedule=true 的组合
% 导致 misfitType 被完全覆盖，属于逻辑冗余。
% 现在直接用三个 misfitStage_* 变量按反演大阶段调度：
%   misfitStage_SoSOnly      ：Stage1 SoS-only 阶段（fDATA_SoS 段）的 SoS 更新
%   misfitStage_SoSAtten_SoS ：Stage2 SoS+Atten 阶段中的 SoS 子迭代
%   misfitStage_SoSAtten_Att ：Stage2 SoS+Atten 阶段中的 Atten 子迭代
%
% 若想全程统一 misfit，将三变量设为相同值即可，无需额外开关。
% 例：全程 L2    → 三变量均设 'L2'
%     全程 PolarPhase → 三变量均设 'PolarPhase'
misfitStage_SoSOnly      = 'PolarPhase';  % Stage1：SoS-only，PolarPhase改善声速
misfitStage_SoSAtten_SoS = 'L2';          % Stage2：SoS+Atten里的SoS更新默认L2
misfitStage_SoSAtten_Att = 'L2';          % Stage2：Atten更新默认L2
% -----------------------------------------------------------------------

% ------------------ fws: 分阶段任意失配函数频率调度（频点 / 频段） ------------------
% 目的：把"局部频率使用 L2"的硬编码，升级为"局部频率使用指定 misfit"的通用调度表。
%       默认方案仍由上方 misfitStage_* 三变量决定；只有命中下方规则时，才临时覆盖。
%
% 支持两类覆盖：
%   1) 频率段覆盖：misfitFreqRanges_* 的每一行写成 {[f_min, f_max], '方案名'}
%   2) 频率点覆盖：misfitFreqPoints_* 的每一行写成 {[f1, f2, ...], '方案名'}
%
% 支持的方案名：'L2' | 'PolarPhase'
%
% 优先级：misfitStage_* 默认 < 频率段表（后写覆盖前写）< 频率点表（最高优先级）
misfitFreqTol_Hz = 1.0;  % 频率匹配容差 [Hz]，避免浮点边界误差

% Stage1 SoS-only
enableMisfitFreqSchedule_SoSOnly = false;
misfitFreqRanges_SoSOnly = {[0.300e6, 0.400e6], 'L2';};
misfitFreqPoints_SoSOnly = {};

% Stage2 SoS子迭代 
enableMisfitFreqSchedule_SoSAtten_SoS = false;
misfitFreqRanges_SoSAtten_SoS = {[0.325e6, 0.625e6], 'PolarPhase';};
misfitFreqPoints_SoSAtten_SoS = {};

% Stage2 Atten子迭代
enableMisfitFreqSchedule_SoSAtten_Att = false;
misfitFreqRanges_SoSAtten_Att = {[0.775e6, 0.950e6], 'PolarPhase';};
misfitFreqPoints_SoSAtten_Att = {};

validateMisfitScheduleTable(misfitFreqRanges_SoSOnly,      'misfitFreqRanges_SoSOnly',      true);
validateMisfitScheduleTable(misfitFreqPoints_SoSOnly,      'misfitFreqPoints_SoSOnly',      false);
validateMisfitScheduleTable(misfitFreqRanges_SoSAtten_SoS, 'misfitFreqRanges_SoSAtten_SoS', true);
validateMisfitScheduleTable(misfitFreqPoints_SoSAtten_SoS, 'misfitFreqPoints_SoSAtten_SoS', false);
validateMisfitScheduleTable(misfitFreqRanges_SoSAtten_Att, 'misfitFreqRanges_SoSAtten_Att', true);
validateMisfitScheduleTable(misfitFreqPoints_SoSAtten_Att, 'misfitFreqPoints_SoSAtten_Att', false);
% -----------------------------------------------------------------------

% ------------------ fws: 按 misfit 类型选择完整迭代链 ------------------
% PolarPhase：当前流程（平方权重ΔA/Δφ自洽线搜索 + TrustGate/CrossFreqTrust/EdgeGuard等门控）
% L2        ：VE+加权L2流程（原CG + 加权自洽线搜索 + 原慢度更新）
% VE通道降权（beta_ve_cur）两条链均生效，并同步进入对应线搜索分母
useOriginalVEChainForL2 = true;

% ------------------ fws: PolarPhase 伴随源与线搜索权重设置 ------------------
% alpha_hv_src ：SoS阶段伴随源的相位权重（控制梯度方向的相位/幅度比例）
%   = 1.0 → 纯相位伴随源；∈ (0,1) → 混合伴随源
%   建议消融测试 [0.6, 0.8, 1.0]；当前默认 1.0（纯相位）
%
% alpha_hv_ls  ：PolarPhase 自洽线搜索分母权重（控制步长大小，与梯度方向解耦）
alpha_hv_src = 0.8;   % SoS阶段：伴随源相位权重（1.0=纯相位；<1.0=混合）
alpha_hv_ls  = 0.6;   % PolarPhase 自洽线搜索目标权重（控制步长，与梯度方向解耦）

% alpha_hv_atten_src：衰减更新阶段的伴随源相位权重（控制Atten梯度方向）
%   衰减主要由幅度约束，应保持较低相位权重；
%   = 0.00 → 纯幅度方向；0.05~0.15 → 幅度主导 + 少量相位定位
alpha_hv_atten_src = 0.2;

% alpha_hv_atten_ls：衰减更新阶段的PolarPhase线搜索权重（控制Atten步长曲率估计）
alpha_hv_atten_ls  = 0.4;

% polarLogAmpFloorRel：Atten阶段 log-amplitude 残差的幅度下限比例
%   避免REC_DATA中被outlier处理置零的通道导致 log(0) 产生极大残差
polarLogAmpFloorRel = 1e-6;
% -------------------------------------------------------

% ------------------ fws: 方向A3：Atten梯度角度软加权 ------------------
% 物理依据：SoS 与 Atten 应从同一套数据 A 求偏导（保持目标函数一致性），
%           但 Atten 主要受幅度残差控制，对近端直达波、界面散射和虚拟通道幅度插值误差更敏感。
%           因此不再为 Atten 单独更换 hard receiver mask，而是在同一 elemInclude_cur 内
%           对 Atten 的每个 Tx-Rx 通道施加基于角向间隔 dtheta 的软权重。
%
% 角度权重定义：
%   dtheta = |wrap(theta_rx - theta_tx)|,  dtheta ∈ [0, π]
%   w_raw  = sin(dtheta/2)^attenAngleWeightExp
%   w_att  = w_min + (1 - w_min) * w_raw
%
% 物理含义：
%   dtheta 小       ：接收器靠近发射器，近端直达/界面散射/幅度建模误差更容易污染 Atten，权重接近 w_min
%   dtheta 接近 π   ：远端透射路径，幅度信息更接近有效衰减约束，权重接近 1
%
% 参数说明：
%   enableAttenAngleWeight      = true  → 开启软加权；false → 等权（原版消融对比）
%   attenAngleWeightMode        = 'dtheta' → 使用 Tx-Rx 角向间隔权重（当前推荐）
%   attenAngleWeightExp         = 2     → 平滑二次递增；=4 更强调远端透射路径；=0 等权
%   attenAngleWeightMinSchedule = [低频, 中频, 高频] 最小保留权重，避免 Atten 梯度被压得过淡
%
% 注意：该权重在 ADJ_SRC 构建时逐 Rx 通道施加，只作用于 Atten 梯度；
%       SoS 梯度不施加该角度加权，仍保留宽孔径相位约束。
%       SoS 和 Atten 仍共用同一套前向场、同一套 elemInclude_cur 和同一套 VE 几何。
enableAttenAngleWeight      = true;           % true=开启；false=等权（原版消融对比）
attenAngleWeightMode        = 'dtheta';       % 'dtheta'=按 Tx-Rx 角向间隔加权（当前推荐）
attenAngleWeightExp         = 2;              % sin(dtheta/2)^n，建议从 2 开始消融
attenAngleWeightMinSchedule = [0.70, 0.80, 0.90];  % Atten角度软权重下限：[低频, 中频, 高频]

if numel(attenAngleWeightMinSchedule) ~= 3
    error('attenAngleWeightMinSchedule 必须为 [低频, 中频, 高频] 三个值。');
end
if any(attenAngleWeightMinSchedule < 0) || any(attenAngleWeightMinSchedule > 1)
    error('attenAngleWeightMinSchedule 取值应满足 0 <= w_min <= 1。');
end
if attenAngleWeightExp < 0
    error('attenAngleWeightExp 必须 >= 0。');
end
if ~strcmpi(attenAngleWeightMode, 'dtheta')
    error('当前版本仅支持 attenAngleWeightMode=''dtheta''。');
end
% -----------------------------------------------------------------------
% ------------------ fws: 优化器设置 ------------------
% 1) L2：使用与 ADJ_SRC 中 w_channel .* residual 一致的加权线性化线搜索，
%    alpha_ls = -<g,d> / sum(w_sr * |Jd|²)，使梯度方向与步长曲率自洽。
% 2) PolarPhase：保留文本1的平方权重 ΔA / Δφ 自洽线搜索，
%    即残差层面加权后取平方范数，分子/分母均使用 (1-α)^2 与 α^2，
%    且通道权重 w_sr 与 ADJ_SRC 中 w_channel 保持一致。
% 3) [修改] 删除 max|Δv| 人工夹紧：
%    不再使用 enableStepCap / target_dv_per_iter / auto_dv_ratio / alpha_cap / alpha_floor；
%    max|Δv| 仅作为日志诊断量，不参与步长控制。
% 4) 若 alpha_ls 非法或非正，则本轮 alpha=0，跳过模型更新，不再强行给 alpha_floor。
optimizerType      = 'CG_PR_FR';  % 目前支持: 'CG_PR_FR'
% --------------------------------------------------------------

% ------------------ fws: VE通道置信度加权（方向A2）------------------
% 物理依据：虚拟Rx通道由插值生成，不是独立测量，其残差包含插值误差引入的偏差。
%           对真实Rx保持单位权重，对虚拟Rx施加 beta_ve < 1 的降权，
%           使伴随源注入中虚假相关成分减少，梯度更准确地反映真实物理约束。
% [修改] beta_ve_cur 现对 SoS / Atten 两阶段均真正生效：
%        两条链（L2 / PolarPhase）均通过 w_channel 施加 beta_ve_cur 加权。
%        beta_ve_cur=1.0 时等价于原版等权注入。
% [冗余修复] 删除单独固定 beta_ve_sos / beta_ve_atten 变量：
%        现在 SoS 与 Atten 都直接使用三段 schedule；若想固定某阶段 beta，
%        将对应 schedule 三个值填成相同值即可，例如 [0.50,0.50,0.50]。
% 取值范围：1.0 = 等权（与原版等价）；0.0 = 仅真实Rx参与梯度（最保守）
% 建议：SoS 可先固定 [1.00,1.00,1.00]；Atten 可测试 [0.50,0.25,0.10]
%
% 频率段边界：
%   低频：f <= beta_ve_edges(1)
%   中频：beta_ve_edges(1) < f <= beta_ve_edges(2)
%   高频：f > beta_ve_edges(2)
beta_ve_edges          = [0.450e6, 0.850e6];     % <=0.45MHz低频，0.45~0.85MHz中频，>0.85MHz高频
beta_ve_sos_schedule   = [1.00, 1.00, 1.00];    % SoS阶段：[低频, 中频, 高频]；三段相同=固定beta
beta_ve_atten_schedule = [0.90, 0.50, 0.30];    % Atten阶段：[低频, 中频, 高频]；三段相同=固定beta

if numel(beta_ve_edges) ~= 2 || numel(beta_ve_sos_schedule) ~= 3 || numel(beta_ve_atten_schedule) ~= 3
    error('beta_ve_edges 必须为2个边界，beta_ve_*_schedule 必须为 [低频, 中频, 高频] 三个值。');
end
% -------------------------------------------------------------------

% [修改] ------------------ fws: 三阶段低频先验正则化（方向E）------------------
% 物理依据：对应师兄时域FWI"均匀分三个阶段进行低频复用，效果较好"的频域等价。
%           阶段1（低频）：无先验约束，自由重建背景声速
%           阶段2（中频）：约束拉向 VEL_stage1_ref（阶段1末尾）
%           阶段3（高频）：约束拉向 VEL_stage2_ref（阶段2末尾）
%           约束强度：lambda_k = lambda_stage · (f_stageN_cutoff / fDATA(f_idx))²
%           相对归一化：lambda_stage 为占数据梯度的比例（无量纲）。
enableLFPrior    = true;     % true=开启三阶段LF先验；false=关闭
lambda_stage     = 5e-3;     % 各阶段约束强度（消融可测试 [5e-4, 1e-3, 5e-3]）
f_stage1_cutoff  = 0.475e6;  % 阶段1/2分界 [Hz]
f_stage2_cutoff  = 0.850e6;  % 阶段2/3分界 [Hz]
% -----------------------------------------------------------------------

% ------------------ fws: 高频更新可信度监控与门控 ------------------
enableTrustGate        = true;
phi_small_drop_ratio   = 1e-3;  % "失配下降很小"的相对阈值
tau_rise_ratio         = 1.01;  % τ_k > 1.01*τ_{k-1} 视为上升
ring_rise_ratio        = 1.02;  % R_k > 1.02*R_{k-1} 视为环伪影增强
alpha_gate_badfit      = 0.50;  % 拟合退化时步长缩放
alpha_gate_ringgrowth  = 0.70;  % 环伪影增长时额外缩放

% ------------------ fws: 跨频一致性门控 ------------------
enableCrossFreqTrust   = true;
clf_warn               = 0.10;  % C_lf 警戒阈值
clf_rise_ratio         = 1.05;  % C_lf 较上一步上升>5% 视为恶化
sigma_bg_warn          = 5.0;   % 纯水外环标准差警戒 [m/s]
cos_grad_warn          = 0.30;  % 梯度相似度警戒阈值
alpha_gate_xf_min      = 0.35;  % 跨频门控最小保留步长比例

water_bg_inner_ratio   = 0.93;  % 纯水外环统计起点
water_bg_outer_ratio   = 0.97;  % 纯水外环统计终点
% -----------------------------------------------------------------------

% ------------------ fws: 高频边界保护（抑制"中心污染外圈"） ------------------
enableEdgeGuard    = false;
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
saveGIF        = true;   % true=保存迭代 GIF；false=跳过
gif_delay      = 0.20;    % 每帧停留时间 [s]
gif_save_every = 5;       % 每隔多少次迭代保存一帧
% -----------------------------------------------

% Red-Blue colormap (blue=low, red=high)
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

% fws-Modified: Multi-stage step size with same overall frequency range
fDATA_SoS = [ (0.30:0.025:0.45), (0.50:0.05:0.70), (0.775:0.075:1.225), 1.250 ] * 1e6;
fDATA_SoSAtten = fDATA_SoS + 25e3;  % same +25 kHz shift as original
fDATA = [fDATA_SoS, fDATA_SoSAtten]; % all frequencies [Hz]

% [调试] 频率表校验
fprintf('[CHECK] numel(fDATA_SoS)=%d | SoS总迭代=%d | Niter预估=%d\n', ...
    numel(fDATA_SoS), numel(fDATA_SoS)*3, numel(fDATA_SoS)*3 + numel(fDATA_SoSAtten)*6);

% ------------------ fws: VE几何偏移（统一单套，SoS/Atten共用） ------------------
% [修改] 删除 frac_shift_atten_* 及 enableFracShiftSchedule_Att：
%        SoS 与 Atten 共用同一套 frac_shift 几何，保持目标函数一致性。
%        Atten 对外围路径的差异化需求由 enableAttenAngleWeight 软加权实现。
%
% SoS阶段 frac_shift（SoS/Atten共用）：低频/高频两档
frac_shift_sos_lo      = 0.10;     % 低频段偏移（f < frac_shift_sos_f_split）
frac_shift_sos_hi      = 0.20;     % 高频段偏移（f >= frac_shift_sos_f_split）
frac_shift_sos_f_split = 0.450e6;  % 低/高频分界 [Hz]

if ~isscalar(frac_shift_sos_lo) || ~isscalar(frac_shift_sos_hi) || ~isscalar(frac_shift_sos_f_split)
    error('frac_shift_sos_* 必须为标量。');
end

% ------------------ VE 频率开关 ------------------
VE_freqThresh  = NaN;    % 仅 threshold 模式使用
VE_freqMode    = 'always';
VE_manualFlags = mod(0:numel(fDATA)-1, 2) == 1;
% -------------------------------------------------

% Iteration schedule
niterSoSPerFreq   = [3*ones(size(fDATA_SoS)),   3*ones(size(fDATA_SoSAtten))];
niterAttenPerFreq = [0*ones(size(fDATA_SoS)),   3*ones(size(fDATA_SoSAtten))];

% DTFT (not FFT)
DTFT = exp(-1i*2*pi*fDATA'*time)*mean(diff(time));

% Geometric TOF Based on Sound Speed
c_geom = 1540; % Sound Speed [m/s]
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
dxi  = 0.3e-3; xmax = 120e-3;
xi   = -xmax:dxi:xmax; yi = xi;
Nxi  = numel(xi); Nyi = numel(yi);
[Xi, Yi] = meshgrid(xi, yi);
R_grid = sqrt(Xi.^2 + Yi.^2);
% Discretize transducer ring onto grid
x_circ = transducerPositionsX;
y_circ = transducerPositionsY;
x_idx  = dsearchn(xi(:), x_circ(:));
y_idx  = dsearchn(yi(:), y_circ(:));
ind    = sub2ind([Nyi, Nxi], y_idx, x_idx);
msk    = zeros(Nyi, Nxi); msk(ind) = 1;

% 备份原始粗网格
xi_original = xi; yi_original = yi; dxi_original = dxi;


%% ===================== Parameters =====================
h = dxi; % [m]
g = 1;
alphaDB = 0.0;
alphaNp = (log(10)/20)*alphaDB*((1e3)/(1e6));
ATTEN = alphaNp*ones(Nyi,Nxi);
sign_conv = -1;
a0 = 10;
L_PML = 9.0e-3;

Np2dB = 20/log(10);
slow2atten = (1e6)/(1e2);

crange = [1350, 1600];
attenrange = 10*[-1,1];

%% ===================== Phase Screen =====================
x_circ_disc = xi(x_idx); y_circ_disc = yi(y_idx);
[x_circ_disc_tx, x_circ_disc_rx] = meshgrid(x_circ_disc, x_circ_disc);
[y_circ_disc_tx, y_circ_disc_rx] = meshgrid(y_circ_disc, y_circ_disc);
geomTOFs_disc  = sqrt((x_circ_disc_tx-x_circ_disc_rx).^2 + ...
                      (y_circ_disc_tx-y_circ_disc_rx).^2)/c_geom;
geomTOFs_error = geomTOFs_disc - geomTOFs;

PS = zeros(numElements, numElements, numel(fDATA));
for f_idx = 1:numel(fDATA)
    PS(:,:,f_idx) = exp(1i*sign_conv*2*pi*fDATA(f_idx)*geomTOFs_error);
end
REC_DATA = REC_DATA .* PS;
clearvars PS geomTOFs_error geomTOFs_disc geomTOFs;

REC_DATA(isnan(REC_DATA)) = 0;

%% ===================== dwnsmp =====================
dwnsmp = 1;

%% ===================== 预计算几何套装 =====================
%  geo_orig : 原始128阵元套装
%  geo_ve_lo / geo_ve_hi : VE扩展，按 frac_shift_sos_lo / frac_shift_sos_hi 两档
%  [修改] 删除 geo_ve_atten_lo/hi（Atten与SoS共用几何套装）
%  SoS/Atten 共用同一套 frac_shift 几何，差异化需求由 Atten梯度角度软加权实现

numElem_Ali     = 512;
numExclSide_Ali = 63;
theta_excl      = 1.2 * numExclSide_Ali * (2*pi/numElem_Ali);

% 预清理使用最大开放度对应的 frac_excl
if enableFracExclStageSchedule && strcmpi(exclMode, 'frac')
    frac_excl_precompute = min(frac_excl_schedule);
else
    frac_excl_precompute = frac_excl;
end

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
theta_all_orig   = atan2(transducerPositionsY(:), transducerPositionsX(:));
[elemInclude_orig, numElemLeftRightExcl_orig, keepRx_orig] = buildNearEndExclusionMask( ...
    orig_numElements, exclMode, numElemLeftRightExcl, frac_excl_precompute, theta_excl, theta_all_orig);

switch lower(exclMode)
    case 'manual'
        fprintf('[EXCL预清理] mode=manual | exclude per side=%d | keep=%d/%d receivers\n', ...
            numElemLeftRightExcl_orig, keepRx_orig, orig_numElements);
    case 'angle'
        fprintf('[EXCL预清理] mode=angle | theta_excl=%.2f deg | keep min=%d/%d receivers\n', ...
            theta_excl*180/pi, keepRx_orig, orig_numElements);
    case 'frac'
        fprintf('[EXCL预清理] mode=frac | frac_excl_precompute=%.3f | exclude per side=%d | keep=%d/%d receivers\n', ...
            frac_excl_precompute, numElemLeftRightExcl_orig, keepRx_orig, orig_numElements);
    otherwise
        error('Unknown exclMode: %s.', exclMode);
end
geo_orig.elemInclude = elemInclude_orig;
geo_orig.tx_include  = 1:dwnsmp:orig_numElements;

% Remove Outliers
perc_outliers = 0.99;
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
%  套装B：Virtual Elements Scheme（统一两档 frac_shift，SoS/Atten共用）
%  [修改] 只预计算 SoS frac_shift 两档；Atten 直接复用相同套装。
% ============================================================
if enableVirtualElements

    virtualInterpKernel = 'zoh';   % 'zoh' | 'linear' | 'cubic'
    keys_a = -0.5;

    fs_all_unique = unique([frac_shift_sos_lo, frac_shift_sos_hi]);
    fprintf('[VE预计算] 需要构建 %d 套不同frac_shift几何（lo=%.2f，hi=%.2f）\n', ...
        numel(fs_all_unique), frac_shift_sos_lo, frac_shift_sos_hi);

    geo_ve_cache = cell(1, numel(fs_all_unique));

    for fsi = 1:numel(fs_all_unique)
        fs_cur = fs_all_unique(fsi);

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
                upsample_factor = 2;
                sgn = ones(size(theta_base)); sgn(2:2:end) = -1;
                dtheta_use = dtheta_next; dtheta_use(sgn<0) = dtheta_prev(sgn<0);
                theta_virtual = theta_base + sgn.*fs_cur.*dtheta_use;
                x_virtual = r_base.*cos(theta_virtual);
                y_virtual = r_base.*sin(theta_virtual);
                tx_ve = reshape([x_base; x_virtual],1,[]);
                ty_ve = reshape([y_base; y_virtual],1,[]);

            case 'bilateral'
                upsample_factor = 3;
                theta_left  = theta_base - fs_cur.*dtheta_prev;
                theta_right = theta_base + fs_cur.*dtheta_next;
                x_left  = r_base.*cos(theta_left);  y_left  = r_base.*sin(theta_left);
                x_right = r_base.*cos(theta_right); y_right = r_base.*sin(theta_right);
                tx_ve = reshape([x_left; x_base; x_right],1,[]);
                ty_ve = reshape([y_left; y_base; y_right],1,[]);

            otherwise
                error('Unknown virtualElementsMode. Use ''unilateral'' or ''bilateral''.');
        end

        L_ve = upsample_factor;
        N_ve = N0;
        numElements_ve = L_ve * N_ve;
        base_map = repelem(1:N_ve, L_ve);

        x_idx_ve = dsearchn(xi(:), tx_ve(:));
        y_idx_ve = dsearchn(yi(:), ty_ve(:));
        ind_ve   = sub2ind([Nyi, Nxi], y_idx_ve, x_idx_ve);

        % 构造 H_theta
        H_theta = zeros(L_ve*N_ve, N_ve, 'single');
        for n = 1:N_ve
            switch lower(virtualElementsMode)
                case 'unilateral'
                    H_theta(2*n-1, n) = 1;
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
                    H_theta(3*n-1, n) = 1;
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

        % 角向卷积重采样
        REC_DATA_ve_full = zeros(L_ve*N_ve, L_ve*N_ve, numel(fDATA), 'like', REC_DATA);
        for f_idx = 1:numel(fDATA)
            D0 = REC_DATA(:,:,f_idx);
            REC_DATA_ve_full(:,:,f_idx) = H_theta * D0 * H_theta.';
        end

        % 离散几何TOF差相位搬运
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

        % 真实Rx掩码
        isRealRx_ve = false(1, numElements_ve);
        switch lower(virtualElementsMode)
            case 'unilateral'; isRealRx_ve(1:L_ve:end) = true;
            case 'bilateral';  isRealRx_ve(2:L_ve:end) = true;
        end

        % tx_include
        if useVirtualTx
            tx_include_ve = 1:dwnsmp:numElements_ve;                   %#ok<UNRCH>
        else
            switch lower(virtualElementsMode)
                case 'unilateral'; tx_include_ve = 1:L_ve:numElements_ve;
                case 'bilateral';  tx_include_ve = 2:L_ve:numElements_ve;
            end
            tx_include_ve = tx_include_ve(1:dwnsmp:end);
        end

        % Near-end exclusion（VE版本）
        theta_all_ve   = atan2(ty_ve(:), tx_ve(:));
        [elemInclude_ve, numElemLeftRightExcl_ve, keepRx_ve] = buildNearEndExclusionMask( ...
            numElements_ve, exclMode, numElemLeftRightExcl, frac_excl_precompute, theta_excl, theta_all_ve);

        switch lower(exclMode)
            case 'manual'
                fprintf('  [EXCL预清理][VE %.2f] mode=manual | exclude per side=%d | keep=%d/%d receivers\n', ...
                    fs_cur, numElemLeftRightExcl_ve, keepRx_ve, numElements_ve);
            case 'angle'
                fprintf('  [EXCL预清理][VE %.2f] mode=angle | theta_excl=%.2f deg | keep min=%d/%d receivers\n', ...
                    fs_cur, theta_excl*180/pi, keepRx_ve, numElements_ve);
            case 'frac'
                fprintf('  [EXCL预清理][VE %.2f] mode=frac | frac_excl_precompute=%.3f | exclude per side=%d | keep=%d/%d receivers\n', ...
                    fs_cur, frac_excl_precompute, numElemLeftRightExcl_ve, keepRx_ve, numElements_ve);
            otherwise
                error('Unknown exclMode: %s.', exclMode);
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
        g_tmp.frac_shift           = fs_cur;

        geo_ve_cache{fsi} = g_tmp;
        fprintf('  [VE预计算] frac_shift=%.2f 完成（套装 %d/%d）\n', fs_cur, fsi, numel(fs_all_unique));
    end

    % 从缓存取两档套装
    [~, idx_sos_lo] = min(abs(fs_all_unique - frac_shift_sos_lo));
    [~, idx_sos_hi] = min(abs(fs_all_unique - frac_shift_sos_hi));
    geo_ve_lo = geo_ve_cache{idx_sos_lo};   % frac_shift = frac_shift_sos_lo
    geo_ve_hi = geo_ve_cache{idx_sos_hi};   % frac_shift = frac_shift_sos_hi

    % [修改] Atten 与 SoS 共用同一套几何，无需额外 geo_ve_atten_*
    fprintf('[VE分配] SoS/Atten统一共用: lo=%.2f, hi=%.2f（分界=%.3f MHz）\n', ...
        frac_shift_sos_lo, frac_shift_sos_hi, frac_shift_sos_f_split/1e6);
    fprintf('[VE分配] Atten差异化需求由 enableAttenAngleWeight=%d（mode=%s, exp=%.1f, w_min=[%.2f %.2f %.2f]）实现\n', ...
        enableAttenAngleWeight, attenAngleWeightMode, attenAngleWeightExp, ...
        attenAngleWeightMinSchedule(1), attenAngleWeightMinSchedule(2), attenAngleWeightMinSchedule(3));
end

%% ===================== 构造 VE_flags 向量 =====================
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

% [新增] 打印逐阶开放接收阵元预览
fprintf('\n=== Near-end exclusion 动态调度预览（SoS/Atten统一共用）===\n');
fprintf('mode=%s | enableFracExclStageSchedule=%d | boundary=[%.3f, %.3f] MHz | schedule=[%.2f, %.2f, %.2f]\n', ...
    exclMode, enableFracExclStageSchedule, stage_freq_boundary(1)/1e6, stage_freq_boundary(2)/1e6, ...
    frac_excl_schedule(1), frac_excl_schedule(2), frac_excl_schedule(3));
fprintf('%-10s  %-10s  %-10s  %-12s  %-12s  %-12s\n', ...
    'freq(MHz)', 'stage', 'frac', 'excl/side', 'keepRx', 'numRx');
for fi = 1:numel(fDATA)
    if VE_flags(fi) && enableVirtualElements
        numElements_preview = geo_ve_lo.numElements;
    else
        numElements_preview = geo_orig.numElements;
    end
    [frac_excl_preview, frac_excl_stage_idx_preview, frac_excl_stage_preview] = getFracExclByFreq( ...
        fDATA(fi), enableFracExclStageSchedule, frac_excl, frac_excl_schedule, stage_freq_boundary);
    if strcmpi(exclMode, 'frac')
        numElemLeftRightExcl_preview = round(numElements_preview * frac_excl_preview / 2);
        keepRx_preview = numElements_preview - (2*numElemLeftRightExcl_preview + 1);
        frac_str_preview = sprintf('%.2f', frac_excl_preview);
        excl_side_str_preview = sprintf('%d', numElemLeftRightExcl_preview);
        keep_str_preview = sprintf('%d', keepRx_preview);
    else
        frac_str_preview = '-';
        excl_side_str_preview = '-';
        keep_str_preview = '-';
    end
    fprintf('  %-8.3f  %-10s  %-10s  %-12s  %-12s  %-12d\n', ...
        fDATA(fi)/1e6, frac_excl_stage_preview, frac_str_preview, ...
        excl_side_str_preview, keep_str_preview, numElements_preview);
end
fprintf('======================\n\n');

% [修改] 打印VE分配预览（统一套装，SoS/Atten共用）
fprintf('\n=== VE 频率分配预览（SoS/Atten统一共用同一套几何）===\n');
fprintf('模式: %s | frac_shift: %.2f(<%.2fMHz) / %.2f(>=%.2fMHz)\n', ...
    VE_freqMode, frac_shift_sos_lo, frac_shift_sos_f_split/1e6, ...
    frac_shift_sos_hi, frac_shift_sos_f_split/1e6);
fprintf('Atten差异化: enableAttenAngleWeight=%d | mode=%s | exp=%.1f | w_min=[%.2f, %.2f, %.2f]\n', ...
    enableAttenAngleWeight, attenAngleWeightMode, attenAngleWeightExp, ...
    attenAngleWeightMinSchedule(1), attenAngleWeightMinSchedule(2), attenAngleWeightMinSchedule(3));
fprintf('beta SoS schedule: [%.2f, %.2f, %.2f] | beta Atten schedule: [%.2f, %.2f, %.2f]\n', ...
    beta_ve_sos_schedule(1), beta_ve_sos_schedule(2), beta_ve_sos_schedule(3), ...
    beta_ve_atten_schedule(1), beta_ve_atten_schedule(2), beta_ve_atten_schedule(3));
fprintf('%-10s  %-10s  %-6s  %-12s  %-12s  %-12s\n', ...
    'freq(MHz)','type','VE','fs','betaSoS','betaAtt');

for fi = 1:numel(fDATA)
    if niterAttenPerFreq(fi) == 0, it = 'SoS'; else, it = 'SoS+Att'; end

    if VE_flags(fi) && enableVirtualElements
        if fDATA(fi) < frac_shift_sos_f_split
            fs_disp = frac_shift_sos_lo;
        else
            fs_disp = frac_shift_sos_hi;
        end
        if fDATA(fi) <= beta_ve_edges(1)
            beta_stage_idx_disp = 1;
        elseif fDATA(fi) <= beta_ve_edges(2)
            beta_stage_idx_disp = 2;
        else
            beta_stage_idx_disp = 3;
        end
        beta_sos_str = sprintf('%.2f', beta_ve_sos_schedule(beta_stage_idx_disp));
        if niterAttenPerFreq(fi) > 0
            beta_att_str = sprintf('%.2f', beta_ve_atten_schedule(beta_stage_idx_disp));
        else
            beta_att_str = '-';
        end
    else
        fs_disp      = 0;
        beta_sos_str = '-';
        beta_att_str = '-';
    end

    fprintf('  %-8.3f  %-10s  %d     %-12s  %-12s  %-12s\n', ...
        fDATA(fi)/1e6, it, VE_flags(fi), ...
        sprintf('%.2f', fs_disp), beta_sos_str, beta_att_str);
end
fprintf('======================\n\n');

% 打印失配函数调度预览
fprintf('\n=== Misfit 调度预览 ===\n');
fprintf('Stage默认方案: SoSOnly=%s | SoSAtten-SoS=%s | SoSAtten-Att=%s\n', ...
    misfitStage_SoSOnly, misfitStage_SoSAtten_SoS, misfitStage_SoSAtten_Att);
fprintf('频率覆盖调度: Stage1 enable=%d | Stage2-SoS enable=%d | Stage2-Att enable=%d | tol=%.1f Hz\n', ...
    enableMisfitFreqSchedule_SoSOnly, enableMisfitFreqSchedule_SoSAtten_SoS, ...
    enableMisfitFreqSchedule_SoSAtten_Att, misfitFreqTol_Hz);
fprintf('%-10s  %-14s  %-42s\n', 'freq(MHz)', 'stage', 'misfit');
for fi = 1:numel(fDATA)
    if niterAttenPerFreq(fi) == 0
        stage_preview = 'Stage1-SoS';
        [misfit_preview_cur, override_preview, override_mode_preview, override_reason_preview] = ...
            applyMisfitFreqSchedule(fDATA(fi), misfitStage_SoSOnly, ...
                enableMisfitFreqSchedule_SoSOnly, misfitFreqRanges_SoSOnly, ...
                misfitFreqPoints_SoSOnly, misfitFreqTol_Hz);
        misfit_preview = formatMisfitPreview(misfit_preview_cur, override_preview, ...
            override_mode_preview, override_reason_preview);
    else
        stage_preview = 'Stage2-SoS/Att';
        [misfit_preview_sos_cur, override_sos_preview, override_sos_mode_preview, override_sos_reason_preview] = ...
            applyMisfitFreqSchedule(fDATA(fi), misfitStage_SoSAtten_SoS, ...
                enableMisfitFreqSchedule_SoSAtten_SoS, misfitFreqRanges_SoSAtten_SoS, ...
                misfitFreqPoints_SoSAtten_SoS, misfitFreqTol_Hz);
        misfit_preview_sos = formatMisfitPreview(misfit_preview_sos_cur, override_sos_preview, ...
            override_sos_mode_preview, override_sos_reason_preview);
        [misfit_preview_att_cur, override_att_preview, override_att_mode_preview, override_att_reason_preview] = ...
            applyMisfitFreqSchedule(fDATA(fi), misfitStage_SoSAtten_Att, ...
                enableMisfitFreqSchedule_SoSAtten_Att, misfitFreqRanges_SoSAtten_Att, ...
                misfitFreqPoints_SoSAtten_Att, misfitFreqTol_Hz);
        misfit_preview_att = formatMisfitPreview(misfit_preview_att_cur, override_att_preview, ...
            override_att_mode_preview, override_att_reason_preview);
        misfit_preview = sprintf('SoS=%s, Att=%s', misfit_preview_sos, misfit_preview_att);
    end
    fprintf('  %-8.3f  %-14s  %-42s\n', fDATA(fi)/1e6, stage_preview, misfit_preview);
end
fprintf('======================\n\n');

%% ===================== Initialize Model =====================
c_init = 1480;
VEL_INIT   = c_init*ones(Nyi,Nxi);
ATTEN_INIT = 0*alphaNp*ones(Nyi,Nxi);

search_dir = zeros(Nyi,Nxi);
gradient_img_prev = zeros(Nyi,Nxi);

VEL_ESTIM   = VEL_INIT;
ATTEN_ESTIM = ATTEN_INIT;
SLOW_ESTIM  = 1./VEL_ESTIM + 1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi);

c0 = mean(VEL_ESTIM(:)); cutoff = 0.75; ord = Inf;
Niter = sum(niterSoSPerFreq) + sum(niterAttenPerFreq);

% 背景环带（环伪影比 R_k 统计）
Rmax = max(R_grid(:));
ring_bg_inner_ratio = 0.65;
ring_bg_outer_ratio = 0.92;
ring_bg_mask = (R_grid >= ring_bg_inner_ratio*Rmax) & (R_grid <= ring_bg_outer_ratio*Rmax);
% 外环保护软掩膜
edge_guard_mask = (R_grid - edge_r_inner_ratio*Rmax) ./ ...
                  max((edge_r_outer_ratio - edge_r_inner_ratio)*Rmax, eps);
edge_guard_mask = min(max(edge_guard_mask, 0), 1);
edge_outer_mask = edge_guard_mask > 0.95;
pure_water_mask = (R_grid >= water_bg_inner_ratio*Rmax) & ...
                  (R_grid <= water_bg_outer_ratio*Rmax);

% 迭代历史
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
LAMBDA_EDGE_EFF_ITER  = nan(1, Niter);
EDGE_DAMP_EFF_ITER    = nan(1, Niter);
EDGE_CONTAM_ITER      = nan(1, Niter);
EDGE_BLEND_EFF_ITER   = nan(1, Niter);
C_LF_ITER             = nan(1, Niter);
MU_BG_ITER            = nan(1, Niter);
SIGMA_BG_ITER         = nan(1, Niter);
COS_GRAD_ITER         = nan(1, Niter);
ALPHA_GATE_XF_ITER    = ones(1, Niter);
BETA_VE_ITER                 = nan(1, Niter);
BETA_VE_STAGE_IDX_ITER       = nan(1, Niter);
FRAC_SHIFT_ITER              = nan(1, Niter);
FRAC_SHIFT_STAGE_IDX_ITER    = nan(1, Niter);
FRAC_EXCL_ITER               = nan(1, Niter);
FRAC_EXCL_STAGE_IDX_ITER     = nan(1, Niter);
NUM_RX_EXCL_SIDE_ITER        = nan(1, Niter);
NUM_RX_KEEP_ITER             = nan(1, Niter);
ATTEN_ANGLE_WEIGHT_ITER      = nan(1, Niter);   % [新增] 本轮 Atten 梯度角度权重的均值（监控用）
ATTEN_ANGLE_WEIGHT_MIN_ITER  = nan(1, Niter);   % [新增] 本轮 Atten 角度软权重下限（监控用）
MISFIT_OVERRIDE_ITER         = false(1, Niter);
MISFIT_OVERRIDE_STAGE_ITER   = cell(1, Niter);
MISFIT_OVERRIDE_MODE_ITER    = cell(1, Niter);
MISFIT_OVERRIDE_REASON_ITER  = cell(1, Niter);
MISFIT_TYPE_ITER             = cell(1, Niter);

%% ---------- 用户控制：迭代执行选项 ----------
runAllIterations    = true;
requestedIterations = 60;
if requestedIterations >= Niter, runAllIterations = true; end
savedIters   = [];
stopNow      = false;
prev_VE_flag = false;
stage1_ref_saved = false;
stage2_ref_saved = false;
prev_phi_k       = nan;
prev_tau_k       = nan;
prev_ring_k      = nan;
prev_clf_k       = nan;
prev_pred_drop   = nan;
prev_edge_contam = nan;
% -----------------------------------------------

% perc_step_size 移至主循环外
perc_step_size = 1;

%% ---------- 路径 & GIF 初始化 ----------
result_dir = 'D:\Document_ING_fws\WaveformInversionUST\Results\start20260513\';
if ~exist(result_dir,'dir'), mkdir(result_dir); end

gif_filepath = '';
if saveGIF
    gif_filepath = [result_dir, filename, '_', verTag, '.gif'];
    fprintf('[GIF] 将保存至: %s\n', gif_filepath);
    gif_initialized = false;
    gif_iter_count  = 0;

    hgif = figure(99);
    set(hgif, 'Visible','on', 'Position',[80 80 1180 520]);
    clf(hgif);

    tgif = tiledlayout(hgif, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    axgif_vel = nexttile(tgif, 1);
    hgif_vel = imagesc(axgif_vel, xi_original, yi_original, VEL_ESTIM, crange);
    axis(axgif_vel, 'image'); colorbar(axgif_vel); colormap(axgif_vel, cmap_rb);
    xlabel(axgif_vel, 'Lateral [m]'); ylabel(axgif_vel, 'Axial [m]');
    title(axgif_vel, 'SoS [m/s]', 'Interpreter', 'none');

    axgif_att = nexttile(tgif, 2);
    hgif_att = imagesc(axgif_att, xi_original, yi_original, Np2dB*slow2atten*ATTEN_ESTIM, attenrange);
    axis(axgif_att, 'image'); colorbar(axgif_att); colormap(axgif_att, cmap_rb);
    xlabel(axgif_att, 'Lateral [m]'); ylabel(axgif_att, 'Axial [m]');
    title(axgif_att, 'Attenuation [dB/(cm MHz)]', 'Interpreter', 'none');
end

mainLoopTimer = tic;
% -------------------------------------------------------

%% ===================== Main Inversion Loop =====================
for f_idx = 1:numel(fDATA)
    for iter_f_idx = 1:(niterSoSPerFreq(f_idx)+niterAttenPerFreq(f_idx))
        tic;
        iter = iter_f_idx + sum(niterSoSPerFreq(1:f_idx-1)) + sum(niterAttenPerFreq(1:f_idx-1));

        if ~runAllIterations && (numel(savedIters) >= requestedIterations)
            iter    = numel(savedIters);
            stopNow = true;
            disp(['Stopped early after ', num2str(numel(savedIters)), ' iters.']);
            break;
        end

        % 提前计算 updateAttenuation（SoS/Atten阶段判断）
        updateAttenuation = (iter_f_idx > niterSoSPerFreq(f_idx));
        isSoSOnlyStage    = (niterAttenPerFreq(f_idx) == 0);

        % ---- 按频率动态选择几何套装（SoS/Atten统一共用）----
        % [修改] 删除 SoS/Atten 分离逻辑：两者均使用相同的 frac_shift 套装
        if VE_flags(f_idx) && enableVirtualElements
            if fDATA(f_idx) < frac_shift_sos_f_split
                frac_shift_stage_idx_cur = 1;
                frac_shift_stage_cur     = 'lo';
                frac_shift_cur           = frac_shift_sos_lo;
                geo_cur                  = geo_ve_lo;
                cur_fs_tag               = sprintf('VE_lo(%.2f)', frac_shift_cur);
            else
                frac_shift_stage_idx_cur = 2;
                frac_shift_stage_cur     = 'hi';
                frac_shift_cur           = frac_shift_sos_hi;
                geo_cur                  = geo_ve_hi;
                cur_fs_tag               = sprintf('VE_hi(%.2f)', frac_shift_cur);
            end
        else
            geo_cur                  = geo_orig;
            cur_fs_tag               = 'orig';
            frac_shift_cur           = 0;
            frac_shift_stage_idx_cur = 0;
            frac_shift_stage_cur     = 'orig';
        end
        numElements_cur = geo_cur.numElements;
        tx_include_cur  = geo_cur.tx_include;
        ind_cur         = geo_cur.ind;
        isRealRx_cur    = geo_cur.isRealRx;
        REC_DATA_cur    = geo_cur.REC_DATA;

        % [新增] 逐阶开放接收阵元（SoS/Atten统一共用同一套 frac_excl）
        [frac_excl_cur, frac_excl_stage_idx_excl_cur, frac_excl_stage_cur] = getFracExclByFreq( ...
            fDATA(f_idx), enableFracExclStageSchedule, frac_excl, frac_excl_schedule, stage_freq_boundary);
        theta_all_cur = atan2(geo_cur.transducerPositionsY(:), geo_cur.transducerPositionsX(:));
        [elemInclude_cur, numElemLeftRightExcl_cur, numRxKeep_cur] = buildNearEndExclusionMask( ...
            numElements_cur, exclMode, numElemLeftRightExcl, frac_excl_cur, theta_excl, theta_all_cur);

        if strcmpi(exclMode, 'frac')
            cur_fs_tag = sprintf('%s | excl=%s(%.2f)', cur_fs_tag, frac_excl_stage_cur, frac_excl_cur);
        end

        % Reset CG at Each Frequency (SoS and Attenuation)
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
            prev_clf_k        = nan;
            prev_pred_drop    = nan;
            prev_edge_contam  = nan;
        end

        if iter_f_idx == 1
            prev_VE_flag = VE_flags(f_idx);
        end

        % ---- 三阶段LF参考模型保存（方向E）----
        if enableLFPrior && ~stage1_ref_saved && fDATA(f_idx) > f_stage1_cutoff
            VEL_stage1_ref   = VEL_ESTIM;
            grad_stage1_ref  = gradient_img_prev;
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

        gradient_img = zeros(Nyi, Nxi);

        % Sources
        SRC = zeros(Nyi, Nxi, numel(tx_include_cur));
        for elmt_idx = 1:numel(tx_include_cur)
            x_idx_src = geo_cur.x_idx(tx_include_cur(elmt_idx));
            y_idx_src = geo_cur.y_idx(tx_include_cur(elmt_idx));
            SRC(y_idx_src, x_idx_src, elmt_idx) = 1;
        end

        HS = HelmholtzSolver(xi_original, yi_original, VEL_ESTIM, ATTEN_ESTIM, ...
            fDATA(f_idx), sign_conv, a0, L_PML);
        [WVFIELD, VIRT_SRC] = HS.solve(SRC, false);

        if updateAttenuation
            VIRT_SRC = 1i*sign(sign_conv)*VIRT_SRC;
        end

        % 当前 pass 线搜索与伴随源权重
        if updateAttenuation
            alpha_hv_ls_cur  = alpha_hv_atten_ls;
            alpha_hv_src_cur = alpha_hv_atten_src;
        else
            alpha_hv_ls_cur  = alpha_hv_ls;
            alpha_hv_src_cur = alpha_hv_src;
        end

        % VE通道置信度加权（方向A2）
        if fDATA(f_idx) <= beta_ve_edges(1)
            beta_ve_stage_idx_cur = 1;
            beta_ve_stage_cur     = 'low';
        elseif fDATA(f_idx) <= beta_ve_edges(2)
            beta_ve_stage_idx_cur = 2;
            beta_ve_stage_cur     = 'mid';
        else
            beta_ve_stage_idx_cur = 3;
            beta_ve_stage_cur     = 'high';
        end

        if updateAttenuation
            beta_ve_cur = beta_ve_atten_schedule(beta_ve_stage_idx_cur);
        else
            beta_ve_cur = beta_ve_sos_schedule(beta_ve_stage_idx_cur);
        end

        % 失配函数阶段调度
        if isSoSOnlyStage
            misfitType_cur = misfitStage_SoSOnly;
        elseif updateAttenuation
            misfitType_cur = misfitStage_SoSAtten_Att;
        else
            misfitType_cur = misfitStage_SoSAtten_SoS;
        end

        % 分阶段任意 misfit 频率覆盖规则
        isMisfitOverrideCur    = false;
        misfitOverrideStageCur = 'none';
        misfitOverrideModeCur  = 'default';
        misfitOverrideReasonCur = 'default';

        if isSoSOnlyStage && ~updateAttenuation
            [misfitType_cur, isMisfitOverrideCur, misfitOverrideModeCur, misfitOverrideReasonCur] = ...
                applyMisfitFreqSchedule(fDATA(f_idx), misfitType_cur, ...
                    enableMisfitFreqSchedule_SoSOnly, misfitFreqRanges_SoSOnly, ...
                    misfitFreqPoints_SoSOnly, misfitFreqTol_Hz);
            if isMisfitOverrideCur, misfitOverrideStageCur = 'Stage1-SoSOnly'; end

        elseif ~isSoSOnlyStage && ~updateAttenuation
            [misfitType_cur, isMisfitOverrideCur, misfitOverrideModeCur, misfitOverrideReasonCur] = ...
                applyMisfitFreqSchedule(fDATA(f_idx), misfitType_cur, ...
                    enableMisfitFreqSchedule_SoSAtten_SoS, misfitFreqRanges_SoSAtten_SoS, ...
                    misfitFreqPoints_SoSAtten_SoS, misfitFreqTol_Hz);
            if isMisfitOverrideCur, misfitOverrideStageCur = 'Stage2-SoS'; end

        elseif ~isSoSOnlyStage && updateAttenuation
            [misfitType_cur, isMisfitOverrideCur, misfitOverrideModeCur, misfitOverrideReasonCur] = ...
                applyMisfitFreqSchedule(fDATA(f_idx), misfitType_cur, ...
                    enableMisfitFreqSchedule_SoSAtten_Att, misfitFreqRanges_SoSAtten_Att, ...
                    misfitFreqPoints_SoSAtten_Att, misfitFreqTol_Hz);
            if isMisfitOverrideCur, misfitOverrideStageCur = 'Stage2-Atten'; end
        end

        % 完整迭代链选择
        useL2OriginalChain = useOriginalVEChainForL2 && strcmpi(misfitType_cur, 'L2');
        usePolarPhaseChain = strcmpi(misfitType_cur, 'PolarPhase');

        % 当前显示标签
        if usePolarPhaseChain
            if ~updateAttenuation && (alpha_hv_src_cur >= 1.0 - eps)
                misfitLabelCur = sprintf('%s(α_src=1.0/PurePhase,α_ls=%.2f)', ...
                    misfitType_cur, alpha_hv_ls_cur);
            else
                misfitLabelCur = sprintf('%s(α_src=%.2f,α_ls=%.2f)', ...
                    misfitType_cur, alpha_hv_src_cur, alpha_hv_ls_cur);
            end
            if isMisfitOverrideCur
                misfitLabelCur = sprintf('%s-override[%s,%s]', ...
                    misfitLabelCur, misfitOverrideStageCur, misfitOverrideModeCur);
            end
        else
            if isMisfitOverrideCur
                misfitLabelCur = sprintf('%s-override[%s,%s](β=%.2f,%s)', ...
                    misfitType_cur, misfitOverrideStageCur, misfitOverrideModeCur, beta_ve_cur, beta_ve_stage_cur);
            else
                misfitLabelCur = sprintf('%s(β=%.2f,%s)', misfitType_cur, beta_ve_cur, beta_ve_stage_cur);
            end
        end

        if isSoSOnlyStage
            stageLabelCur = 'Stage1 SoS-only';
        elseif updateAttenuation
            stageLabelCur = 'Stage2 Atten update';
        else
            stageLabelCur = 'Stage2 SoS update';
        end
        if isMisfitOverrideCur
            stageLabelCur = [stageLabelCur, ' | misfit override'];
        end

        % ================================================================
        % Build Adjoint Sources
        % ================================================================
        scaling = zeros(numel(tx_include_cur), 1);
        ADJ_SRC = zeros(Nyi, Nxi, numel(tx_include_cur));
        phi_k_accum = 0;
        tau_k_accum = 0;
        tau_k_count = 0;

        % fws: 通道权重缓存（VE置信度 + Atten角度软权重）
        % 只在当前 elemInclude_cur 对应通道填入权重。
        % [修改] REC_WEIGHT_CUR 与 ADJ_SRC 中的 w_channel 保持一致：
        %        - L2链线搜索分母使用同一权重，形成加权L2目标函数；
        %        - PolarPhase链线搜索分子/分母使用同一权重，保持平方权重PP目标函数自洽。
        % 注意：w_channel 作为目标函数通道权重，只使用一次，不在分母中再次平方。
        REC_WEIGHT_CUR = zeros(numel(tx_include_cur), numElements_cur, 'single');

        if strcmpi(misfitType_cur, 'PolarPhase')
            phi_dir_store = zeros(numel(tx_include_cur), numElements_cur, 'like', single(1i));
            amp_sim_store = zeros(numel(tx_include_cur), numElements_cur, 'single');
            amp_res_store = zeros(numel(tx_include_cur), numElements_cur, 'single');
            dphi_store    = zeros(numel(tx_include_cur), numElements_cur, 'single');
        end

        % ================================================================
        % [新增] 方向A3：Atten梯度角度软加权预计算（dtheta版本）
        % 在 Build Adjoint Sources 之前计算每个 Tx-Rx 对的角向间隔权重矩阵。
        %
        % 权重定义：
        %   dtheta = |wrap(theta_rx - theta_tx)|,  dtheta ∈ [0, π]
        %   w_raw  = sin(dtheta/2)^attenAngleWeightExp
        %   w_att  = w_min + (1 - w_min) * w_raw
        %
        % 解释：
        %   dtheta 小     ：Rx靠近Tx，近端直达/界面散射/幅度建模误差更容易污染 Atten，权重接近 w_min
        %   dtheta 接近π ：远端透射路径，幅度信息更接近有效衰减约束，权重接近 1
        %
        % 注意：
        %   该权重不改变 elemInclude_cur，不等价于为 Atten 换一套 hard mask；
        %   它只在已保留通道内部改变 Atten 梯度贡献大小。
        % ================================================================
        atten_angle_weight_mean = 1.0;  % 监控量，仅 Atten 阶段有意义
        atten_angle_weight_min_cur = 1.0;
        if updateAttenuation && enableAttenAngleWeight
            atten_angle_weight_min_cur = attenAngleWeightMinSchedule(beta_ve_stage_idx_cur);

            theta_rx_all = atan2(geo_cur.transducerPositionsY(:).', geo_cur.transducerPositionsX(:).');
            atten_w_matrix = ones(numel(tx_include_cur), numElements_cur, 'single');

            active_w_sum = 0;
            active_w_cnt = 0;

            for elmt_idx = 1:numel(tx_include_cur)
                tx_id = tx_include_cur(elmt_idx);
                theta_tx = theta_rx_all(tx_id);

                dtheta_rx = atan2(sin(theta_rx_all - theta_tx), cos(theta_rx_all - theta_tx));
                dtheta_rx = abs(dtheta_rx);  % [0, pi]

                w_raw = sin(0.5 * dtheta_rx) .^ attenAngleWeightExp;
                w_raw = min(max(w_raw, 0), 1);

                w_att = atten_angle_weight_min_cur + ...
                    (1 - atten_angle_weight_min_cur) .* w_raw;

                atten_w_matrix(elmt_idx, :) = single(w_att);

                rx_mask_preview = elemInclude_cur(tx_id, :);
                active_w_sum = active_w_sum + sum(w_att(rx_mask_preview), 'omitnan');
                active_w_cnt = active_w_cnt + nnz(rx_mask_preview);
            end

            atten_angle_weight_mean = active_w_sum / max(active_w_cnt, 1);

            fprintf('  [AttenAngleWeight] mode=%s | exp=%.1f | w_min=%.2f | mean_w(active)=%.4f\n', ...
                attenAngleWeightMode, attenAngleWeightExp, atten_angle_weight_min_cur, atten_angle_weight_mean);
        end
        % ================================================================

        for elmt_idx = 1:numel(tx_include_cur)
            WVFIELD_elmt = WVFIELD(:,:,elmt_idx);
            tx_id = tx_include_cur(elmt_idx);

            rx_mask_all = elemInclude_cur(tx_id, :);
            idx_all     = ind_cur(rx_mask_all);
            REC_SIM_all = WVFIELD_elmt(idx_all);
            REC_all     = REC_DATA_cur(elmt_idx, rx_mask_all, f_idx);
            REC_all     = REC_all(:);

            % scaling 仅在真实Rx上计算
            rx_mask_scale = rx_mask_all & isRealRx_cur;
            idx_scale     = ind_cur(rx_mask_scale);
            REC_SIM_scale = WVFIELD_elmt(idx_scale);
            REC_scale     = REC_DATA_cur(elmt_idx, rx_mask_scale, f_idx);
            REC_scale     = REC_scale(:);

            scaling(elmt_idx) = (REC_SIM_scale(:)'*REC_scale(:)) / ...
                                (REC_SIM_scale(:)'*REC_SIM_scale(:));

            % 计算失配残差
            p_sim_sc = scaling(elmt_idx) * REC_SIM_all;
            switch lower(misfitType_cur)
                case 'l2'
                    residual = p_sim_sc - REC_all;

                case 'polarphase'
                    amp_sim_raw = abs(p_sim_sc);
                    amp_obs_raw = abs(REC_all);
                    amp_sim     = amp_sim_raw + eps;
                    amp_obs     = amp_obs_raw + eps;
                    phi_dir     = p_sim_sc ./ amp_sim;
                    dphi        = angle(p_sim_sc .* conj(REC_all));

                    if updateAttenuation
                        amp_ref = median(amp_obs_raw(amp_obs_raw > 0), 'omitnan');
                        if ~isfinite(amp_ref) || amp_ref <= 0
                            amp_ref = median(amp_sim_raw(amp_sim_raw > 0), 'omitnan');
                        end
                        if ~isfinite(amp_ref) || amp_ref <= 0, amp_ref = 1; end
                        amp_floor = max(polarLogAmpFloorRel * amp_ref, eps);
                        amp_res   = log(max(amp_sim_raw, amp_floor)) - log(max(amp_obs_raw, amp_floor));
                    else
                        amp_res = amp_sim - amp_obs;
                    end

                    residual = (1 - alpha_hv_src_cur) .* amp_res .* phi_dir + ...
                                alpha_hv_src_cur      .* dphi    .* (1i * phi_dir);

                    phi_dir_store(elmt_idx, rx_mask_all) = phi_dir(:).';
                    amp_sim_store(elmt_idx, rx_mask_all) = amp_sim(:).';
                    amp_res_store(elmt_idx, rx_mask_all) = amp_res(:).';
                    dphi_store(elmt_idx, rx_mask_all)    = dphi(:).';

                otherwise
                    error('Unknown misfitType: %s. Use ''L2'' or ''PolarPhase''.', misfitType_cur);
            end

            phi_k_accum = phi_k_accum + sum(abs(residual).^2);
            dphi_tau    = angle(p_sim_sc .* conj(REC_all));
            tau_k_accum = tau_k_accum + sum(abs(dphi_tau) / (2*pi*fDATA(f_idx)));
            tau_k_count = tau_k_count + numel(dphi_tau);

            ADJ_SRC_elmt = zeros(Nyi, Nxi);

            % ================================================================
            % VE通道置信度加权（方向A2） + Atten角向间隔软加权（方向A3）
            % 两条链（L2 / PolarPhase）均统一通过 w_channel 施加 beta_ve_cur 降权；
            % Atten阶段在同一 elemInclude_cur 内额外乘以 dtheta 软权重，不改变 hard mask。
            % ================================================================
            w_channel = ones(sum(rx_mask_all), 1);
            w_channel(~isRealRx_cur(rx_mask_all)) = beta_ve_cur;

            if updateAttenuation && enableAttenAngleWeight
                w_ang_elmt = atten_w_matrix(elmt_idx, rx_mask_all).';
                w_channel  = w_channel .* w_ang_elmt;
            end

            REC_WEIGHT_CUR(elmt_idx, rx_mask_all) = single(w_channel(:).');
            ADJ_SRC_elmt(idx_all) = w_channel .* residual;
            % ================================================================

            ADJ_SRC(:,:,elmt_idx) = ADJ_SRC_elmt;
        end

        % 可信度监控量
        phi_k = real(phi_k_accum) / max(tau_k_count,1);
        tau_k = tau_k_accum / max(tau_k_count,1);

        if ~isnan(prev_phi_k)
            delta_phi_k = prev_phi_k - phi_k;
        else
            delta_phi_k = nan;
        end

        if ~isnan(delta_phi_k) && ~isnan(prev_pred_drop) && abs(prev_pred_drop) > eps
            gamma_k = delta_phi_k / prev_pred_drop;
        else
            gamma_k = nan;
        end

        % Backproject error
        ADJ_WVFIELD = HS.solve(ADJ_SRC, true);
        SCALING     = repmat(reshape(scaling,[1,1,numel(scaling)]), [Nyi,Nxi,1]);
        BACKPROJ    = -real(conj(SCALING.*VIRT_SRC).*ADJ_WVFIELD);

        % ================================================================
        % [修改] 方向A3：Atten角度软权重已经在 ADJ_SRC 构建时逐 Rx 通道施加。
        % 因此这里的 BACKPROJ/gradient_img 已经是加权后的 Atten 梯度；
        % SoS 阶段则保持原版等权宽孔径梯度。
        % ================================================================

        for elmt_idx = 1:numel(tx_include_cur)
            gradient_img = gradient_img + BACKPROJ(:,:,elmt_idx);
        end

        % Remove Ringing from Data Gradient
        gradient_img = ringingRemovalFilt(xi_original, yi_original, ...
            gradient_img, c0, fDATA(f_idx), cutoff, ord);

        % ---- 三阶段低频先验正则化梯度项（方向E）----
        if enableLFPrior && ~updateAttenuation && ~useL2OriginalChain
            lf_ref_cur   = [];
            f_cutoff_cur = [];
            if stage2_ref_saved && fDATA(f_idx) > f_stage2_cutoff
                lf_ref_cur   = VEL_stage2_ref;
                f_cutoff_cur = f_stage2_cutoff;
            elseif stage1_ref_saved && fDATA(f_idx) > f_stage1_cutoff
                lf_ref_cur   = VEL_stage1_ref;
                f_cutoff_cur = f_stage1_cutoff;
            end

            if ~isempty(lf_ref_cur)
                lambda_k      = lambda_stage * (f_cutoff_cur / fDATA(f_idx))^2;
                lf_diff       = VEL_ESTIM - lf_ref_cur;
                grad_scale_lf = max(abs(gradient_img(:))) + eps;
                lf_diff_scale = max(abs(lf_diff(:))) + eps;
                lf_prior_reg  = -lambda_k * grad_scale_lf * (lf_diff / lf_diff_scale);
                gradient_img  = gradient_img + lf_prior_reg;
            end
        end

        % ---- 高频外环锚定（边界保护）----
        lambda_edge_eff = 0;
        edge_contam_k   = nan;
        edge_ref_cur    = [];
        if enableEdgeGuard && ~updateAttenuation && ~useL2OriginalChain && (fDATA(f_idx) >= f_edge_guard_start)
            if stage2_ref_saved
                edge_ref = VEL_stage2_ref;
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

        % ================================================================
        % [修改] 方向A3：Atten角度软权重已在 ADJ_SRC 构建阶段逐 Rx 精确施加。
        % 这里不再二次重算 Atten 梯度，避免重复加权并节省一次伴随求解。
        % ================================================================

        % 环伪影比监控
        [gx_ring, gy_ring] = gradient(gradient_img, dxi_original, dxi_original);
        R_pol  = sqrt(Xi.^2 + Yi.^2) + eps;
        er_x   = Xi ./ R_pol;
        er_y   = Yi ./ R_pol;
        grad_radial = gx_ring .* er_x + gy_ring .* er_y;
        grad_mag    = sqrt(gx_ring.^2 + gy_ring.^2);
        num_ring    = mean(abs(grad_radial(ring_bg_mask)), 'omitnan');
        den_ring    = mean(grad_mag(ring_bg_mask), 'omitnan') + eps;
        ring_k      = num_ring / den_ring;

        % Step 2: Compute New Conjugate Gradient Search Direction
        if useL2OriginalChain
            if do_reset_CG
                beta = 0;
            else
                betaPR = (gradient_img(:)'*(gradient_img(:)-gradient_img_prev(:))) / ...
                         (gradient_img_prev(:)'*gradient_img_prev(:));
                betaFR = (gradient_img(:)'*gradient_img(:)) / ...
                         (gradient_img_prev(:)'*gradient_img_prev(:));
                beta = min(max(betaPR,0), betaFR);
            end
        else
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
        end
        search_dir        = beta*search_dir - gradient_img;
        gradient_img_prev = gradient_img;

        % Step 3: Forward Projection of Current Search Direction
        PERTURBED_WVFIELD = HS.solve(VIRT_SRC.*search_dir, false);
        dREC_SIM = zeros(numel(tx_include_cur), numElements_cur);

        for elmt_idx = 1:numel(tx_include_cur)
            PERTURBED_WVFIELD_elmt = PERTURBED_WVFIELD(:,:,elmt_idx);
            mask_elmt = elemInclude_cur(tx_include_cur(elmt_idx),:);
            dREC_SIM(elmt_idx, mask_elmt) = -permute( ...
                scaling(elmt_idx) * PERTURBED_WVFIELD_elmt(ind_cur(mask_elmt)), [2,1]);
        end

        % Step 4: Line Search
        if useL2OriginalChain
            % ================================================================
            % [修改] fws: L2加权自洽线搜索
            % ================================================================
            % 由于 ADJ_SRC_elmt(idx_all) = w_channel .* residual，
            % L2梯度方向已经对应加权目标函数：
            %   J_L2 = 1/2 * sum(w_sr * |p_sim - p_obs|^2)
            % 因此线搜索分母也同步使用同一个 REC_WEIGHT_CUR：
            %   alpha_ls = -<g,d> / sum(w_sr * |Jd|^2)
            % 这样避免"梯度按加权目标函数、曲率按未加权目标函数"的不自洽。
            % 注意：w_sr 作为目标函数通道权重，只使用一次，不再平方。
            w_l2 = double(REC_WEIGHT_CUR);
            if isempty(w_l2) || all(w_l2(:) == 0)
                w_l2 = ones(size(dREC_SIM));
            end
            den_l2_weighted = sum(w_l2(:) .* abs(dREC_SIM(:)).^2);
            num_l2_weighted = -real(gradient_img(:)'*search_dir(:));
            alpha = num_l2_weighted / (den_l2_weighted + eps);

            if ~isfinite(alpha) || (alpha <= 0)
                alpha = 0;
            end
            % ================================================================
            alpha_ls = alpha;
            alpha_gate = 1.0;
            alpha_gate_xf = 1.0;
            C_lf_k = nan;
            mu_bg_k = nan;
            sigma_bg_k = nan;
            cos_sim_grad = nan;
            edge_damp_eff = 0;
            delta_s_max  = perc_step_size * alpha * max(abs(search_dir(:)));
            delta_v_max  = mean(VEL_ESTIM(:))^2 * delta_s_max;
            if updateAttenuation && enableAttenAngleWeight
                fprintf(['  [L2加权VE链 | %s | β=%.2f + A3] Φ=%.3e | ΔΦ=% .3e | τ=%.3eus | R=%.3f' ...
                         ' | alpha=%.2e | max|Δv|≈%.3f m/s | AttenW_min=%.2f | AttenW_mean=%.4f\n'], ...
                    stageLabelCur, beta_ve_cur, phi_k, delta_phi_k, tau_k*1e6, ring_k, alpha, delta_v_max, ...
                    atten_angle_weight_min_cur, atten_angle_weight_mean);
            else
                fprintf(['  [L2加权VE链 | %s | β=%.2f] Φ=%.3e | ΔΦ=% .3e | τ=%.3eus | R=%.3f' ...
                         ' | alpha=%.2e | max|Δv|≈%.3f m/s\n'], ...
                    stageLabelCur, beta_ve_cur, phi_k, delta_phi_k, tau_k*1e6, ring_k, alpha, delta_v_max);
            end
        else

        % 统一线搜索
        den_ls   = real(dREC_SIM(:)'*dREC_SIM(:));
        num_ls   = -real(gradient_img(:)'*search_dir(:));
        alpha_ls = num_ls / (den_ls + eps);

        % PolarPhase 自洽线搜索
        if strcmpi(misfitType_cur, 'PolarPhase')
            num_pp = 0;
            den_pp = 0;
            for elmt_idx = 1:numel(tx_include_cur)
                tx_id = tx_include_cur(elmt_idx);
                mask  = elemInclude_cur(tx_id, :);

                dpsi   = dREC_SIM(elmt_idx, mask).';
                phi_d  = phi_dir_store(elmt_idx, mask).';
                A_sim  = amp_sim_store(elmt_idx, mask).';
                dA     = amp_res_store(elmt_idx, mask).';
                dphi_v = dphi_store(elmt_idx, mask).';

                a = real(conj(phi_d) .* dpsi);
                b = imag(conj(phi_d) .* dpsi);

                if updateAttenuation
                    a_amp = a ./ (A_sim + eps);
                else
                    a_amp = a;
                end

                % fws: PolarPhase线搜索通道权重（与ADJ_SRC中的 w_channel 保持一致）
                %      SoS阶段：体现VE通道置信度 beta_ve_cur；
                %      Atten阶段：体现VE通道置信度 beta_ve_cur + A3角度软权重。
                %      注意：w_sr 作为目标函数通道权重，只使用一次，不再平方。
                w_ls = double(REC_WEIGHT_CUR(elmt_idx, mask).');
                if isempty(w_ls) || all(w_ls == 0)
                    w_ls = ones(size(a_amp));
                end

                num_pp = num_pp + sum(w_ls .* ( ...
                    (1 - alpha_hv_ls_cur)^2 .* dA     .* a_amp + ...
                     alpha_hv_ls_cur^2      .* dphi_v .* b ./ (A_sim + eps) ));

                den_pp = den_pp + sum(w_ls .* ( ...
                    (1 - alpha_hv_ls_cur)^2 .* a_amp.^2 + ...
                     alpha_hv_ls_cur^2      .* b.^2 ./ (A_sim.^2 + eps) ));
            end
            alpha_ls = -num_pp / (den_pp + eps);
        end

        % 步长处理（无 StepCap）
        % [修改] 按文本2逻辑删除 max|Δv| 人工夹紧：
        %        alpha 直接采用当前线搜索结果 alpha_ls；
        %        max|Δv| 仅用于日志诊断，不再反向限制 alpha。
        alpha = alpha_ls;

        if ~isfinite(alpha) || (alpha <= 0)
            alpha = 0;
        end

        % 跨频一致性度量
        C_lf_k        = nan;
        mu_bg_k       = nan;
        sigma_bg_k    = nan;
        cos_sim_grad  = nan;
        alpha_gate_xf = 1.0;
        if enableCrossFreqTrust && ~updateAttenuation
            if any(pure_water_mask(:))
                mu_bg_k    = mean(VEL_ESTIM(pure_water_mask), 'omitnan');
                sigma_bg_k = std(VEL_ESTIM(pure_water_mask), 0, 'omitnan');
            end

            if stage1_ref_saved
                if stage2_ref_saved && fDATA(f_idx) > f_stage2_cutoff
                    vel_ref_xf = VEL_stage2_ref;
                else
                    vel_ref_xf = VEL_stage1_ref;
                end
                sigma_lp = max(1.0, c0 / max(2*f_stage1_cutoff*dxi_original, eps));
                r_lp = max(2, ceil(3*sigma_lp));
                [xx_lp, yy_lp] = meshgrid(-r_lp:r_lp, -r_lp:r_lp);
                ker_lp = exp(-(xx_lp.^2 + yy_lp.^2) / (2*sigma_lp^2));
                ker_lp = ker_lp / (sum(ker_lp(:)) + eps);
                vel_lp_cur = conv2(VEL_ESTIM,  ker_lp, 'same');
                vel_lp_ref = conv2(vel_ref_xf, ker_lp, 'same');
                C_lf_k = norm(vel_lp_cur(:) - vel_lp_ref(:)) / (norm(vel_lp_ref(:)) + eps);
            end

            if exist('grad_stage1_ref', 'var')
                g_ref = grad_stage1_ref(:);
                g_cur = gradient_img(:);
                cos_sim_grad = real(g_ref' * g_cur) / ((norm(g_ref) + eps) * (norm(g_cur) + eps));
            end
        end

        % 可信度门控
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
        if enableCrossFreqTrust && ~updateAttenuation
            if ~isnan(C_lf_k)
                clf_rise = ~isnan(prev_clf_k) && (C_lf_k > clf_rise_ratio * prev_clf_k);
                if (C_lf_k > clf_warn) || clf_rise
                    alpha_gate_xf = min(alpha_gate_xf, ...
                        max(alpha_gate_xf_min, 1 - 5*max(C_lf_k - clf_warn, 0)));
                end
            end
            if ~isnan(sigma_bg_k) && (sigma_bg_k > sigma_bg_warn)
                alpha_gate_xf = min(alpha_gate_xf, ...
                    max(alpha_gate_xf_min, 1 - (sigma_bg_k - sigma_bg_warn)/20));
            end
            if ~isnan(cos_sim_grad) && (cos_sim_grad < cos_grad_warn)
                alpha_gate_xf = min(alpha_gate_xf, ...
                    max(alpha_gate_xf_min, (cos_sim_grad + 1)/2));
            end
        end
        alpha_gate = min(alpha_gate, alpha_gate_xf);
        alpha = alpha * alpha_gate;

        % 高频外环步长阻尼
        edge_damp_eff = 0;
        if enableEdgeGuard && ~updateAttenuation && (fDATA(f_idx) >= f_edge_guard_start)
            edge_damp_eff = edge_step_damp_max;
            if ~isnan(edge_contam_k) && ~isnan(prev_edge_contam) && ...
                    (edge_contam_k > edge_contam_rise_ratio * prev_edge_contam)
                edge_damp_eff = min(0.90, edge_damp_eff * 1.20);
            end
            search_dir = search_dir .* (1 - edge_damp_eff * edge_guard_mask);
        end

        % 调试打印
        delta_s_max  = perc_step_size * alpha * max(abs(search_dir(:)));
        delta_v_max  = mean(VEL_ESTIM(:))^2 * delta_s_max;
        if updateAttenuation && enableAttenAngleWeight
            fprintf(['  [%s] Φ=%.3e | ΔΦ=% .3e | τ=%.3eus | R=%.3f | γ=% .3f' ...
                     ' | C_lf=%.3f σbg=%.2f cosg=%.2f' ...
                     ' | EdgeContam=%.3f λedge=%.2e damp=%.2f' ...
                     ' | alpha_ls=%.2e' ...
                     ' | alpha=%.2e(gate=%.2f,xf=%.2f) | max|Δv|≈%.3f m/s' ...
                     ' | AttenW_min=%.2f | AttenW_mean=%.4f\n'], ...
                stageLabelCur, phi_k, delta_phi_k, tau_k*1e6, ring_k, gamma_k, ...
                C_lf_k, sigma_bg_k, cos_sim_grad, edge_contam_k, ...
                lambda_edge_eff, edge_damp_eff, alpha_ls, ...
                alpha, alpha_gate, alpha_gate_xf, delta_v_max, ...
                atten_angle_weight_min_cur, atten_angle_weight_mean);
        else
            fprintf(['  [%s] Φ=%.3e | ΔΦ=% .3e | τ=%.3eus | R=%.3f | γ=% .3f' ...
                     ' | C_lf=%.3f σbg=%.2f cosg=%.2f' ...
                     ' | EdgeContam=%.3f λedge=%.2e damp=%.2f' ...
                     ' | alpha_ls=%.2e' ...
                     ' | alpha=%.2e(gate=%.2f,xf=%.2f) | max|Δv|≈%.3f m/s\n'], ...
                stageLabelCur, phi_k, delta_phi_k, tau_k*1e6, ring_k, gamma_k, ...
                C_lf_k, sigma_bg_k, cos_sim_grad, edge_contam_k, ...
                lambda_edge_eff, edge_damp_eff, alpha_ls, ...
                alpha, alpha_gate, alpha_gate_xf, delta_v_max);
        end

        end

        % 本步预测下降
        pred_drop_cur = -perc_step_size * alpha * real(gradient_img(:)'*search_dir(:));

        % Update slowness
        if updateAttenuation
            SI = sign(sign_conv) * imag(SLOW_ESTIM) + perc_step_size * alpha * search_dir;
            SLOW_ESTIM = real(SLOW_ESTIM) + 1i * sign(sign_conv) * SI;
        else
            SLOW_ESTIM = SLOW_ESTIM + perc_step_size * alpha * search_dir;
        end
        VEL_ESTIM   = 1./real(SLOW_ESTIM);

        % 外环后处理回拉（高频SoS阶段）
        edge_blend_eff = 0;
        if enableEdgeGuard && ~updateAttenuation && ~useL2OriginalChain && (fDATA(f_idx) >= f_edge_guard_start) && ~isempty(edge_ref_cur)
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
        LAMBDA_EDGE_EFF_ITER(iter) = lambda_edge_eff;
        EDGE_DAMP_EFF_ITER(iter)   = edge_damp_eff;
        EDGE_CONTAM_ITER(iter)     = edge_contam_k;
        EDGE_BLEND_EFF_ITER(iter)  = edge_blend_eff;
        C_LF_ITER(iter)            = C_lf_k;
        MU_BG_ITER(iter)           = mu_bg_k;
        SIGMA_BG_ITER(iter)        = sigma_bg_k;
        COS_GRAD_ITER(iter)        = cos_sim_grad;
        ALPHA_GATE_XF_ITER(iter)   = alpha_gate_xf;
        BETA_VE_ITER(iter)         = beta_ve_cur;
        BETA_VE_STAGE_IDX_ITER(iter) = beta_ve_stage_idx_cur;
        FRAC_SHIFT_ITER(iter)        = frac_shift_cur;
        FRAC_SHIFT_STAGE_IDX_ITER(iter) = frac_shift_stage_idx_cur;
        FRAC_EXCL_ITER(iter)         = frac_excl_cur;
        FRAC_EXCL_STAGE_IDX_ITER(iter) = frac_excl_stage_idx_excl_cur;
        NUM_RX_EXCL_SIDE_ITER(iter)  = numElemLeftRightExcl_cur;
        NUM_RX_KEEP_ITER(iter)       = numRxKeep_cur;
        ATTEN_ANGLE_WEIGHT_ITER(iter)= atten_angle_weight_mean;
        ATTEN_ANGLE_WEIGHT_MIN_ITER(iter)= atten_angle_weight_min_cur;
        MISFIT_OVERRIDE_ITER(iter)         = isMisfitOverrideCur;
        MISFIT_OVERRIDE_STAGE_ITER{iter}   = misfitOverrideStageCur;
        MISFIT_OVERRIDE_MODE_ITER{iter}    = misfitOverrideModeCur;
        MISFIT_OVERRIDE_REASON_ITER{iter}  = misfitOverrideReasonCur;
        MISFIT_TYPE_ITER{iter}             = misfitType_cur;
        savedIters(end+1)          = iter; %#ok<SAGROW>

        % 迭代状态滚动
        prev_phi_k       = phi_k;
        prev_tau_k       = tau_k;
        prev_ring_k      = ring_k;
        prev_clf_k       = C_lf_k;
        prev_pred_drop   = pred_drop_cur;
        prev_edge_contam = edge_contam_k;

        % Visualize coarse
        hfig = figure(1);
        set(hfig, 'Name', sprintf('f=%.3fMHz | %s | iter=%d/%d | %s | misfit=%s', ...
            fDATA(f_idx)/1e6, cur_fs_tag, iter, Niter, stageLabelCur, misfitLabelCur));

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

        % GIF 帧保存
        if saveGIF
            gif_iter_count = gif_iter_count + 1;
            if mod(gif_iter_count - 1, gif_save_every) == 0
                set(hgif_vel, 'CData', VEL_ESTIM);
                set(hgif_att, 'CData', Np2dB*slow2atten*ATTEN_ESTIM);

                title(axgif_vel, sprintf('SoS | iter=%d  f=%.3fMHz  %s  %s  misfit=%s', ...
                    iter, fDATA(f_idx)/1e6, cur_fs_tag, stageLabelCur, misfitLabelCur), ...
                    'FontSize', 10, 'Interpreter', 'none');

                if updateAttenuation
                    title(axgif_att, sprintf('Attenuation | iter=%d  %s  misfit=%s  β_ve=%.2f  AttenW=%.3f', ...
                        iter, stageLabelCur, misfitLabelCur, beta_ve_cur, atten_angle_weight_mean), ...
                        'FontSize', 10, 'Interpreter', 'none');
                else
                    title(axgif_att, sprintf('Attenuation | iter=%d  %s  not updated', ...
                        iter, stageLabelCur), ...
                        'FontSize', 10, 'Interpreter', 'none');
                end
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

        disp(['Iteration ', num2str(iter), ...
              ' | f=', num2str(fDATA(f_idx)/1e6,'%.3f'), 'MHz', ...
              ' | ', cur_fs_tag, ...
              ' | keepRx=', num2str(numRxKeep_cur), '/', num2str(numElements_cur), ...
              ' | stage=', stageLabelCur, ...
              ' | Atten=', num2str(updateAttenuation), ...
              ' | misfit=', misfitLabelCur]);
        toc;
    end

    % 邻频高斯平滑衔接：仅在三阶段边界触发
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

        R_trans = max(sqrt(transducerPositionsX.^2 + transducerPositionsY.^2));
        smooth_mask = (R_grid <= R_trans * 0.98);

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


%% ===================== Metrics ROI Definition =====================
enableMetricsROIMask  = true;
metrics_roi_radius_mode = 'median';
metrics_roi_ratio_main  = 0.85;

switch lower(metrics_roi_radius_mode)
    case 'median'
        R_trans_metrics = median(sqrt(transducerPositionsX.^2 + transducerPositionsY.^2), 'omitnan');
    case 'max'
        R_trans_metrics = max(sqrt(transducerPositionsX.^2 + transducerPositionsY.^2));
    otherwise
        error('Unknown metrics_roi_radius_mode: %s.', metrics_roi_radius_mode);
end
if ~isfinite(R_trans_metrics) || R_trans_metrics <= 0
    error('R_trans_metrics 非法，请检查 transducerPositionsX/Y。');
end

if enableMetricsROIMask
    roi_mask_metrics = (R_grid <= R_trans_metrics * metrics_roi_ratio_main);
else
    roi_mask_metrics = true(size(R_grid));
end

roi_mask = roi_mask_metrics;
R_trans  = R_trans_metrics;

roi_radius_metrics = R_trans_metrics * metrics_roi_ratio_main;
theta_roi_metrics  = linspace(0, 2*pi, 720);
x_roi_metrics      = roi_radius_metrics * cos(theta_roi_metrics);
y_roi_metrics      = roi_radius_metrics * sin(theta_roi_metrics);

metrics_roi_info = struct();
metrics_roi_info.enableMetricsROIMask   = enableMetricsROIMask;
metrics_roi_info.radius_mode            = metrics_roi_radius_mode;
metrics_roi_info.R_trans_m              = R_trans_metrics;
metrics_roi_info.roi_ratio_main         = metrics_roi_ratio_main;
metrics_roi_info.num_roi_pixels_main    = nnz(roi_mask_metrics);
metrics_roi_info.num_grid_pixels        = numel(R_grid);


%% ===================== X/Y Cutting-Line Layout + Main Profiles =====================
cut_x_m = -0.01;
cut_y_m = 0.01;
showStageRefsInMainProfile = false;
showStageRefsStandalone    = true;
lineColorX = [1 0 0];
lineColorY = [0 0 1];

[~, ix_coarse] = min(abs(xi_original - cut_x_m));
[~, ix_true]   = min(abs(xi_orig     - cut_x_m));
[~, iy_coarse] = min(abs(yi_original - cut_y_m));
[~, iy_true]   = min(abs(yi_orig     - cut_y_m));

vel_prof_x_est  = VEL_ESTIM(iy_coarse, :);
vel_prof_x_true = C(iy_true, :);
vel_prof_y_est  = VEL_ESTIM(:, ix_coarse);
vel_prof_y_true = C(:, ix_true);

has_stage1 = exist('VEL_stage1_ref', 'var') == 1;
has_stage2 = exist('VEL_stage2_ref', 'var') == 1;

if has_stage1
    vel_prof_x_stage1 = VEL_stage1_ref(iy_coarse, :);
    vel_prof_y_stage1 = VEL_stage1_ref(:, ix_coarse);
end
if has_stage2
    vel_prof_x_stage2 = VEL_stage2_ref(iy_coarse, :);
    vel_prof_y_stage2 = VEL_stage2_ref(:, ix_coarse);
end

figure(5); clf;
set(gcf, 'Color', 'w', 'Position', [80 80 1300 760]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t, [2 1]);
imagesc(ax1, xi_original, yi_original, VEL_ESTIM, crange);
axis(ax1, 'image'); colorbar(ax1); colormap(ax1, cmap_rb);
hold(ax1, 'on');
plot(ax1, [xi_original(1), xi_original(end)], ...
          [yi_original(iy_coarse), yi_original(iy_coarse)], ...
          '-', 'Color', lineColorX, 'LineWidth', 1.8, ...
          'DisplayName', sprintf('X-profile line: y = %.1f mm', yi_original(iy_coarse)*1e3));
plot(ax1, [xi_original(ix_coarse), xi_original(ix_coarse)], ...
          [yi_original(1), yi_original(end)], ...
          '-', 'Color', lineColorY, 'LineWidth', 1.8, ...
          'DisplayName', sprintf('Y-profile line: x = %.1f mm', xi_original(ix_coarse)*1e3));
plot(ax1, x_roi_metrics, y_roi_metrics, 'k--', 'LineWidth', 1.8, ...
          'DisplayName', sprintf('Metrics ROI: R <= %.1f mm', roi_radius_metrics*1e3));
hold(ax1, 'off');
xlabel(ax1, 'Lateral [m]'); ylabel(ax1, 'Axial [m]');
title(ax1, 'Estimated SOS with profile lines and metrics ROI');
legend(ax1, 'Location', 'best');

ax2 = nexttile(t, 2);
plot(ax2, xi_orig, vel_prof_x_true, 'k-', 'LineWidth', 2.0, 'DisplayName', 'True'); hold(ax2, 'on');
plot(ax2, xi_original, vel_prof_x_est, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Estimated');
if showStageRefsInMainProfile
    if has_stage1, plot(ax2, xi_original, vel_prof_x_stage1, 'b-.', 'LineWidth', 1.3, 'DisplayName', 'Stage1 ref'); end
    if has_stage2, plot(ax2, xi_original, vel_prof_x_stage2, 'm-.', 'LineWidth', 1.3, 'DisplayName', 'Stage2 ref'); end
end
grid(ax2, 'on');
xlabel(ax2, 'Lateral position [m]'); ylabel(ax2, 'Sound Speed [m/s]');
title(ax2, sprintf('X-direction SOS profile (y = %.1f mm)', yi_original(iy_coarse)*1e3));
legend(ax2, 'Location', 'best');
xlim(ax2, [xi_original(1), xi_original(end)]); ylim(ax2, crange);

ax3 = nexttile(t, 4);
plot(ax3, vel_prof_y_true, yi_orig, 'k-', 'LineWidth', 2.0, 'DisplayName', 'True'); hold(ax3, 'on');
plot(ax3, vel_prof_y_est, yi_original, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Estimated');
if showStageRefsInMainProfile
    if has_stage1, plot(ax3, vel_prof_y_stage1, yi_original, 'b-.', 'LineWidth', 1.3, 'DisplayName', 'Stage1 ref'); end
    if has_stage2, plot(ax3, vel_prof_y_stage2, yi_original, 'm-.', 'LineWidth', 1.3, 'DisplayName', 'Stage2 ref'); end
end
set(ax3, 'YDir', 'reverse'); grid(ax3, 'on');
xlabel(ax3, 'Sound Speed [m/s]'); ylabel(ax3, 'Axial position [m]');
title(ax3, sprintf('Y-direction SOS profile (x = %.1f mm, iter=%d)', xi_original(ix_coarse)*1e3, iter));
legend(ax3, 'Location', 'best');
xlim(ax3, crange); ylim(ax3, [yi_original(1), yi_original(end)]);

if showStageRefsStandalone
    if has_stage1 || has_stage2
        figure(6); clf;
        set(gcf, 'Color', 'w', 'Position', [120 120 1250 520]);
        t_stage = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        ax_stage_x = nexttile(t_stage, 1);
        plot(ax_stage_x, xi_orig, vel_prof_x_true, 'k-', 'LineWidth', 2.0, 'DisplayName', 'True'); hold(ax_stage_x, 'on');
        plot(ax_stage_x, xi_original, vel_prof_x_est, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Estimated');
        if has_stage1, plot(ax_stage_x, xi_original, vel_prof_x_stage1, 'b-.', 'LineWidth', 1.4, 'DisplayName', 'Stage1 ref'); end
        if has_stage2, plot(ax_stage_x, xi_original, vel_prof_x_stage2, 'm-.', 'LineWidth', 1.4, 'DisplayName', 'Stage2 ref'); end
        grid(ax_stage_x, 'on');
        xlabel(ax_stage_x, 'Lateral position [m]'); ylabel(ax_stage_x, 'Sound Speed [m/s]');
        title(ax_stage_x, sprintf('Stage-reference X profile (y = %.1f mm)', yi_original(iy_coarse)*1e3));
        legend(ax_stage_x, 'Location', 'best');
        xlim(ax_stage_x, [xi_original(1), xi_original(end)]); ylim(ax_stage_x, crange);

        ax_stage_y = nexttile(t_stage, 2);
        plot(ax_stage_y, vel_prof_y_true, yi_orig, 'k-', 'LineWidth', 2.0, 'DisplayName', 'True'); hold(ax_stage_y, 'on');
        plot(ax_stage_y, vel_prof_y_est, yi_original, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Estimated');
        if has_stage1, plot(ax_stage_y, vel_prof_y_stage1, yi_original, 'b-.', 'LineWidth', 1.4, 'DisplayName', 'Stage1 ref'); end
        if has_stage2, plot(ax_stage_y, vel_prof_y_stage2, yi_original, 'm-.', 'LineWidth', 1.4, 'DisplayName', 'Stage2 ref'); end
        set(ax_stage_y, 'YDir', 'reverse'); grid(ax_stage_y, 'on');
        xlabel(ax_stage_y, 'Sound Speed [m/s]'); ylabel(ax_stage_y, 'Axial position [m]');
        title(ax_stage_y, sprintf('Stage-reference Y profile (x = %.1f mm)', xi_original(ix_coarse)*1e3));
        legend(ax_stage_y, 'Location', 'best');
        xlim(ax_stage_y, crange); ylim(ax_stage_y, [yi_original(1), yi_original(end)]);

        title(t_stage, sprintf('Two-stage SOS profile comparison | iter=%d', iter), 'Interpreter', 'none');
    else
        fprintf('[Profile] 未检测到 VEL_stage1_ref / VEL_stage2_ref，跳过figure(6)两阶段参考profile。\n');
    end
end

scriptElapsedSec = toc(scriptTimer);
fprintf('整个脚本总时长: %.2f 秒 (%.2f 分钟)\n', scriptElapsedSec, scriptElapsedSec/60);

%% ===================== Quantitative Metrics (SoS / Attenuation) =====================
targetDevPct     = 85;
bgDevPct         = 40;
bg_ref_tol       = 5.0;    % [m/s]
bg_ref_tol_atten = 0.05;   % [dB/(cm MHz)]

if ~exist('roi_mask_metrics', 'var') || isempty(roi_mask_metrics)
    error('roi_mask_metrics 未定义，请检查 Metrics ROI Definition 代码块是否被执行。');
end
if ~exist('metrics_roi_info', 'var') || isempty(metrics_roi_info)
    metrics_roi_info = struct();
end
metrics_roi_info.num_roi_pixels_main = nnz(roi_mask_metrics);
metrics_roi_info.num_grid_pixels     = numel(R_grid);

fprintf('\n[METRICS ROI] enable=%d | R_trans=%.3f mm (%s) | ratio_main=%.2f | ROI pixels=%d/%d\n', ...
    enableMetricsROIMask, R_trans_metrics*1e3, metrics_roi_radius_mode, ...
    metrics_roi_ratio_main, nnz(roi_mask_metrics), numel(R_grid));

[X_true, Y_true]     = meshgrid(xi_orig, yi_orig);
[X_recon, Y_recon]   = meshgrid(xi_original, yi_original);
C_true_on_reconGrid  = interp2(X_true, Y_true, C,     X_recon, Y_recon, 'linear', nan);
atten_true_on_reconGrid = interp2(X_true, Y_true, atten, X_recon, Y_recon, 'linear', nan);

ATTEN_ESTIM_dB = Np2dB * slow2atten * ATTEN_ESTIM;

[PSNR_SoS_dB, RMSE_SoS, SSIM_SoS, RD_SoS_percent, ...
    target_mask_SoS, bg_mask_SoS, valid_mask_SoS, metrics_SoS] = ...
    computeReconMetrics(VEL_ESTIM, C_true_on_reconGrid, ...
                        targetDevPct, bgDevPct, bg_ref_tol, 'SoS', roi_mask_metrics);

[PSNR_Att_dB, RMSE_Att, SSIM_Att, RD_Att_percent, ...
    target_mask_Att, bg_mask_Att, valid_mask_Att, metrics_Att] = ...
    computeReconMetrics(ATTEN_ESTIM_dB, atten_true_on_reconGrid, ...
                        targetDevPct, bgDevPct, bg_ref_tol_atten, 'Attenuation', roi_mask_metrics);

fprintf('\n================ 定量指标 (SoS) ================\n');
fprintf('PSNR : %.4f dB\n', PSNR_SoS_dB);
fprintf('RMSE : %.6f m/s\n', RMSE_SoS);
fprintf('SSIM : %.6f\n', SSIM_SoS);
fprintf('RD   : %.4f %%\n', RD_SoS_percent);
fprintf('目标像素数: %d, 背景像素数: %d\n', nnz(target_mask_SoS), nnz(bg_mask_SoS));
fprintf('===============================================\n\n');

fprintf('\n============ 定量指标 (Attenuation) ============\n');
fprintf('PSNR : %.4f dB\n', PSNR_Att_dB);
fprintf('RMSE : %.6f dB/(cm MHz)\n', RMSE_Att);
fprintf('SSIM : %.6f\n', SSIM_Att);
fprintf('RD   : %.4f %%\n', RD_Att_percent);
fprintf('目标像素数: %d, 背景像素数: %d\n', nnz(target_mask_Att), nnz(bg_mask_Att));
fprintf('===============================================\n\n');

%% ===================== Save =====================
suffix = 'WaveformInversionResults';
filename_results = [result_dir, filename, '_', verTag, '_', suffix, '.mat'];

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
    LAMBDA_EDGE_EFF_ITER  = LAMBDA_EDGE_EFF_ITER(savedIters);
    EDGE_DAMP_EFF_ITER    = EDGE_DAMP_EFF_ITER(savedIters);
    EDGE_CONTAM_ITER      = EDGE_CONTAM_ITER(savedIters);
    EDGE_BLEND_EFF_ITER   = EDGE_BLEND_EFF_ITER(savedIters);
    C_LF_ITER             = C_LF_ITER(savedIters);
    MU_BG_ITER            = MU_BG_ITER(savedIters);
    SIGMA_BG_ITER         = SIGMA_BG_ITER(savedIters);
    COS_GRAD_ITER         = COS_GRAD_ITER(savedIters);
    ALPHA_GATE_XF_ITER    = ALPHA_GATE_XF_ITER(savedIters);
    BETA_VE_ITER               = BETA_VE_ITER(savedIters);
    BETA_VE_STAGE_IDX_ITER     = BETA_VE_STAGE_IDX_ITER(savedIters);
    FRAC_SHIFT_ITER            = FRAC_SHIFT_ITER(savedIters);
    FRAC_SHIFT_STAGE_IDX_ITER  = FRAC_SHIFT_STAGE_IDX_ITER(savedIters);
    FRAC_EXCL_ITER             = FRAC_EXCL_ITER(savedIters);
    FRAC_EXCL_STAGE_IDX_ITER   = FRAC_EXCL_STAGE_IDX_ITER(savedIters);
    NUM_RX_EXCL_SIDE_ITER      = NUM_RX_EXCL_SIDE_ITER(savedIters);
    NUM_RX_KEEP_ITER           = NUM_RX_KEEP_ITER(savedIters);
    ATTEN_ANGLE_WEIGHT_ITER    = ATTEN_ANGLE_WEIGHT_ITER(savedIters);
    ATTEN_ANGLE_WEIGHT_MIN_ITER = ATTEN_ANGLE_WEIGHT_MIN_ITER(savedIters);
    MISFIT_OVERRIDE_ITER            = MISFIT_OVERRIDE_ITER(savedIters);
    MISFIT_OVERRIDE_STAGE_ITER      = MISFIT_OVERRIDE_STAGE_ITER(savedIters);
    MISFIT_OVERRIDE_MODE_ITER       = MISFIT_OVERRIDE_MODE_ITER(savedIters);
    MISFIT_OVERRIDE_REASON_ITER     = MISFIT_OVERRIDE_REASON_ITER(savedIters);
    MISFIT_TYPE_ITER                = MISFIT_TYPE_ITER(savedIters);
end

xi = xi_original; yi = yi_original;
VEL_ESTIM_FINAL_PPSTAGE2ATTEN   = VEL_ESTIM;
ATTEN_ESTIM_FINAL_PPSTAGE2ATTEN = ATTEN_ESTIM;

save(filename_results, '-v7.3', ...
    'xi','yi','xi_original','yi_original', ...
    'VEL_ESTIM_FINAL_PPSTAGE2ATTEN','ATTEN_ESTIM_FINAL_PPSTAGE2ATTEN', ...
    'fDATA','niterAttenPerFreq','niterSoSPerFreq', ...
    'VEL_ESTIM_ITER','ATTEN_ESTIM_ITER','GRAD_IMG_ITER','SEARCH_DIR_ITER', ...
    'PHI_ITER','DELTA_PHI_ITER','TAU_SHIFT_ITER','R_RING_ITER', ...
    'GAMMA_TRUST_ITER','ALPHA_GATE_ITER', ...
    'LAMBDA_EDGE_EFF_ITER','EDGE_DAMP_EFF_ITER', ...
    'EDGE_CONTAM_ITER','EDGE_BLEND_EFF_ITER', ...
    'C_LF_ITER','MU_BG_ITER','SIGMA_BG_ITER','COS_GRAD_ITER','ALPHA_GATE_XF_ITER', ...
    'savedIters', 'enableVirtualElements', 'useVirtualTx', 'orig_numElements', ...
    'virtualElementsMode', 'exclMode','numElemLeftRightExcl','frac_excl', ...
    'enableFracExclStageSchedule','frac_excl_schedule','stage_freq_boundary','frac_excl_precompute', ...
    'FRAC_EXCL_ITER','FRAC_EXCL_STAGE_IDX_ITER','NUM_RX_EXCL_SIDE_ITER','NUM_RX_KEEP_ITER', ...
    'VE_freqMode','VE_freqThresh','VE_manualFlags','VE_flags', ...
    'frac_shift_sos_lo','frac_shift_sos_hi','frac_shift_sos_f_split', ...
    'beta_ve_edges','beta_ve_sos_schedule','beta_ve_atten_schedule', ...
    'BETA_VE_ITER','BETA_VE_STAGE_IDX_ITER','FRAC_SHIFT_ITER','FRAC_SHIFT_STAGE_IDX_ITER', ...
    'enableAttenAngleWeight','attenAngleWeightMode','attenAngleWeightExp','attenAngleWeightMinSchedule','ATTEN_ANGLE_WEIGHT_ITER','ATTEN_ANGLE_WEIGHT_MIN_ITER', ...
    'misfitStage_SoSOnly','misfitStage_SoSAtten_SoS','misfitStage_SoSAtten_Att', ...
    'enableMisfitFreqSchedule_SoSOnly','misfitFreqRanges_SoSOnly','misfitFreqPoints_SoSOnly', ...
    'enableMisfitFreqSchedule_SoSAtten_SoS','misfitFreqRanges_SoSAtten_SoS','misfitFreqPoints_SoSAtten_SoS', ...
    'enableMisfitFreqSchedule_SoSAtten_Att','misfitFreqRanges_SoSAtten_Att','misfitFreqPoints_SoSAtten_Att', ...
    'misfitFreqTol_Hz','MISFIT_OVERRIDE_ITER','MISFIT_OVERRIDE_STAGE_ITER','MISFIT_OVERRIDE_MODE_ITER','MISFIT_OVERRIDE_REASON_ITER','MISFIT_TYPE_ITER', ...
    'alpha_hv_src', 'alpha_hv_ls', 'alpha_hv_atten_src', 'alpha_hv_atten_ls', 'polarLogAmpFloorRel', ...
    'useOriginalVEChainForL2', ...
    'enableLFPrior', 'lambda_stage', 'f_stage1_cutoff', 'f_stage2_cutoff', ...
    'stage1_ref_saved', 'stage2_ref_saved', ...
    'enableTrustGate','phi_small_drop_ratio','tau_rise_ratio','ring_rise_ratio', ...
    'enableCrossFreqTrust','clf_warn','clf_rise_ratio','sigma_bg_warn', ...
    'cos_grad_warn','alpha_gate_xf_min','water_bg_inner_ratio','water_bg_outer_ratio', ...
    'enableEdgeGuard','f_edge_guard_start', ...
    'edge_r_inner_ratio','edge_r_outer_ratio', ...
    'lambda_edge_anchor','edge_step_damp_max', ...
    'edge_blend_base','edge_blend_max','edge_contam_rise_ratio', ...
    'ring_bg_inner_ratio','ring_bg_outer_ratio', ...
    'mainLoopElapsedSec', 'scriptElapsedSec', 'gif_filepath', ...
    'PSNR_SoS_dB','RMSE_SoS','SSIM_SoS','RD_SoS_percent', ...
    'PSNR_Att_dB','RMSE_Att','SSIM_Att','RD_Att_percent', ...
    'metrics_SoS','metrics_Att', ...
    'enableMetricsROIMask','metrics_roi_radius_mode','metrics_roi_ratio_main', ...
    'R_trans','R_trans_metrics','roi_mask','roi_mask_metrics','metrics_roi_info', ...
    'roi_radius_metrics','x_roi_metrics','y_roi_metrics', ...
    'cut_x_m','cut_y_m','showStageRefsInMainProfile','showStageRefsStandalone', ...
    'C_true_on_reconGrid','atten_true_on_reconGrid','ATTEN_ESTIM_dB', ...
    'target_mask_SoS','bg_mask_SoS','valid_mask_SoS', ...
    'target_mask_Att','bg_mask_Att','valid_mask_Att', ...
    'targetDevPct','bgDevPct','bg_ref_tol','bg_ref_tol_atten');

if exist('VEL_stage1_ref', 'var')
    save(filename_results, '-append', 'VEL_stage1_ref');
end
if exist('grad_stage1_ref', 'var')
    save(filename_results, '-append', 'grad_stage1_ref');
end
if exist('VEL_stage2_ref', 'var')
    save(filename_results, '-append', 'VEL_stage2_ref');
end

disp(['Save completed: ', filename_results]);


%% ============================ Local Functions ============================
function [frac_excl_cur, stage_idx, stage_tag] = getFracExclByFreq(freqHz, enableSchedule, frac_excl_default, frac_excl_schedule, stage_freq_boundary)
% 根据当前频率返回动态 frac_excl
% SoS 与 Atten 共用同一套调度，保持目标函数一致性
    if enableSchedule
        if freqHz <= stage_freq_boundary(1)
            stage_idx = 1;
            stage_tag = 'Stage1';
        elseif freqHz <= stage_freq_boundary(2)
            stage_idx = 2;
            stage_tag = 'Stage2';
        else
            stage_idx = 3;
            stage_tag = 'Stage3';
        end
        frac_excl_cur = frac_excl_schedule(stage_idx);
    else
        stage_idx = 0;
        stage_tag = 'fixed';
        frac_excl_cur = frac_excl_default;
    end
end

function [elemInclude, numElemLeftRightExcl_eff, keepRx_eff] = buildNearEndExclusionMask(numElements_eff, exclMode, numElemLeftRightExcl, frac_excl_eff, theta_excl, theta_all)
% 统一生成 Near-end exclusion 掩码
% SoS/Atten 共用同一套掩码（保持目标函数一致性）
    elemInclude = true(numElements_eff, numElements_eff);
    numElemLeftRightExcl_eff = nan;

    switch lower(exclMode)
        case 'manual'
            numElemLeftRightExcl_eff = numElemLeftRightExcl;
            elemLeftRightExcl_eff = -numElemLeftRightExcl_eff:numElemLeftRightExcl_eff;
            for tx_element = 1:numElements_eff
                elemLeftRightExclCurrent = elemLeftRightExcl_eff + tx_element;
                elemLeftRightExclCurrent(elemLeftRightExclCurrent<1) = numElements_eff + ...
                    elemLeftRightExclCurrent(elemLeftRightExclCurrent<1);
                elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements_eff) = ...
                    elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements_eff) - numElements_eff;
                elemInclude(tx_element, elemLeftRightExclCurrent) = false;
            end

        case 'angle'
            for tx_element = 1:numElements_eff
                dtheta = theta_all - theta_all(tx_element);
                dtheta = atan2(sin(dtheta), cos(dtheta));
                elemInclude(tx_element, abs(dtheta) <= theta_excl) = false;
            end

        case 'frac'
            numElemLeftRightExcl_eff = round(numElements_eff * frac_excl_eff / 2);
            elemLeftRightExcl_eff = -numElemLeftRightExcl_eff:numElemLeftRightExcl_eff;
            for tx_element = 1:numElements_eff
                elemLeftRightExclCurrent = elemLeftRightExcl_eff + tx_element;
                elemLeftRightExclCurrent(elemLeftRightExclCurrent<1) = numElements_eff + ...
                    elemLeftRightExclCurrent(elemLeftRightExclCurrent<1);
                elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements_eff) = ...
                    elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements_eff) - numElements_eff;
                elemInclude(tx_element, elemLeftRightExclCurrent) = false;
            end

        otherwise
            error('Unknown exclMode: %s. Use ''manual'', ''angle'' or ''frac''.', exclMode);
    end

    keepRx_eff = min(sum(elemInclude, 2));
end

function validateMisfitScheduleTable(ruleTable, tableName, isRangeTable)
% 校验 misfit 频率调度表
    if isempty(ruleTable)
        return;
    end
    if ~iscell(ruleTable) || size(ruleTable,2) ~= 2
        error('%s 必须是 N×2 cell：第1列为频率点/频率段，第2列为 ''L2'' 或 ''PolarPhase''。', tableName);
    end
    for rr = 1:size(ruleTable,1)
        freqSpec = ruleTable{rr,1};
        misfitName = ruleTable{rr,2};
        if ~isnumeric(freqSpec) || isempty(freqSpec) || any(~isfinite(freqSpec(:)))
            error('%s 第 %d 行频率必须是有限数值，单位 Hz。', tableName, rr);
        end
        if isRangeTable
            if numel(freqSpec) ~= 2
                error('%s 第 %d 行是频率段规则，必须写成 [f_min, f_max]。', tableName, rr);
            end
            if freqSpec(1) > freqSpec(2)
                error('%s 第 %d 行 f_min 不能大于 f_max。', tableName, rr);
            end
        end
        validateMisfitName(misfitName, sprintf('%s 第 %d 行', tableName, rr));
    end
end

function validateMisfitName(misfitName, whereStr)
    if ~(ischar(misfitName) || isstring(misfitName))
        error('%s 的方案名必须是字符：''L2'' 或 ''PolarPhase''。', whereStr);
    end
    misfitName = char(misfitName);
    if ~(strcmpi(misfitName, 'L2') || strcmpi(misfitName, 'PolarPhase'))
        error('%s 的方案名为 %s，当前仅支持 ''L2'' 或 ''PolarPhase''。', whereStr, misfitName);
    end
end

function [misfitCur, isOverride, overrideMode, overrideReason] = applyMisfitFreqSchedule(freqHz, defaultMisfit, enableSchedule, rangeTable, pointTable, tolHz)
% 根据"频率段 + 频率点"调度表，返回当前频率最终使用的 misfit
% 优先级：默认 misfitStage_* < 频率段表 < 频率点表
    misfitCur = char(defaultMisfit);
    isOverride = false;
    overrideMode = 'default';
    overrideReason = 'default';

    if ~enableSchedule
        return;
    end
    if nargin < 6 || isempty(tolHz)
        tolHz = 1.0;
    end

    if ~isempty(rangeTable)
        for rr = 1:size(rangeTable,1)
            fRange = rangeTable{rr,1};
            if freqHz >= fRange(1) - tolHz && freqHz <= fRange(2) + tolHz
                misfitCur = char(rangeTable{rr,2});
                isOverride = true;
                overrideMode = 'range';
                overrideReason = sprintf('range[%.3f,%.3f]MHz', fRange(1)/1e6, fRange(2)/1e6);
            end
        end
    end

    if ~isempty(pointTable)
        for rr = 1:size(pointTable,1)
            fPoints = pointTable{rr,1};
            [minDist, hitIdx] = min(abs(freqHz - fPoints(:)));
            if minDist <= tolHz
                hitPoint = fPoints(hitIdx);
                misfitCur = char(pointTable{rr,2});
                isOverride = true;
                overrideMode = 'point';
                overrideReason = sprintf('point %.3fMHz', hitPoint/1e6);
            end
        end
    end
end

function txt = formatMisfitPreview(misfitCur, isOverride, overrideMode, overrideReason)
    if isOverride
        txt = sprintf('%s override[%s:%s]', char(misfitCur), overrideMode, overrideReason);
    else
        txt = char(misfitCur);
    end
end


function [PSNR_dB, RMSE_val, SSIM_val, RD_percent, target_mask, bg_mask, valid_mask, metrics] = computeReconMetrics(recon_img, true_img, targetDevPct, bgDevPct, bg_ref_tol, metricName, roi_mask)
    % 通用定量指标计算函数（SoS / Attenuation 共用）
    if nargin < 7 || isempty(roi_mask)
        roi_mask = true(size(true_img));
    end
    if ~isequal(size(recon_img), size(true_img))
        error('computeReconMetrics: recon_img 与 true_img 尺寸必须一致。');
    end
    if ~isequal(size(roi_mask), size(true_img))
        error('computeReconMetrics: roi_mask 尺寸必须与 true_img / recon_img 一致。');
    end

    roi_mask = logical(roi_mask);
    valid_mask = isfinite(true_img) & isfinite(recon_img) & roi_mask;
    x_true  = true_img(valid_mask);
    x_recon = recon_img(valid_mask);

    if isempty(x_true)
        PSNR_dB = nan; RMSE_val = nan; SSIM_val = nan; RD_percent = nan;
        target_mask = false(size(true_img)); bg_mask = false(size(true_img));
        metrics = struct('name', metricName, 'PSNR_dB', PSNR_dB, 'RMSE', RMSE_val, ...
                         'SSIM', SSIM_val, 'RD_percent', RD_percent, 'num_valid', 0, ...
                         'num_roi', nnz(roi_mask), 'num_target', 0, 'num_background', 0);
        return;
    end

    % RMSE / PSNR
    mse_val = mean((x_recon - x_true).^2);
    RMSE_val = sqrt(mse_val);
    peakVal = max(abs(x_true));
    if mse_val > 0 && peakVal > eps
        PSNR_dB = 10 * log10((peakVal^2) / mse_val);
    elseif mse_val <= eps
        PSNR_dB = inf;
    else
        PSNR_dB = nan;
    end

    % SSIM（全局统计形式）
    L_dyn = max(x_true) - min(x_true);
    if L_dyn <= eps
        if RMSE_val <= eps; SSIM_val = 1.0; else; SSIM_val = nan; end
    else
        C1_ssim = (0.01 * L_dyn)^2;
        C2_ssim = (0.03 * L_dyn)^2;
        mu_t = mean(x_true); mu_r = mean(x_recon);
        var_t = var(x_true, 1); var_r = var(x_recon, 1);
        cov_tr = mean((x_true - mu_t) .* (x_recon - mu_r));
        SSIM_val = ((2*mu_t*mu_r + C1_ssim) * (2*cov_tr + C2_ssim)) / ...
                   ((mu_t^2 + mu_r^2 + C1_ssim) * (var_t + var_r + C2_ssim));
    end

    % RD
    bg_ref = median(x_true, 'omitnan');
    dev_true = abs(true_img - bg_ref);
    dev_valid = dev_true(valid_mask);

    target_mask = false(size(true_img));
    bg_mask     = false(size(true_img));

    if isempty(dev_valid) || max(dev_valid) <= eps
        bg_mask = valid_mask;
        RD_percent = nan;
    else
        dev_target_th = prctile(dev_valid, targetDevPct);
        dev_bg_th     = prctile(dev_valid, bgDevPct);

        target_mask = valid_mask & (dev_true >= dev_target_th);
        bg_mask     = valid_mask & (dev_true <= dev_bg_th);

        if nnz(target_mask) < 10
            target_mask = valid_mask & (abs(true_img - bg_ref) > bg_ref_tol);
        end
        if nnz(bg_mask) < 10
            bg_mask = valid_mask & (abs(true_img - bg_ref) <= bg_ref_tol);
        end

        if nnz(target_mask) > 1 && nnz(bg_mask) > 1
            target_recon = recon_img(target_mask);
            target_true  = true_img(target_mask);
            x_bkgnd_scalar = mean(recon_img(bg_mask), 'omitnan');
            numerator_RD   = norm(target_recon - target_true);
            denominator_RD = norm(x_bkgnd_scalar - target_true);
            RD_percent = numerator_RD / max(denominator_RD, eps) * 100;
        else
            RD_percent = nan;
        end
    end

    metrics = struct('name', metricName, 'PSNR_dB', PSNR_dB, 'RMSE', RMSE_val, ...
                     'SSIM', SSIM_val, 'RD_percent', RD_percent, ...
                     'num_valid', nnz(valid_mask), 'num_roi', nnz(roi_mask), ...
                     'num_target', nnz(target_mask), 'num_background', nnz(bg_mask));
end

function w = keysCubicWeights(alpha, a)
% Keys cubic convolution kernel weights
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
