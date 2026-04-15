clear
clc
scriptTimer = tic; % Total script runtime (s)
% Add Functions to Path
addpath(genpath('Functions'));

% ------------------ Iteration/GIF/Display Controls ------------------
runAllIterations = false; % true=run all iterations, false=stop at requestedIterations
requestedIterations = 60; % Only effective when runAllIterations=false
saveGIF = true;           % Save per-iteration velocity animation GIF
gif_delay = 0.20;         % Delay time [s]
gif_save_every = 1;       % Save every N-th iteration
% Ground-truth options (Ali paper 512-element GIF)
useGifAsGroundTruth = true;
groundTruthGifPath = 'D:\Document_ING_fws\WaveformInversionUST\Results\start20260303\BenignCyst_WaveformInversionVEL_origin.gif';
groundTruthFrameIter = 60; % choose frame by iteration index in GIF
groundTruthCropRows = [];   % e.g., 60:500 if GIF has labels/colorbar margins
groundTruthCropCols = [];   % e.g., 40:560 if GIF has labels/colorbar margins
groundTruthFlipUD = false;  % set true when GIF vertical direction is opposite
% Soft red-blue colormap (blue=low, red=high)
n_cmap = 256; n2 = n_cmap/2;
cmap_rb = [
    linspace(0.17, 1.00, n2)', linspace(0.51, 1.00, n2)', linspace(0.83, 1.00, n2)';
    linspace(1.00, 0.84, n2)', linspace(1.00, 0.18, n2)', linspace(1.00, 0.15, n2)'
];
% --------------------------------------------------------------------
% ------------------ 稀疏阵元重建模块（参考 20260413 的 VE 思路） ------------------
% 模块目标：
%   1) 针对 sparse256 / 稀疏阵元采样，构造“虚拟阵元(VE)”以补足角向采样密度；
%   2) 在不引入正则化（TV/LF/Polar）与照明补偿（IllumComp）的前提下，
%      仅通过数据/几何侧改造提升梯度稳定性与可用视角覆盖；
%   3) 维持原 KCI 反演主流程（Helmholtz + NLCG + line search）不变。
%
% 关键开关说明：
%   enableVirtualElements : 是否启用 VE；false 时完全回退原始阵列流程。
%   virtualElementsMode   : VE 排布方式；当前脚本实现 unilateral（N->2N）。
%   frac_shift_lo/hi      : 低/高频 VE 偏移比例；低频保守、高频标准。
%   VE_freqMode           : 按频率启停 VE（always/threshold/manual）。
%   beta_ve               : 虚拟Rx残差权重，抑制插值通道偏差污染伴随源。
% -------------------------------------------------------------------------
enableVirtualElements = true;          % true=enable virtual elements for sparse-array compensation
useVirtualTx          = false;         % true=virtual elements also transmit, false=real Tx only
virtualElementsMode   = 'unilateral';  % 'unilateral' => N -> 2N
virtualInterpKernel   = 'zoh';         % 'zoh' or 'linear'
frac_shift_lo         = 0.10;          % conservative shift for lower frequencies
frac_shift_hi         = 0.20;          % default shift for higher frequencies
frac_shift_f_split    = 0.45e6;        % switch point between lo/hi shift [Hz]
VE_freqMode           = 'always';      % 'always' | 'threshold' | 'manual'
VE_freqThresh         = NaN;           % valid only when VE_freqMode='threshold'
beta_ve               = 0.50;          % virtual-Rx confidence weight (real=1, virtual=beta_ve)

% ------------------ 失配函数/步长/LF先验（从 V3_5 移植） ------------------
misfitType        = 'PolarPhase';  % 'L2' | 'PolarPhase'
alpha_hv          = 0.6;           % SoS阶段相位权重
alpha_hv_atten    = 0.1;           % Atten阶段相位权重
alpha_floor       = 1e-8;          % 步长下限，防止零步
auto_dv_ratio     = 3e-3;          % alpha_cap 自动步长上限比例

% 三阶段低频先验（真实数据默认开启）
enableLFPrior     = true;
lambda_stage      = 5e-4;
f_stage1_cutoff   = 0.475e6;
f_stage2_cutoff   = 0.850e6;
% -------------------------------------------------------------------------

% Load Extracted Dataset
filename = 'BenignCyst_sparse256'; % 'BenignCyst', 'Malignancy'
dataPath = ['SampleData/', filename, '.mat'];
load(dataPath, 'time', 'transducerPositionsXY', 'full_dataset');

% Optional ground truth for curve comparison / quantitative metrics
hasGroundTruth = false;
if exist(dataPath, 'file') == 2
    varsInFile = who('-file', dataPath);
    if all(ismember({'C', 'xi_orig', 'yi_orig'}, varsInFile))
        load(dataPath, 'C', 'xi_orig', 'yi_orig');
        hasGroundTruth = true;
    end
end

% Assemble Full Dataset
numElements = size(transducerPositionsXY,2);
fsa_rf_data = single(full_dataset);
clearvars full_dataset;

% Extract the Desired Frequency for Waveform Inversion
fDATA_SoS = (0.3:0.05:1.25)*(1e6); % Frequencies for SoS-only Iterations [Hz]
fDATA_SoSAtten = (0.325:0.05:1.275)*(1e6); % Frequencies for SoS/Attenuation Iterations [Hz]
fDATA = [fDATA_SoS, fDATA_SoSAtten]; % All Frequencies [Hz]

% Attenuation Iterations Always Happen After SoS Iterations
niterSoSPerFreq = [3*ones(size(fDATA_SoS)), 3*ones(size(fDATA_SoSAtten))]; % Number of Sound Speed Iterations Per Frequency
niterAttenPerFreq = [0*ones(size(fDATA_SoS)), 3*ones(size(fDATA_SoSAtten))]; % Number of Attenuation Iterations Per Frequency

% Discrete-Time Fourier Transform (DTFT) - Not an FFT though!
DTFT = exp(-1i*2*pi*fDATA'*time)*mean(diff(time));

% Geometric TOF Based on Sound Speed
c_geom = 1540; % Sound Speed [mm/us] 
transducerPositionsX = transducerPositionsXY(1,:);
transducerPositionsY = transducerPositionsXY(2,:);
geomTOFs = zeros(size(transducerPositionsXY,2));
for col = 1:size(transducerPositionsXY,2)
    geomTOFs(:,col) = ...
        sqrt((transducerPositionsX-transducerPositionsX(col)).^2 + ...
        (transducerPositionsY-transducerPositionsY(col)).^2)/c_geom;
end
clearvars col transducerPositionsXY;

% Window Time Traces - Extract the Frequencies for Waveform Inversion
REC_DATA = zeros(numElements, numElements, numel(fDATA), 'single');
sos_perc_change_pre = 0.05; sos_perc_change_post = Inf; 
twinpre = sos_perc_change_pre*max(geomTOFs(:));
twinpost = sos_perc_change_post*max(geomTOFs(:));
for tx_element = 1:numElements
    times_tx = geomTOFs(tx_element,:);
    [TIMES_TX, TIME] = meshgrid(times_tx, time);
    window = exp(-(1/2)*(subplus(TIME-TIMES_TX)/twinpost + ...
        subplus(TIMES_TX-TIME)/twinpre).^2);
    REC_DATA(tx_element,:,:) = ...
        permute(DTFT*(window.*fsa_rf_data(:,:,tx_element)),[2,1]);
end
clearvars twin window fsa_rf_data time times_tx TIME TIMES_TX DTFT;

%% Waveform Inversion

% Create Sound Speed Map and Transducer Ring
dxi = 0.3e-3; xmax = 120e-3;
xi = -xmax:dxi:xmax; yi = xi;
Nxi = numel(xi); Nyi = numel(yi);
[Xi, Yi] = meshgrid(xi, yi);
x_circ = transducerPositionsX;
y_circ = transducerPositionsY;
x_idx = dsearchn(xi(:), x_circ(:));
y_idx = dsearchn(yi(:), y_circ(:));
ind = sub2ind([Nyi, Nxi], y_idx, x_idx);
msk = zeros(Nyi, Nxi); msk(ind) = 1;

% Parameters
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

% Compute Phase Screen to Correct for Discretization of Element Positions
% Times for Discretized Positions Assuming Constant Sound Speed
x_circ_disc = xi(x_idx); y_circ_disc = yi(y_idx);
[x_circ_disc_tx, x_circ_disc_rx] = meshgrid(x_circ_disc, x_circ_disc);
[y_circ_disc_tx, y_circ_disc_rx] = meshgrid(y_circ_disc, y_circ_disc);
geomTOFs_disc = sqrt((x_circ_disc_tx-x_circ_disc_rx).^2 + ...
    (y_circ_disc_tx-y_circ_disc_rx).^2)/c_geom;
% Time-of-Flight Error for Discretized Positions
geomTOFs_error = geomTOFs_disc-geomTOFs;

% Phase Screen and Gain Correction
PS = zeros(numElements, numElements, numel(fDATA));
for f_idx = 1:numel(fDATA)
    PS(:,:,f_idx) = exp(1i*sign_conv*2*pi*fDATA(f_idx)*geomTOFs_error);
end
REC_DATA = REC_DATA .* PS; % Phase Screen and Gain Corrections
clearvars geom_distance col PS geomTOFs_error geomTOFs_disc geomTOFs;

% ================== 稀疏阵元重建：几何/数据预处理阶段 ==================
% 这里构造“几何套装（geometry package）”：
%   geo_orig : 原始稀疏阵列套装；
%   geo_ve_lo / geo_ve_hi : 两套 VE 套装（低频/高频 frac_shift）。
% 在主循环里按频率选择 geo_cur，实现“同一反演框架、不同采样几何”。
dwnsmp = 1; % can be 1, 2, or 4
REC_DATA(isnan(REC_DATA)) = 0; % Eliminate Blank Channel

% ---- 原始阵列近端接收排除（near-end exclusion）----
% 与 baseline 一致：每个 Tx 近邻 Rx 通道剔除，减小直达波/近场强耦合影响。
numElemLeftRightExcl = 63;
elemLeftRightExcl = -numElemLeftRightExcl:numElemLeftRightExcl;
elemInclude_orig = true(numElements, numElements);
for tx_element = 1:numElements 
    elemLeftRightExclCurrent = elemLeftRightExcl + tx_element;
    elemLeftRightExclCurrent(elemLeftRightExclCurrent<1) = numElements + ...
         elemLeftRightExclCurrent(elemLeftRightExclCurrent<1);
    elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements) = ...
        elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements) - numElements;
    elemInclude_orig(tx_element,elemLeftRightExclCurrent) = false;
end

% ---- 套装A：原始稀疏阵列 geo_orig ----
% 字段语义与 20260413 保持一致，便于主循环按 geo_cur 统一索引。
geo_orig.numElements          = numElements;
geo_orig.transducerPositionsX = transducerPositionsX;
geo_orig.transducerPositionsY = transducerPositionsY;
geo_orig.x_idx                = x_idx;
geo_orig.y_idx                = y_idx;
geo_orig.ind                  = ind;
geo_orig.elemInclude          = elemInclude_orig;
geo_orig.isRealRx             = true(1, numElements);
geo_orig.tx_include           = 1:dwnsmp:numElements;
geo_orig.REC_DATA             = REC_DATA(geo_orig.tx_include,:,:);

% ---- 套装A观测异常值抑制 ----
% 在每个频点按幅值分位剔除极端通道，降低坏道/异常回波对反演的干扰。
perc_outliers = 0.99; % Confidence Interval Cutoff
for f_idx = 1:numel(fDATA)
    REC_DATA_SINGLE_FREQ = geo_orig.REC_DATA(:,:,f_idx);
    signalMagnitudes = geo_orig.elemInclude(geo_orig.tx_include,:).*abs(REC_DATA_SINGLE_FREQ);
    num_outliers = ceil((1-perc_outliers)*numel(signalMagnitudes));
    [~,idx_outliers] = maxk(signalMagnitudes(:),num_outliers);
    REC_DATA_SINGLE_FREQ(idx_outliers) = 0;
    geo_orig.REC_DATA(:,:,f_idx) = REC_DATA_SINGLE_FREQ;
end

% ---- 套装B/C：虚拟阵元 geo_ve_lo / geo_ve_hi ----
% 处理流程（参考 20260413）：
%   (i)  基于环阵角度构造虚拟阵元坐标（unilateral 交替偏移）；
%   (ii) 用 H_theta 做角向重采样，得到扩展观测矩阵；
%   (iii) 用离散几何 TOF 差执行相位搬运，补偿网格离散误差；
%   (iv) 构造 isRealRx 掩码，后续用于 scaling 与伴随源加权。
if enableVirtualElements
    frac_shift_build_list = [frac_shift_lo, frac_shift_hi];
    geo_ve_built = cell(1,2);

    for fsi = 1:2
        fs_cur = frac_shift_build_list(fsi);
        x_base = transducerPositionsX(:).';
        y_base = transducerPositionsY(:).';
        theta_base = unwrap(atan2(y_base, x_base));
        dtheta_next = circshift(theta_base,-1) - theta_base;
        dtheta_next(dtheta_next<=0) = dtheta_next(dtheta_next<=0) + 2*pi;
        dtheta_prev = theta_base - circshift(theta_base,1);
        dtheta_prev(dtheta_prev<=0) = dtheta_prev(dtheta_prev<=0) + 2*pi;
        r_base = sqrt(x_base.^2 + y_base.^2);

        switch lower(virtualElementsMode)
            case 'unilateral'
                L_ve = 2;
                sgn = ones(size(theta_base)); sgn(2:2:end) = -1;
                dtheta_use = dtheta_next; dtheta_use(sgn<0) = dtheta_prev(sgn<0);
                theta_virtual = theta_base + sgn.*fs_cur.*dtheta_use;
                x_virtual = r_base.*cos(theta_virtual);
                y_virtual = r_base.*sin(theta_virtual);
                tx_ve = reshape([x_base; x_virtual],1,[]);
                ty_ve = reshape([y_base; y_virtual],1,[]);
            otherwise
                error('Unknown virtualElementsMode. Use ''unilateral''.');
        end

        N_ve = numElements;
        numElements_ve = L_ve * N_ve;
        base_map = repelem(1:N_ve, L_ve);
        x_idx_ve = dsearchn(xi(:), tx_ve(:));
        y_idx_ve = dsearchn(yi(:), ty_ve(:));
        ind_ve   = sub2ind([Nyi, Nxi], y_idx_ve, x_idx_ve);

        % H_theta: 角向重采样算子（虚拟阵元由真实阵元插值生成）
        H_theta = zeros(L_ve*N_ve, N_ve, 'single');
        for n = 1:N_ve
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
                otherwise
                    error('Unknown virtualInterpKernel. Use ''zoh'' or ''linear''.');
            end
        end

        % 角向重采样后生成扩展阵列观测（Tx/Rx 两侧同时扩展）
        REC_DATA_ve_full = zeros(L_ve*N_ve, L_ve*N_ve, numel(fDATA), 'like', REC_DATA);
        for f_idx = 1:numel(fDATA)
            D0 = REC_DATA(:,:,f_idx);
            REC_DATA_ve_full(:,:,f_idx) = H_theta * D0 * H_theta.';
        end

        % 相位搬运（copy phase transport）：
        %   在插值后的观测上，按 dTOF 施加 exp(i*2*pi*f*dTOF) 修正相位。
        x_base_disc = xi(x_idx); y_base_disc = yi(y_idx);
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

        % isRealRx_ve: 标记“真实接收通道”与“虚拟接收通道”
        % unilateral 排布下：奇数位为真实Rx，偶数位为虚拟Rx
        isRealRx_ve = false(1, numElements_ve);
        isRealRx_ve(1:L_ve:end) = true;
        if useVirtualTx
            tx_include_ve = 1:dwnsmp:numElements_ve;
        else
            tx_include_ve = 1:L_ve:numElements_ve;
            tx_include_ve = tx_include_ve(1:dwnsmp:end);
        end

        % VE 套装的 near-end exclusion（按扩展后阵元数同比例扩展）
        elemInclude_ve = true(numElements_ve, numElements_ve);
        numElemLeftRightExcl_ve = round(numElemLeftRightExcl * (numElements_ve / numElements));
        elemLeftRightExcl_ve = -numElemLeftRightExcl_ve:numElemLeftRightExcl_ve;
        for tx_element = 1:numElements_ve
            ec = elemLeftRightExcl_ve + tx_element;
            ec(ec<1) = numElements_ve + ec(ec<1);
            ec(ec>numElements_ve) = ec(ec>numElements_ve) - numElements_ve;
            elemInclude_ve(tx_element, ec) = false;
        end

        REC_DATA_ve_full(isnan(REC_DATA_ve_full)) = 0;
        for f_idx = 1:numel(fDATA)
            slice = REC_DATA_ve_full(tx_include_ve,:,f_idx);
            mag = elemInclude_ve(tx_include_ve,:) .* abs(slice);
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
        geo_ve_built{fsi} = g_tmp;
    end

    geo_ve_lo = geo_ve_built{1};
    geo_ve_hi = geo_ve_built{2};
end

% ---- VE 频率调度向量 VE_flags ----
% VE_flags(f_idx)=true  : 当前频率使用 VE 套装（lo/hi 再由分界频率决定）
% VE_flags(f_idx)=false : 使用原始稀疏阵列套装 geo_orig
VE_flags = false(1, numel(fDATA));
if enableVirtualElements
    switch lower(VE_freqMode)
        case 'always'
            VE_flags(:) = true;
        case 'threshold'
            VE_flags = (fDATA >= VE_freqThresh);
        case 'manual'
            VE_flags = mod(0:numel(fDATA)-1, 2) == 1;
        otherwise
            error('Unknown VE_freqMode');
    end
end

% Initial Constant Sound Speed Map [m/s]
c_init = 1480; % Initial Homogeneous Sound Speed [m/s] Guess
VEL_INIT = c_init*ones(Nyi,Nxi); 

% Initial Constant Attenuation [Np/(Hz m)]
ATTEN_INIT = 0*alphaNp*ones(Nyi,Nxi);

% (Nonlinear) Conjugate Gradient
search_dir = zeros(Nyi,Nxi); % Conjugate Gradient Direction
gradient_img_prev = zeros(Nyi,Nxi); % Previous Gradient Image
VEL_ESTIM = VEL_INIT; ATTEN_ESTIM = ATTEN_INIT;
SLOW_ESTIM = 1./VEL_ESTIM + ...
    1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi); % Initial Slowness Image [s/m]
VEL_stage1_ref = [];
VEL_stage2_ref = [];
stage1_ref_saved = false;
stage2_ref_saved = false;
crange = [1350, 1600]; % For reconstruction display [m/s]
attenrange = 10*[-1,1]; % For reconstruction display [dB/(cm MHz)]
c0 = mean(VEL_ESTIM(:)); cutoff = 0.75; ord = Inf; % Parameters for Ringing Removal Filter
if useGifAsGroundTruth
    if exist(groundTruthGifPath, 'file') == 2
        [C, xi_orig, yi_orig, hasGroundTruth] = readTruthFromGif( ...
            groundTruthGifPath, groundTruthFrameIter, xi, yi, crange, cmap_rb, ...
            groundTruthCropRows, groundTruthCropCols, groundTruthFlipUD);
    else
        warning('GT:GIFNotFound', 'GIF not found: %s', groundTruthGifPath);
    end
end
% Values to Save at Each Iteration
Niter = sum(niterSoSPerFreq)+sum(niterAttenPerFreq);
if requestedIterations >= Niter
    runAllIterations = true;
end
VEL_ESTIM_ITER = zeros(Nyi,Nxi,Niter);
ATTEN_ESTIM_ITER = zeros(Nyi, Nxi, Niter);
GRAD_IMG_ITER = zeros(Nyi,Nxi,Niter);
SEARCH_DIR_ITER = zeros(Nyi,Nxi,Niter);

% GIF output config
result_dir = 'D:\Document_ING_fws\WaveformInversionUST\Results\start20260303\';
gif_filepath = [result_dir, filename, '_WaveformInversionVEL_anim.gif'];
gif_initialized = false;
gif_iter_count = 0;
if saveGIF
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    fprintf('[GIF] will save to: %s\n', gif_filepath);
    hgif = figure(99); % GIF only tracks iterative wave-velocity figure
    set(hgif, 'Visible', 'on', 'Position', [100 100 640 560]);
    axgif = axes('Parent', hgif);
    hgif_im = imagesc(axgif, xi, yi, VEL_ESTIM, crange);
    axis(axgif, 'image');
    colorbar(axgif);
    colormap(axgif, cmap_rb);
    xlabel(axgif, 'Lateral [m]');
    ylabel(axgif, 'Axial [m]');
end
hmain = figure(1); % Main 2x2 iterative display (same as original workflow)
prev_VE_flag = VE_flags(1);
for f_idx = 1:numel(fDATA)
    % Iterations at Each Frequency
    for iter_f_idx = 1:(niterSoSPerFreq(f_idx)+niterAttenPerFreq(f_idx))
        iter = iter_f_idx + sum(niterSoSPerFreq(1:f_idx-1)) + ...
            sum(niterAttenPerFreq(1:f_idx-1));
        if ~runAllIterations && (iter > requestedIterations)
            fprintf('[STOP] reached requestedIterations=%d, stop inversion loop.\n', requestedIterations);
            break;
        end
        tic;
        % Step 1: Accumulate Backprojection Over Each Element
        % Reset CG at Each Frequency (SoS and Attenuation)
        if ((iter_f_idx == 1) || (iter_f_idx == 1+niterSoSPerFreq(f_idx)))
            search_dir = zeros(Nyi, Nxi); % Conjugate Gradient Direction
            gradient_img_prev = zeros(Nyi, Nxi); % Previous Gradient Image
        end
        % Attenuation Iterations Always Happen After SoS Iterations
        if iter_f_idx > niterSoSPerFreq(f_idx)
            updateAttenuation = true;
        else
            updateAttenuation = false; 
        end
        % ---- 按频率选择当前几何套装 geo_cur ----
        % 低频: geo_ve_lo（更保守偏移，抑制插值相位偏差）
        % 高频: geo_ve_hi（较强角向补偿）
        if VE_flags(f_idx) && enableVirtualElements
            if fDATA(f_idx) < frac_shift_f_split
                geo_cur = geo_ve_lo;
            else
                geo_cur = geo_ve_hi;
            end
        else
            geo_cur = geo_orig;
        end
        tx_include_cur  = geo_cur.tx_include;
        ind_cur         = geo_cur.ind;
        elemInclude_cur = geo_cur.elemInclude;
        isRealRx_cur    = geo_cur.isRealRx;
        REC_DATA_cur    = geo_cur.REC_DATA;

        % 若跨频率发生 VE 状态切换（orig <-> VE），重置 CG 历史量，
        % 避免不同几何下的梯度共轭关系被错误继承。
        if iter_f_idx == 1 && f_idx > 1 && (VE_flags(f_idx) ~= prev_VE_flag)
            search_dir = zeros(Nyi, Nxi);
            gradient_img_prev = zeros(Nyi, Nxi);
        end
        if iter_f_idx == 1
            prev_VE_flag = VE_flags(f_idx);
        end

        if enableLFPrior && ~stage1_ref_saved && fDATA(f_idx) > f_stage1_cutoff
            VEL_stage1_ref   = VEL_ESTIM;
            stage1_ref_saved = true;
            fprintf('[LF 3-stage] saved stage-1 ref at f=%.3f MHz\n', fDATA(f_idx)/1e6);
        end
        if enableLFPrior && ~stage2_ref_saved && fDATA(f_idx) > f_stage2_cutoff
            VEL_stage2_ref   = VEL_ESTIM;
            stage2_ref_saved = true;
            fprintf('[LF 3-stage] saved stage-2 ref at f=%.3f MHz\n', fDATA(f_idx)/1e6);
        end

        gradient_img = zeros(Nyi,Nxi); 
        % Generate Sources
        SRC = zeros(Nyi, Nxi, numel(tx_include_cur));
        for elmt_idx = 1:numel(tx_include_cur)
            % Single Element Source
            x_idx_src = geo_cur.x_idx(tx_include_cur(elmt_idx)); 
            y_idx_src = geo_cur.y_idx(tx_include_cur(elmt_idx)); 
            SRC(y_idx_src, x_idx_src, elmt_idx) = 1; 
        end
        % Forward Solve Helmholtz Equation
        HS = HelmholtzSolver(xi, yi, VEL_ESTIM, ATTEN_ESTIM, ...
            fDATA(f_idx), sign_conv, a0, L_PML);
        [WVFIELD, VIRT_SRC] = HS.solve(SRC,false);
        % Virtual Sources and Virtual (Sound Speed-Only) Wavefield for Attenuation
        if updateAttenuation % Virtual Sources for Attenuation
            % Modify Virtual Source for Attenuation
            VIRT_SRC = 1i*sign(sign_conv)*VIRT_SRC; 
        end
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
            rx_mask_all = elemInclude_cur(tx_id,:);
            rx_mask_scale = rx_mask_all & isRealRx_cur;
            REC_SIM_all = WVFIELD_elmt(ind_cur(rx_mask_all));
            REC_all = REC_DATA_cur(elmt_idx, rx_mask_all, f_idx); REC_all = REC_all(:);
            REC_SIM_scale = WVFIELD_elmt(ind_cur(rx_mask_scale));
            REC_scale = REC_DATA_cur(elmt_idx, rx_mask_scale, f_idx); REC_scale = REC_scale(:);
            % scaling 仅在真实Rx上估计，降低虚拟Rx插值偏差导致的幅相失真。
            scaling(elmt_idx) = (REC_SIM_scale(:)'*REC_scale(:)) / ...
                (REC_SIM_scale(:)'*REC_SIM_scale(:)); % Source Scaling with real-Rx only
            ADJ_SRC_elmt = zeros(Nyi, Nxi);
            p_sim_sc = scaling(elmt_idx) * REC_SIM_all;
            switch lower(misfitType)
                case 'l2'
                    residual = p_sim_sc - REC_all;
                case 'polarphase'
                    amp_sim  = abs(p_sim_sc) + eps;
                    phi_dir  = p_sim_sc ./ amp_sim;
                    amp_res  = amp_sim - abs(REC_all);
                    dphi     = angle(p_sim_sc .* conj(REC_all));
                    residual = (1 - alpha_hv_cur) .* amp_res .* phi_dir + ...
                               alpha_hv_cur       .* dphi    .* (1i * phi_dir);
                otherwise
                    error('Unknown misfitType: %s. Use ''L2'' or ''PolarPhase''.', misfitType);
            end
            phi_k_accum = phi_k_accum + sum(abs(residual).^2);
            % 伴随源通道加权：真实Rx=1，虚拟Rx=beta_ve
            % 目的：保留 VE 的角向补偿收益，同时抑制“非实测通道”过度主导梯度。
            w_channel = ones(sum(rx_mask_all), 1);
            w_channel(~isRealRx_cur(rx_mask_all)) = beta_ve;
            ADJ_SRC_elmt(ind_cur(rx_mask_all)) = w_channel .* residual;
            ADJ_SRC(:,:,elmt_idx) = ADJ_SRC_elmt;
        end
        % Backproject Error
        ADJ_WVFIELD = HS.solve(ADJ_SRC,true);
        SCALING = repmat(reshape(scaling, [1,1,numel(scaling)]), [Nyi, Nxi, 1]);
        BACKPROJ = -real(conj(SCALING.*VIRT_SRC).*ADJ_WVFIELD);
        % Accumulate Gradient Over Each Element
        for elmt_idx = 1:numel(tx_include_cur)
            % Accumulate Backprojection
            gradient_img = gradient_img + BACKPROJ(:,:,elmt_idx);
            showAnimation = false; % Show Backprojection Animation
            if showAnimation
                % Visualize Numerical Solution
                imagesc(xi,yi,gradient_img)
                xlabel('Lateral [m]'); ylabel('Axial [m]'); axis image;
                title(['Backprojection up to Transmit ', num2str(elmt_idx)]); 
                colorbar; colormap gray; drawnow;
                % Source and Receiver Positions Included
                hold on; plot(geo_cur.transducerPositionsX(tx_include_cur(elmt_idx)), ...
                    geo_cur.transducerPositionsY(tx_include_cur(elmt_idx)), 'r*')
                plot(geo_cur.transducerPositionsX(elemInclude_cur(tx_include_cur(elmt_idx),:)), ...
                    geo_cur.transducerPositionsY(elemInclude_cur(tx_include_cur(elmt_idx),:)), 'y.'); 
                drawnow; clf;
            end
        end
        % Remove Ringing from Gradient Image
        gradient_img = ringingRemovalFilt(xi, yi, gradient_img, c0, fDATA(f_idx), cutoff, ord);
        if enableLFPrior && ~updateAttenuation
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
        % Step 2: Compute New Conjugate Gradient Search Direction from Gradient
        % Conjugate Gradient Direction Scaling Factor for Updates
        if ((iter_f_idx == 1) || (iter_f_idx == 1+niterSoSPerFreq(f_idx)))
            beta = 0; 
        else 
            betaPR = (gradient_img(:)'*...
                (gradient_img(:)-gradient_img_prev(:))) / ...
                (gradient_img_prev(:)'*gradient_img_prev(:));
            betaFR = (gradient_img(:)'*gradient_img(:)) / ...
                (gradient_img_prev(:)'*gradient_img_prev(:));
            beta = min(max(betaPR,0),betaFR);
        end
        search_dir = beta*search_dir-gradient_img;
        gradient_img_prev = gradient_img;
        % Step 3: Compute Forward Projection of Current Search Direction
        PERTURBED_WVFIELD = HS.solve(VIRT_SRC.*search_dir, false); 
        dREC_SIM = zeros(numel(tx_include_cur), geo_cur.numElements);
        for elmt_idx = 1:numel(tx_include_cur)
            % Forward Projection of Search Direction Image
            PERTURBED_WVFIELD_elmt = PERTURBED_WVFIELD(:,:,elmt_idx);
            tx_id = tx_include_cur(elmt_idx);
            dREC_SIM(elmt_idx,elemInclude_cur(tx_id,:)) = -permute(scaling(elmt_idx) * ...
                PERTURBED_WVFIELD_elmt(ind_cur(elemInclude_cur(tx_id,:))),[2,1]);
        end
        % Step 4: Perform a Linear Approximation of Exact Line Search
        perc_step_size = 1; % (<1/2) Introduced to Improve Compliance with Strong Wolfe Conditions 
        den_ls   = real(dREC_SIM(:)'*dREC_SIM(:));
        num_ls   = -real(gradient_img(:)'*search_dir(:));
        alpha_ls = num_ls / (den_ls + eps);
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
        vmean = mean(VEL_ESTIM(:));
        target_dv_eff = auto_dv_ratio * vmean;
        sd_max      = max(abs(search_dir(:))) + eps;
        alpha_cap   = target_dv_eff / ((vmean^2) * sd_max + eps);
        alpha       = min(alpha_ls, alpha_cap);
        alpha       = max(alpha, alpha_floor);
        if updateAttenuation % Update Complex Slowness
            SI = sign(sign_conv) * imag(SLOW_ESTIM) + perc_step_size * alpha * search_dir;
            SLOW_ESTIM = real(SLOW_ESTIM) + 1i * sign(sign_conv) * SI;
        else
            SLOW_ESTIM = SLOW_ESTIM + perc_step_size * alpha * search_dir;
        end 
        VEL_ESTIM = 1./real(SLOW_ESTIM); % Wave Velocity Estimate [m/s]
        ATTEN_ESTIM = 2*pi*imag(SLOW_ESTIM)*sign(sign_conv);
        SLOW_ESTIM = 1./VEL_ESTIM + 1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi);
        % Save Intermediate Results
        VEL_ESTIM_ITER(:,:,iter) = VEL_ESTIM;
        ATTEN_ESTIM_ITER(:,:,iter) = ATTEN_ESTIM;
        GRAD_IMG_ITER(:,:,iter) = gradient_img;
        SEARCH_DIR_ITER(:,:,iter) = search_dir;
        % Visualize Numerical Solution
        figure(hmain);
        subplot(2,2,1); imagesc(xi,yi,VEL_ESTIM,crange);
        title(['Estimated Wave Velocity ', num2str(iter)]); axis image;
        xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap(cmap_rb);
        subplot(2,2,2); imagesc(xi,yi,Np2dB*slow2atten*ATTEN_ESTIM,attenrange);
        title(['Estimated Attenuation ', num2str(iter)]); axis image;
        xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap(cmap_rb);
        subplot(2,2,3); imagesc(xi,yi,search_dir)
        xlabel('Lateral [m]'); ylabel('Axial [m]'); axis image;
        title(['Search Direction Iteration ', num2str(iter)]); colorbar; colormap(cmap_rb); 
        subplot(2,2,4); imagesc(xi,yi,-gradient_img)
        xlabel('Lateral [m]'); ylabel('Axial [m]'); axis image;
        title(['Gradient Iteration ', num2str(iter)]); colorbar; colormap(cmap_rb); 
        drawnow; disp(['Iteration ', num2str(iter)]); toc;

        if saveGIF
            gif_iter_count = gif_iter_count + 1;
            if mod(gif_iter_count - 1, gif_save_every) == 0
                if ~exist('hgif', 'var') || ~isgraphics(hgif)
                    hgif = figure(99);
                    set(hgif, 'Visible', 'on', 'Position', [100 100 640 560]);
                end
                if ~exist('axgif', 'var') || ~isgraphics(axgif) || ...
                        ~exist('hgif_im', 'var') || ~isgraphics(hgif_im)
                    clf(hgif);
                    axgif = axes('Parent', hgif);
                    hgif_im = imagesc(axgif, xi, yi, VEL_ESTIM, crange);
                    axis(axgif, 'image');
                    colorbar(axgif);
                    colormap(axgif, cmap_rb);
                    xlabel(axgif, 'Lateral [m]');
                    ylabel(axgif, 'Axial [m]');
                else
                    set(hgif_im, 'CData', VEL_ESTIM);
                end
                set(get(axgif, 'Title'), 'String', ...
                    sprintf('Estimated Wave Velocity | iter=%d | f=%.3f MHz | %s', ...
                    iter, fDATA(f_idx)/1e6, ternary(updateAttenuation, 'SoS+Atten', 'SoS')), ...
                    'Interpreter', 'none');
                drawnow limitrate nocallbacks;
                frame = getframe(hgif);
                [A, map] = rgb2ind(frame2im(frame), 128);
                if ~gif_initialized
                    imwrite(A, map, gif_filepath, 'gif', 'LoopCount', Inf, 'DelayTime', gif_delay);
                    gif_initialized = true;
                else
                    imwrite(A, map, gif_filepath, 'gif', 'WriteMode', 'append', 'DelayTime', gif_delay);
                end
            end
        end
    end
    if ~runAllIterations && (iter >= requestedIterations)
        break;
    end
end

% Plot Final Reconstructions
figure(hmain);
subplot(1,2,1); imagesc(xi,yi,VEL_ESTIM,crange);
title(['Estimated Wave Velocity ', num2str(iter)]); axis image;
xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap(cmap_rb);
subplot(1,2,2); imagesc(xi,yi,Np2dB*slow2atten*ATTEN_ESTIM,attenrange); 
title(['Estimated Attenuation ', num2str(iter)]); axis image;
xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap(cmap_rb);

%% Curve comparison (blue=initial, red=final)
figure(7); clf;
set(gcf, 'Position', [80 120 1680 700]);

% Axial profile at x=0
cut_x_m = 0.025;
[~, ix_cut] = min(abs(xi - cut_x_m));
vel_axial_init = VEL_INIT(:, ix_cut);
vel_axial_final = VEL_ESTIM(:, ix_cut);

% Lateral profile at y=0
cut_y_m = 0.025;
[~, iy_cut] = min(abs(yi - cut_y_m));
vel_lateral_init = VEL_INIT(iy_cut, :);
vel_lateral_final = VEL_ESTIM(iy_cut, :);

% Schematic: where the two profiles are extracted
subplot(1,3,1);
imagesc(xi, yi, VEL_ESTIM, crange); axis image;
set(gca, 'YDir', 'reverse');
colormap(cmap_rb); colorbar;
hold on;
plot([cut_x_m, cut_x_m], [yi(1), yi(end)], 'b--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Axial cut x = %.1f mm', cut_x_m*1e3));
plot([xi(1), xi(end)], [cut_y_m, cut_y_m], 'r--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Lateral cut y = %.1f mm', cut_y_m*1e3));
xlabel('Lateral [m]'); ylabel('Axial [m]');
title('Profile extraction schematic');
legend('Location', 'southoutside'); grid on;

subplot(1,3,2);
plot(vel_axial_init, yi, 'b--', 'LineWidth', 1.6, 'DisplayName', 'Initial'); hold on;
plot(vel_axial_final, yi, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Final');
if hasGroundTruth
    [~, ix_gt] = min(abs(xi_orig - cut_x_m));
    vel_axial_true = C(:, ix_gt);
    plot(vel_axial_true, yi_orig, 'k-', 'LineWidth', 1.6, 'DisplayName', 'Ground truth');
end
set(gca, 'YDir', 'reverse');
ylabel('Axial [m]'); xlabel('Sound Speed [m/s]');
title(sprintf('Axial profile at x = %.1f mm', cut_x_m*1e3));
legend('Location', 'best'); grid on;
% Use data-driven x-limits so the profile spans the full axis width.
axial_all = [vel_axial_init(:); vel_axial_final(:)];
if hasGroundTruth
    axial_all = [axial_all; vel_axial_true(:)];
end
axial_x_min = min(axial_all);
axial_x_max = max(axial_all);
if axial_x_max > axial_x_min
    pad = 0.03*(axial_x_max - axial_x_min);
    xlim([axial_x_min - pad, axial_x_max + pad]);
else
    xlim([axial_x_min-1, axial_x_max+1]);
end
ylim([yi(1), yi(end)]);

subplot(1,3,3);
plot(xi, vel_lateral_init, 'b--', 'LineWidth', 1.6, 'DisplayName', 'Initial'); hold on;
plot(xi, vel_lateral_final, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Final');
if hasGroundTruth
    [~, iy_gt] = min(abs(yi_orig - cut_y_m));
    vel_lateral_true = C(iy_gt, :);
    plot(xi_orig, vel_lateral_true, 'k-', 'LineWidth', 1.6, 'DisplayName', 'Ground truth');
end
xlabel('Lateral [m]'); ylabel('Sound Speed [m/s]');
title(sprintf('Lateral profile at y = %.1f mm', cut_y_m*1e3));
legend('Location', 'best'); grid on;
xlim([xi(1), xi(end)]); ylim(crange);

%% Quantitative metrics
if hasGroundTruth
    [X_gt, Y_gt] = meshgrid(xi_orig, yi_orig);
    [X_est, Y_est] = meshgrid(xi, yi);
    C_true_on_est = interp2(X_gt, Y_gt, C, X_est, Y_est, 'linear', nan);
    valid_mask = isfinite(C_true_on_est) & isfinite(VEL_ESTIM);
    v_true = C_true_on_est(valid_mask);
    v_est = VEL_ESTIM(valid_mask);

    mse_vel = mean((v_est - v_true).^2);
    RMSE_val = sqrt(mse_vel);
    peakVel = max(v_true);
    if mse_vel > 0
        PSNR_dB = 10*log10((peakVel^2)/mse_vel);
    else
        PSNR_dB = inf;
    end

    L_dyn = max(v_true) - min(v_true);
    if L_dyn <= 0
        SSIM_val = 1.0;
    else
        C1_ssim = (0.01*L_dyn)^2;
        C2_ssim = (0.03*L_dyn)^2;
        mu_t = mean(v_true); mu_e = mean(v_est);
        var_t = var(v_true,1); var_e = var(v_est,1);
        cov_te = mean((v_true-mu_t).*(v_est-mu_e));
        SSIM_val = ((2*mu_t*mu_e + C1_ssim) * (2*cov_te + C2_ssim)) / ...
                   ((mu_t^2 + mu_e^2 + C1_ssim) * (var_t + var_e + C2_ssim));
    end

    bg_ref = median(C_true_on_est(valid_mask), 'omitnan');
    dev_true = abs(C_true_on_est - bg_ref);
    dev_valid = dev_true(valid_mask);
    target_mask = valid_mask & (dev_true >= prctile(dev_valid, 85));
    bg_mask = valid_mask & (dev_true <= prctile(dev_valid, 40));
    if nnz(target_mask) < 10
        target_mask = valid_mask & (abs(C_true_on_est - bg_ref) > 5.0);
    end
    if nnz(bg_mask) < 10
        bg_mask = valid_mask & (abs(C_true_on_est - bg_ref) <= 5.0);
    end
    if nnz(target_mask) > 1 && nnz(bg_mask) > 1
        target_est = VEL_ESTIM(target_mask);
        target_true = C_true_on_est(target_mask);
        bg_scalar = mean(VEL_ESTIM(bg_mask));
        RD_percent = norm(target_est - target_true) / ...
            max(norm(bg_scalar - target_true), eps) * 100;
    else
        RD_percent = nan;
    end

    fprintf('\n================ Quantitative Metrics (with GT) ================\n');
    fprintf('PSNR : %.4f dB\n', PSNR_dB);
    fprintf('RMSE : %.6f m/s\n', RMSE_val);
    fprintf('SSIM : %.6f\n', SSIM_val);
    fprintf('RD   : %.4f %%\n', RD_percent);
    fprintf('===============================================================\n\n');
else
    deltaV = VEL_ESTIM - VEL_INIT;
    RMSE_init = sqrt(mean(deltaV(:).^2));
    MAE_init = mean(abs(deltaV(:)));
    NRMSE_init_pct = RMSE_init / max(mean(VEL_INIT(:)), eps) * 100;
    [gx, gy] = gradient(VEL_ESTIM, h, h);
    TV_norm = mean(sqrt(gx.^2 + gy.^2), 'all');

    fprintf('\n=========== Quantitative Metrics (no GT, relative) ===========\n');
    fprintf('RMSE(final, init) : %.6f m/s\n', RMSE_init);
    fprintf('MAE(final, init)  : %.6f m/s\n', MAE_init);
    fprintf('NRMSE(final, init): %.4f %%\n', NRMSE_init_pct);
    fprintf('TV-norm(final)    : %.6f\n', TV_norm);
    fprintf('==============================================================\n\n');
end

if saveGIF && gif_initialized
    fprintf('[GIF] saved: %s\n', gif_filepath);
    try
        [gifFrame, gifMap] = imread(gif_filepath, 'gif', 'Frames', 1);
        figure(100); clf;
        imshow(gifFrame, gifMap);
        title('Saved GIF preview (first frame)');
        if exist('implay', 'file') == 2
            implay(gif_filepath);
        end
    catch ME
        if isempty(ME.identifier)
            warning('GIF:PreviewFailed', '%s', ME.message);
        else
            warning(ME.identifier, '%s', ME.message);
        end
    end
end
scriptElapsedSec = toc(scriptTimer);
fprintf('[TIME] total script: %.3f s\n', scriptElapsedSec);

% Save the Result to File
filename_results = ['D:\Document_ING_fws\WaveformInversionUST\Results\start20260303\', filename, '_WaveformInversionResults.mat'];
save(filename_results, '-v7.3', 'xi', 'yi', 'fDATA', 'niterAttenPerFreq', ...
    'niterSoSPerFreq', 'VEL_ESTIM_ITER', 'ATTEN_ESTIM_ITER', 'GRAD_IMG_ITER', ...
    'SEARCH_DIR_ITER')

function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end

function [C_true, xi_true, yi_true, isOk] = readTruthFromGif( ...
    gifPath, frameIter, xi_est, yi_est, crange, cmap_rb, cropRows, cropCols, flipUd)
isOk = false;
C_true = [];
xi_true = [];
yi_true = [];
try
    info = imfinfo(gifPath);
    nFrames = numel(info);
    frameIdx = min(max(round(frameIter), 1), nFrames);
    [A, map] = imread(gifPath, 'gif', 'Frames', frameIdx);
    if ndims(A) == 3
        rgb = im2double(A);
    else
        rgb = ind2rgb(A, map);
    end
    if ~isempty(cropRows)
        rgb = rgb(cropRows, :, :);
    end
    if ~isempty(cropCols)
        rgb = rgb(:, cropCols, :);
    end
    if flipUd
        rgb = flipud(rgb);
    end
    C_true = indexedToVelocity(rgb, crange, cmap_rb);
    [nRows, nCols] = size(C_true);
    xi_true = linspace(xi_est(1), xi_est(end), nCols);
    yi_true = linspace(yi_est(1), yi_est(end), nRows);
    isOk = true;
    fprintf('[GT] Use GIF as truth: %s | frame=%d/%d | size=%dx%d\n', ...
        gifPath, frameIdx, nFrames, nRows, nCols);
catch ME
    warning('GT:GIFLoadFailed', ...
        'Failed to load GIF truth (%s). Message: %s', gifPath, ME.message);
end
end

function velMap = indexedToVelocity(rgb, crange, cmap_rb)
rgb2 = reshape(rgb, [], 3);
cmap = double(cmap_rb);
dist2 = zeros(size(rgb2,1), size(cmap,1));
for k = 1:size(cmap,1)
    d = rgb2 - cmap(k,:);
    dist2(:,k) = sum(d.^2, 2);
end
[~, idx] = min(dist2, [], 2);
idx = reshape(idx, size(rgb,1), size(rgb,2));
velMap = crange(1) + (double(idx)-1)/(size(cmap,1)-1) * (crange(2)-crange(1));
end
