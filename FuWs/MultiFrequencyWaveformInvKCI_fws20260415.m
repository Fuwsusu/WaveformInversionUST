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

% Extract the Desired Frequency for Waveform Inversion (segmented steps)
% Low band:  0.30-0.45 MHz, step 0.025 MHz
% Mid band:  0.50-0.70 MHz, step 0.050 MHz
% High band: 0.775-1.25 MHz, step 0.075 MHz (force include 1.25)
fDATA_SoS = [ (0.30:0.025:0.45), (0.50:0.05:0.70), ...
    (0.775:0.075:1.225), 1.250 ]*(1e6); % Frequencies for SoS-only [Hz]
fDATA_SoSAtten = fDATA_SoS + 25e3;      % SoS/Attenuation [Hz]
fDATA = [fDATA_SoS, fDATA_SoSAtten];    % All Frequencies [Hz]

% Virtual element controls
enableVirtualElements = true;         % true: build/use virtual receive aperture
useVirtualTx = false;                 % true: virtual elements can also transmit
virtualElementsMode = 'unilateral';   % 'unilateral' (128->256) / 'bilateral' (128->384)
VE_freqMode = 'always';               % 'always' | 'threshold' | 'manual'
VE_freqThresh = 0.75e6;               % only used in threshold mode
VE_manualFlags = mod(0:numel(fDATA)-1,2)==1; % only used in manual mode
frac_shift = 0.20;                    % virtual element angular offset ratio

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

% Which Subset of Transmits to Use
dwnsmp = 1; % can be 1, 2, or 4 (faster but less quality with more downsampling)
orig_numElements = numElements;

% Build original geometry package
numElemLeftRightExcl = 63;
elemLeftRightExcl = -numElemLeftRightExcl:numElemLeftRightExcl;
elemInclude_orig = true(orig_numElements, orig_numElements);
for tx_element = 1:orig_numElements
    elemLeftRightExclCurrent = elemLeftRightExcl + tx_element;
    elemLeftRightExclCurrent(elemLeftRightExclCurrent<1) = orig_numElements + ...
         elemLeftRightExclCurrent(elemLeftRightExclCurrent<1);
    elemLeftRightExclCurrent(elemLeftRightExclCurrent>orig_numElements) = ...
        elemLeftRightExclCurrent(elemLeftRightExclCurrent>orig_numElements) - orig_numElements;
    elemInclude_orig(tx_element,elemLeftRightExclCurrent) = false;
end
tx_include_orig = 1:dwnsmp:orig_numElements;
REC_DATA_orig = REC_DATA;
REC_DATA_orig(isnan(REC_DATA_orig)) = 0;

% Remove outliers (original geometry)
perc_outliers = 0.99;
for f_idx = 1:numel(fDATA)
    REC_DATA_SINGLE_FREQ = REC_DATA_orig(tx_include_orig,:,f_idx);
    signalMagnitudes = elemInclude_orig(tx_include_orig,:).*abs(REC_DATA_SINGLE_FREQ);
    num_outliers = ceil((1-perc_outliers)*numel(signalMagnitudes));
    [~,idx_outliers] = maxk(signalMagnitudes(:),num_outliers);
    REC_DATA_SINGLE_FREQ(idx_outliers) = 0;
    REC_DATA_orig(tx_include_orig,:,f_idx) = REC_DATA_SINGLE_FREQ;
end

geo_orig.numElements = orig_numElements;
geo_orig.tx_include = tx_include_orig;
geo_orig.elemInclude = elemInclude_orig;
geo_orig.x_idx = x_idx;
geo_orig.y_idx = y_idx;
geo_orig.ind = ind;
geo_orig.x_circ = transducerPositionsX;
geo_orig.y_circ = transducerPositionsY;
geo_orig.REC_DATA = REC_DATA_orig(tx_include_orig,:,:);

% Build virtual element geometry package
geo_ve = [];
if enableVirtualElements
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
            theta_virtual = theta_base + sgn.*frac_shift.*dtheta_use;
            x_virtual = r_base.*cos(theta_virtual);
            y_virtual = r_base.*sin(theta_virtual);
            x_ve = reshape([x_base; x_virtual],1,[]);
            y_ve = reshape([y_base; y_virtual],1,[]);
            isRealRx_ve = false(1, numel(x_ve)); isRealRx_ve(1:2:end) = true; %#ok<NASGU>
        case 'bilateral'
            L_ve = 3;
            theta_left = theta_base - frac_shift.*dtheta_prev;
            theta_right = theta_base + frac_shift.*dtheta_next;
            x_left = r_base.*cos(theta_left); x_right = r_base.*cos(theta_right);
            y_left = r_base.*sin(theta_left); y_right = r_base.*sin(theta_right);
            x_ve = reshape([x_left; x_base; x_right],1,[]);
            y_ve = reshape([y_left; y_base; y_right],1,[]);
            isRealRx_ve = false(1, numel(x_ve)); isRealRx_ve(2:3:end) = true; %#ok<NASGU>
        otherwise
            error('Unknown virtualElementsMode: %s', virtualElementsMode);
    end

    N0 = orig_numElements;
    H_theta = zeros(L_ve*N0, N0, 'single');
    for n = 1:N0
        switch lower(virtualElementsMode)
            case 'unilateral'
                H_theta(2*n-1,n) = 1;
                H_theta(2*n,n) = 1; % zoh virtual channel
            case 'bilateral'
                H_theta(3*n-1,n) = 1;
                H_theta(3*n-2,n) = 1;
                H_theta(3*n,n) = 1;
        end
    end
    numElements_ve = size(H_theta,1);
    REC_DATA_ve = zeros(numElements_ve, numElements_ve, numel(fDATA), 'like', REC_DATA);
    for f_idx = 1:numel(fDATA)
        D0 = REC_DATA(:,:,f_idx);
        REC_DATA_ve(:,:,f_idx) = H_theta * D0 * H_theta.';
    end

    x_idx_ve = dsearchn(xi(:), x_ve(:));
    y_idx_ve = dsearchn(yi(:), y_ve(:));
    ind_ve = sub2ind([Nyi, Nxi], y_idx_ve, x_idx_ve);
    elemInclude_ve = true(numElements_ve, numElements_ve);
    numElemLeftRightExcl_ve = round(numElemLeftRightExcl*L_ve);
    elemLeftRightExcl_ve = -numElemLeftRightExcl_ve:numElemLeftRightExcl_ve;
    for tx_element = 1:numElements_ve
        ecur = elemLeftRightExcl_ve + tx_element;
        ecur(ecur<1) = numElements_ve + ecur(ecur<1);
        ecur(ecur>numElements_ve) = ecur(ecur>numElements_ve) - numElements_ve;
        elemInclude_ve(tx_element,ecur) = false;
    end

    if useVirtualTx
        tx_include_ve = 1:dwnsmp:numElements_ve;
    else
        switch lower(virtualElementsMode)
            case 'unilateral'
                tx_include_ve = 1:2:numElements_ve;
            case 'bilateral'
                tx_include_ve = 2:3:numElements_ve;
        end
        tx_include_ve = tx_include_ve(1:dwnsmp:end);
    end

    REC_DATA_ve(isnan(REC_DATA_ve)) = 0;
    for f_idx = 1:numel(fDATA)
        REC_DATA_SINGLE_FREQ = REC_DATA_ve(tx_include_ve,:,f_idx);
        signalMagnitudes = elemInclude_ve(tx_include_ve,:).*abs(REC_DATA_SINGLE_FREQ);
        num_outliers = ceil((1-perc_outliers)*numel(signalMagnitudes));
        [~,idx_outliers] = maxk(signalMagnitudes(:),num_outliers);
        REC_DATA_SINGLE_FREQ(idx_outliers) = 0;
        REC_DATA_ve(tx_include_ve,:,f_idx) = REC_DATA_SINGLE_FREQ;
    end

    geo_ve.numElements = numElements_ve;
    geo_ve.tx_include = tx_include_ve;
    geo_ve.elemInclude = elemInclude_ve;
    geo_ve.x_idx = x_idx_ve;
    geo_ve.y_idx = y_idx_ve;
    geo_ve.ind = ind_ve;
    geo_ve.x_circ = x_ve;
    geo_ve.y_circ = y_ve;
    geo_ve.REC_DATA = REC_DATA_ve(tx_include_ve,:,:);
end

% Frequency-wise VE schedule
VE_flags = false(1, numel(fDATA));
if enableVirtualElements
    switch lower(VE_freqMode)
        case 'always'
            VE_flags(:) = true;
        case 'threshold'
            VE_flags = (fDATA >= VE_freqThresh);
        case 'manual'
            if numel(VE_manualFlags) ~= numel(fDATA)
                error('VE_manualFlags length must equal numel(fDATA).');
            end
            VE_flags = logical(VE_manualFlags(:).');
        otherwise
            error('Unknown VE_freqMode: %s', VE_freqMode);
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
prevVEState = NaN;
for f_idx = 1:numel(fDATA)
    % Select original / virtual geometry for this frequency
    if VE_flags(f_idx) && ~isempty(geo_ve)
        geo_cur = geo_ve;
        veState = 1;
    else
        geo_cur = geo_orig;
        veState = 0;
    end
    tx_include_cur = geo_cur.tx_include;
    elemInclude_cur = geo_cur.elemInclude;
    x_idx_cur = geo_cur.x_idx;
    y_idx_cur = geo_cur.y_idx;
    ind_cur = geo_cur.ind;
    x_circ_cur = geo_cur.x_circ;
    y_circ_cur = geo_cur.y_circ;
    numElements_cur = geo_cur.numElements;
    REC_DATA_CUR = geo_cur.REC_DATA(:,:,f_idx);

    % If geometry switches between adjacent frequencies, reset CG memory
    if ~isnan(prevVEState) && (veState ~= prevVEState)
        search_dir = zeros(Nyi, Nxi);
        gradient_img_prev = zeros(Nyi, Nxi);
    end
    prevVEState = veState;

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
        gradient_img = zeros(Nyi,Nxi); 
        % Generate Sources
        SRC = zeros(Nyi, Nxi, numel(tx_include_cur));
        for elmt_idx = 1:numel(tx_include_cur)
            % Single Element Source
            x_idx_src = x_idx_cur(tx_include_cur(elmt_idx)); 
            y_idx_src = y_idx_cur(tx_include_cur(elmt_idx)); 
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
        % Build Adjoint Sources
        scaling = zeros(numel(tx_include_cur), 1);
        ADJ_SRC = zeros(Nyi, Nxi, numel(tx_include_cur));
        for elmt_idx = 1:numel(tx_include_cur)
            WVFIELD_elmt = WVFIELD(:,:,elmt_idx);
            REC_SIM = WVFIELD_elmt(ind_cur(elemInclude_cur(tx_include_cur(elmt_idx),:))); 
            REC = REC_DATA_CUR(elmt_idx, ...
                elemInclude_cur(tx_include_cur(elmt_idx),:)); REC = REC(:);
            scaling(elmt_idx) = (REC_SIM(:)'*REC(:)) / ...
                (REC_SIM(:)'*REC_SIM(:)); % Source Scaling
            ADJ_SRC_elmt = zeros(Nyi, Nxi);
            ADJ_SRC_elmt(ind_cur(elemInclude_cur(tx_include_cur(elmt_idx),:))) = ...
                scaling(elmt_idx)*REC_SIM - REC;
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
                hold on; plot(x_circ_cur(tx_include_cur(elmt_idx)), ...
                    y_circ_cur(tx_include_cur(elmt_idx)), 'r*')
                plot(x_circ_cur(elemInclude_cur(tx_include_cur(elmt_idx),:)), ...
                    y_circ_cur(elemInclude_cur(tx_include_cur(elmt_idx),:)), 'y.'); 
                drawnow; clf;
            end
        end
        % Remove Ringing from Gradient Image
        gradient_img = ringingRemovalFilt(xi, yi, gradient_img, c0, fDATA(f_idx), cutoff, ord);
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
        dREC_SIM = zeros(numel(tx_include_cur), numElements_cur);
        for elmt_idx = 1:numel(tx_include_cur)
            % Forward Projection of Search Direction Image
            PERTURBED_WVFIELD_elmt = PERTURBED_WVFIELD(:,:,elmt_idx);
            dREC_SIM(elmt_idx,elemInclude_cur(tx_include_cur(elmt_idx),:)) = -permute(scaling(elmt_idx) * ...
                PERTURBED_WVFIELD_elmt(ind_cur(elemInclude_cur(tx_include_cur(elmt_idx),:))),[2,1]);
        end
        % Step 4: Perform a Linear Approximation of Exact Line Search
        perc_step_size = 1; % (<1/2) Introduced to Improve Compliance with Strong Wolfe Conditions 
        alpha = -(gradient_img(:)'*search_dir(:))/(dREC_SIM(:)'*dREC_SIM(:));
        if updateAttenuation % Update Complex Slowness
            SI = sign(sign_conv) * imag(SLOW_ESTIM) + perc_step_size * alpha * search_dir;
            SLOW_ESTIM = real(SLOW_ESTIM) + 1i * sign(sign_conv) * SI;
        else
            SLOW_ESTIM = SLOW_ESTIM + perc_step_size * alpha * search_dir;
        end 
        VEL_ESTIM = 1./real(SLOW_ESTIM); % Wave Velocity Estimate [m/s]
        ATTEN_ESTIM = 2*pi*imag(SLOW_ESTIM)*sign(sign_conv);
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
