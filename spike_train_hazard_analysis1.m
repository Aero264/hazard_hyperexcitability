%% Spike-train analysis
% Set neuron_no and the analysis options below, then run this script.
% Uses base MATLAB and local functions. Run the file as a whole.

% Code written by Aravind Kumar Kamaraj. 

% INPUTS
% The recording file supplies trace and fs. 
% The function 'detect_spikes' (From Otomo et al 2020) must be on the path.
% 'detect_spikes' output fields:
%   pks_t                 Spike peak times in seconds, in recording order.
%   filtr                 Filtered voltage trace, sampled at Fs samples/s.
%   thrs_t, spkofs_t      Threshold and offset times for 'offset2thresh'.
%   thrs_bug, spkofs_bug  Optional flags for unusable thresholds or offsets.
%   det_params.Fs or .fs Sampling rate when opt.Fs is empty.
% These outputs are needed to run this code.

% OUTPUT OF THIS CODE: 
% Result is a scalar structure. 

% Hazard test decisions are 'reject', 'fail to reject', or 'not tested'.
% The monotonicity null is a nondecreasing hazard. 

% Class-resolved survival and hazard require two Vmin classes and rejection
% of the nondecreasing-hazard null. Classification is then evaluated by
% permuting class labels. Otherwise survival and hazard are estimated from
% all finite positive intervals as a single sample.
% The permutation null assumes exchangeable V_min class labels. 
% An unavailable p-value is NaN.

% Statistics:
% CV is the sample standard deviation divided by the mean. 
% Local variability uses the method selected below. 
% CohenKappa measures partition agreement
% SpearmanRho measures the rank association between V_min and ISI.

% REFERENCE:
% Data and the function 'detect_spikes' originally reported in:
% Otomo, K. et al. In vivo patch-clamp recordings reveal distinct 
% subthreshold signatures and threshold dynamics of midbrain dopamine neurons. 
% Nature Communications 11, 6286 (2020).
% Data repository URL: doi:10.7910/DVN/Q9ZXLZ

%% Analysis inputs and settings
% Load neuron recording

neuron_no = 64;
S = load("ccc/neuron" + neuron_no + ".mat", 'trace', 'fs');

% Analysis settings
opt = struct;
opt.Fs = S.fs;                         % Hz; needs data to be loaded first
opt.VminWindow = 'offset2thresh';     % 'offset2thresh' or 'peak2peak'.
opt.LocalVariabilityMethod = 'LV';   % 'LV' (Shinomoto) or 'CV2' (Holt).
opt.MinISI = 80;                     % Minimum number of usable ISIs.

% Nondecreasing-hazard test. The constrained bootstrap assumes iid ISIs.
opt.HMTF = 0.95;                     % Fraction used to set the test endpoint.
opt.TTTCutoff = [];                  % Fixed endpoint in seconds, if supplied.
opt.TTTBandwidth = [];               % [] uses cutoff*n^(-1/4).
opt.TTTPenalty = [];                 % [] uses 2*n^(-2/3).
opt.TTTBootstrap = 2000;              % 0 skips bootstrap inference.
opt.TTTAlpha = 0.02;                 % Reject when p is below this level.
opt.TTTSeed = [];                    % [] uses the current RNG stream.

% Shared settings for fitting V_min and log10(ISI).
opt.MaxClasses = 2;                 % Maximum classes in each fit (1 to 5).
opt.MinPeakProm = 0.10;              % Minimum prominence / KDE peak height.
opt.NumKdePoints = 200;
opt.MinModeN = 0.03;                 % >=1 is a count; (0,1) is a fraction
                                     % with a minimum of three; 0 disables
                                     % minimum size beyond nonempty classes.
opt.AsymmetryK = 1.5;                % Shoulder cut in units of fitted sigma.
opt.AsymmetryAlpha = 0.02;           % Two-sided tail symmetry threshold.

% Survival/hazard estimation and the V_min class-label permutation test.
opt.ComputeHazard = true;            % Also enables the permutation test.
opt.HazardGrid = [];                 % Optional vector of evaluation times (s).
opt.HazardGridPoints = 200;
opt.HazardGridQuantile = 0.995;       % Upper end of the automatic grid.
opt.HazardBandwidth = [];            % [] uses a scaled Silverman rule (s).
opt.BandwidthMode = 'fixed';          % 'fixed' or 'adaptive'.
opt.BandwidthK = 40;                  % Neighbours for adaptive bandwidth.
opt.HazardMinAtRisk = 10;            % Threshold for hazard tail truncation.
opt.HazardRelativeFloor = 1e-3;      % Relative hazard floor for comparisons.
opt.NumPermutations = 200;           % 0 skips the permutation test.
opt.MaxPermutationExceedances = 4;   % Maximum null statistics >= observed.
opt.PermSeed = 1;                   % [] uses the current RNG stream.

% Figure settings
opt.ShowFigures = true;              % Display interactive MATLAB figures.
opt.NeuronLabel = sprintf('Neuron %d', neuron_no);
opt.HistogramBins = 30;
opt.MaxTracePlotPoints = 200000;     % Decimation affects only the drawn trace.
lvMethod = opt.LocalVariabilityMethod;

validateattributes(opt.MaxClasses, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 5});

%% Main Analysis

% Spike detection from raw trace. Code from Otomo et al. 2020.
spk = detect_spikes(S.trace, S.fs, 'default', 'filtf', 1000);

if ~exist('spk', 'var') || ~isstruct(spk) || ~isscalar(spk)
    error('spike_train_hazard_analysis:BadInput', ...
        'spk must be a single (scalar) spk structure. To analyse several, run this script once per structure.');
end

result = struct('ISI', [], 'CV', NaN, 'LocalVariabilityMethod', lvMethod, ...
    'localVar', NaN, 'localVar_nPairs', 0, ...
    'Vmin', [], 'VminT',[], ...
    'nModes', 0, 'modeLabel', [], 'modeMeanVmin', [], 'modeCV', [], ...
    'modeLocalVar', [], 'modeLocalVar_nPairs', [], 'modeN', [], 'CohenKappa', NaN, 'SpearmanRho', NaN, 'excluded', false, ...
    'nUsableISI', 0, 'asymmetrySplit', false, 'asymmetryCut', NaN, ...
    'asymmetryP', NaN, ...
    'modeMeanISI', [], 'modeVarISI', [], ...
    'modeISI', {{}}, ...
    'ISImeanGlobal', NaN, 'ISIvarGlobal', NaN, 'ISIvarWithin', NaN, 'ISIvarBetween', NaN, ...
    'ISIvarTotal', NaN, 'hazardMonotonicityP', NaN, ...
    'hazardMonotonicityDecision', 'not tested', 'class_validity_p', NaN, ...
    'class_validity_decision', 'not tested');

Fs = opt.Fs;
if isempty(Fs)
    if isfield(spk, 'det_params') && isfield(spk.det_params, 'Fs')
        Fs = spk.det_params.Fs;
    elseif isfield(spk, 'det_params') && isfield(spk.det_params, 'fs')
        Fs = spk.det_params.fs;
    else
        error('spike_train_hazard_analysis:NoFs', ...
            'Set opt.Fs or supply the sampling rate in spk.det_params.Fs or .fs.');
    end
end

pks_t = spk.pks_t(:);
nSpk = numel(pks_t);

if nSpk < 2
    result.excluded = true;
    disp('Check the variable ''result'' for the analysis results');
    return
end

ISI = diff(pks_t);

% ISI statistics
[CV, localVar, localVar_nPairs] = local_isi_stats(ISI, lvMethod);

% Vmin computation
[Vmin, VminT] = local_compute_vmin(spk, Fs, opt.VminWindow);

%% Check for usuable number of intervals
nUsable = nnz(isfinite(ISI) & ISI > 0);
result.nUsableISI = nUsable;
if nUsable < opt.MinISI
    result.excluded = true;
    result.ISI = ISI;
    result.Vmin = Vmin;
    result.VminT = VminT;
    result.CV = CV;
    result.localVar = localVar;
    result.localVar_nPairs = localVar_nPairs;
    disp('Check the variable ''result'' for the analysis results');
    return
end
result.excluded = false;

%% Test for monotone hazard
[result.hazardMonotonicityP, hazardNonMonotone] = test_hazard_monotonicity( ...
    ISI, opt.TTTBootstrap, opt.HMTF, opt.TTTCutoff, opt.TTTBandwidth, ...
    opt.TTTPenalty, opt.TTTAlpha, opt.TTTSeed);
if isfinite(result.hazardMonotonicityP)
    result.hazardMonotonicityDecision = 'fail to reject - monotone';
    if hazardNonMonotone
        result.hazardMonotonicityDecision = 'reject - nonmonotone';
    end
end

%% Fitting Vmin and ISI classes
% Fit each distribution independently, keeping labels in interval order.
valid = isfinite(Vmin) & isfinite(ISI) & ISI > 0;
[modeLabel, nModes, result.asymmetrySplit, result.asymmetryCut, ...
    result.asymmetryP] = fit_distribution_classes(Vmin, valid, hazardNonMonotone, opt);

modeMeanVmin = nan(nModes, 1);
modeCV = nan(nModes, 1);
modeLocalVar = nan(nModes, 1);
modeLocalVar_nPairs = zeros(nModes, 1);
modeN = zeros(nModes, 1);
for m = 1:nModes
    sel = modeLabel == m;
    modeN(m) = sum(sel);
    modeMeanVmin(m) = mean(Vmin(sel));
    [modeCV(m), modeLocalVar(m), modeLocalVar_nPairs(m)] = ...
        local_isi_stats(ISI(sel), lvMethod);
end

% Fit ISI independently for partition agreement and plot boundaries.
isiValid = isfinite(ISI) & ISI > 0;
logISI = nan(size(ISI));
logISI(isiValid) = log10(ISI(isiValid));
[isiModeLabel, nISI] = fit_distribution_classes(logISI, isiValid, hazardNonMonotone, opt);
result.CohenKappa = compute_cohen_kappa(modeLabel, nModes, isiModeLabel, nISI);
result.SpearmanRho = compute_spearman_rho(Vmin(valid), ISI(valid));

%% ISI variance decomposition
% ISI variance decomposition across classes (law of total variance /
% ANOVA identity): ISIvarTotal = ISIvarWithin + ISIvarBetween, using
% population (N-denominator) variances throughout so the identity holds
% exactly. Operates on whichever classified ISIs ended up with a class
% label (same population used for modeCV/modeLocalVar above), so
% ISIvarGlobal (computed directly from the pooled classified ISIs) and
% ISIvarTotal (reconstructed from the per-class pieces) should agree to
% floating-point precision - a useful internal consistency check.
modeMeanISI = [];
modeVarISI = [];
modeISI = {};
ISImeanGlobal = NaN;
ISIvarGlobal = NaN;
ISIvarWithin = NaN;
ISIvarBetween = NaN;
ISIvarTotal = NaN;

if nModes >= 1
    modeMeanISI = nan(nModes, 1);
    modeVarISI = nan(nModes, 1);
    modeISI = cell(nModes, 1);
    for m = 1:nModes
        isiM = ISI(modeLabel == m);
        isiM = isiM(~isnan(isiM));
        modeISI{m} = isiM;
        if ~isempty(isiM)
            modeMeanISI(m) = mean(isiM);
            modeVarISI(m) = mean((isiM - modeMeanISI(m)).^2); % population variance
        end
    end

    validModes = modeN > 0 & ~isnan(modeMeanISI);
    if any(validModes)
        nTot = sum(modeN(validModes));
        ISImeanGlobal = sum(modeN(validModes) .* modeMeanISI(validModes)) / nTot;
        ISIvarWithin = sum(modeN(validModes) .* modeVarISI(validModes)) / nTot;
        ISIvarBetween = sum(modeN(validModes) .* (modeMeanISI(validModes) - ISImeanGlobal).^2) / nTot;
        ISIvarTotal = ISIvarWithin + ISIvarBetween;

        pooledISI = ISI(ismember(modeLabel, find(validModes)));
        pooledISI = pooledISI(~isnan(pooledISI));
        ISIvarGlobal = mean((pooledISI - ISImeanGlobal).^2); % population variance, direct
    end
end

%% Estimate pooled survival and hazard, with class resolution when eligible
hazard = [];
resolveHazardClasses = hazardNonMonotone && nModes == 2;
if opt.ComputeHazard
    [hazard, isiAll, labAll] = estimate_survival_and_hazard( ...
        ISI, modeLabel, nModes, resolveHazardClasses, opt);
    if resolveHazardClasses && ~isempty(hazard)
        [result.class_validity_p, rejectPermNull] = validate_class_fit( ...
            hazard, isiAll, labAll, nModes, opt);
        if isfinite(result.class_validity_p)
            result.class_validity_decision = 'fail to reject null - invalid';
            if rejectPermNull
                result.class_validity_decision = 'reject null - valid';
            end
        end
    end
end

%% Results and Plotting
result.ISI = ISI;
result.CV = CV;
result.LocalVariabilityMethod = lvMethod;
result.localVar = localVar;
result.localVar_nPairs = localVar_nPairs;
result.Vmin = Vmin;
result.VminT = VminT;
result.nModes = nModes;
result.modeLabel = modeLabel;
result.modeMeanVmin = modeMeanVmin;
result.modeCV = modeCV;
result.modeLocalVar = modeLocalVar;
result.modeLocalVar_nPairs = modeLocalVar_nPairs;
result.modeN = modeN;
result.modeMeanISI = modeMeanISI;
result.modeVarISI = modeVarISI;
result.modeISI = modeISI;
result.ISImeanGlobal = ISImeanGlobal;
result.ISIvarGlobal = ISIvarGlobal;
result.ISIvarWithin = ISIvarWithin;
result.ISIvarBetween = ISIvarBetween;
result.ISIvarTotal = ISIvarTotal;

disp('Check the variable ''result'' for the analysis results');

if opt.ShowFigures 
    % Show available estimates for inspection even if the fit was rejected.
    plot_spike_train_analysis(char(opt.NeuronLabel), ISI, Vmin, modeLabel, ...
        nModes, hazard, isiModeLabel, nISI, spk, Fs, VminT, ...
        result.hazardMonotonicityDecision, opt);
end


%% Local analysis functions
function [label, nModes, viaShoulder, cut, pValue] = ...
    fit_distribution_classes(x, valid, hazardNonMonotone, opt)
% Fit V_min or log10(ISI) using valleys between prominent Gaussian KDE
% peaks, with at most opt.MaxClasses peaks. Merge undersized classes.
% If only one eligible class remains, a rejected nondecreasing-hazard null
% enables the shoulder criterion. It checks for excess mass on one side
% of the KDE peak and places the boundary on that side.
% Labels increase with the class mean. Invalid observations retain NaN.
label = nan(size(x));
nModes = 0;
viaShoulder = false;
cut = NaN;
pValue = NaN;
nValid = nnz(valid);
if nValid < 4
    return
end

[validLabels, nModes] = local_partition_modes(x(valid), ...
    opt.MinPeakProm, opt.NumKdePoints, opt.MaxClasses);
label(valid) = validLabels;
minClassN = max(1, local_resolve_min_mode_n(opt.MinModeN, nValid));
if nModes > 1
    [label, nModes] = local_merge_small_modes(label, x, nModes, minClassN);
end

if nModes == 1 && opt.MaxClasses >= 2 && hazardNonMonotone
    x(~valid) = NaN;
    [shoulderLabel, nShoulder, shoulderCut, shoulderP] = ...
        local_asymmetry_split(x, opt.AsymmetryK, opt.AsymmetryAlpha, ...
        minClassN, opt.NumKdePoints);
    if nShoulder == 2 && all(accumarray(shoulderLabel(valid), 1, [2 1]) >= minClassN)
        label = shoulderLabel;
        nModes = 2;
        viaShoulder = true;
        cut = shoulderCut;
        pValue = shoulderP;
    end
end
end

function [CV, localVar, nPairs] = local_isi_stats(ISI, method)
%LOCAL_ISI_STATS  Coefficient of variation and local variability of an ISI vector.
%   Consecutive elements of ISI are treated as adjacent for the local
%   variability calculation (i.e. ISI should already be in the temporal
%   order you want pairs formed from).
%
%   method = 'CV2' (Holt et al. 1996):
%       CV2_i    = 2*|ISI(i)-ISI(i+1)| / (ISI(i)+ISI(i+1))
%       localVar = mean over all valid i.  Bounded in [0, 2).
%   method = 'LV' (Shinomoto et al. 2003):
%       LV_i     = ((ISI(i)-ISI(i+1)) / (ISI(i)+ISI(i+1)))^2
%       localVar = 3/nPairs * sum over all valid i.
%
%   nPairs is the number of adjacent ISI pairs localVar is based on
%   (n-1 for n valid ISIs) - use it to gauge how reliable a given value
%   is; a value based on 1-2 pairs is a legitimate computation but a
%   noisy estimate.
%
%   CV is sample standard deviation divided by mean duration. Both
%   measurements require at least two finite positive intervals.
ISI = ISI(:);
ISI = ISI(isfinite(ISI) & ISI > 0);
CV = NaN;
localVar = NaN;
nPairs = max(0, numel(ISI) - 1);
if nPairs > 0
    CV = std(ISI) / mean(ISI);
    a = ISI(1:end-1);
    b = ISI(2:end);
    switch lower(method)
        case 'cv2'
            localVar = mean(2 * abs(a - b) ./ (a + b));
        case 'lv'
            localVar = 3 / nPairs * sum(((a - b) ./ (a + b)).^2);
    end
end
end


% ------------------------------------------------------------------
function [Vmin, VminT] = local_compute_vmin(spk, Fs, windowMode)
%LOCAL_COMPUTE_VMIN  Minimum filtr voltage between consecutive spikes.
%   Vmin(k) is the trough between spike k and spike k+1, aligned with
%   ISI(k) = pks_t(k+1) - pks_t(k).

pks_t = spk.pks_t(:);
nSpk = numel(pks_t);
Vmin = nan(nSpk - 1, 1);
VminT = nan(nSpk - 1, 1);

filtr = spk.filtr(:);
nSamp = numel(filtr);

thrs_bug = local_get_bug(spk, 'thrs_bug', nSpk);
spkofs_bug = local_get_bug(spk, 'spkofs_bug', nSpk);

for k = 1:nSpk - 1
    switch lower(windowMode)
        case 'offset2thresh'
            if spkofs_bug(k) || thrs_bug(k + 1)
                continue
            end
            tStart = spk.spkofs_t(k);
            tEnd = spk.thrs_t(k + 1);
        case 'peak2peak'
            tStart = spk.pks_t(k);
            tEnd = spk.pks_t(k + 1);
    end

    if isempty(tStart) || isempty(tEnd) || isnan(tStart) || isnan(tEnd) || tEnd <= tStart
        continue
    end

    idxStart = max(1, round(tStart * Fs) + 1);
    idxEnd = min(nSamp, round(tEnd * Fs) + 1);

    if idxEnd < idxStart
        continue
    end

    [Vmin(k), rel] = min(filtr(idxStart:idxEnd));
    VminT(k) = (idxStart + rel - 2) / Fs;
end
end


% ------------------------------------------------------------------
function bug = local_get_bug(spk, fname, nSpk)
%LOCAL_GET_BUG  Fetch a *_bug logical vector, defaulting to all-false.
if isfield(spk, fname) && ~isempty(spk.(fname))
    bug = logical(spk.(fname)(:));
    if numel(bug) < nSpk
        bug(end + 1:nSpk) = false;
    end
else
    bug = false(nSpk, 1);
end
end


% ------------------------------------------------------------------
function [label, nModes] = local_partition_modes(x, minPeakProm, nPts, maxModes)
%LOCAL_PARTITION_MODES  KDE-based multimodality check and partitioning.
%   Builds a Gaussian KDE of x from scratch, counts prominent peaks
%   (KDE modes), and if more than one is found, splits x at the valleys
%   between consecutive peaks. Labels are re-ordered so label 1 has the
%   smallest mean(x), label 2 the next, etc.

x = x(:);
n = numel(x);

bw = 1.06 * std(x) * n^(-1 / 5); % Silverman's rule of thumb
if ~isfinite(bw) || bw <= 0
    bw = (max(x) - min(x)) / 20;
    if ~(bw > 0)
        bw = 1;
    end
end

pad = 3 * bw;
xi = linspace(min(x) - pad, max(x) + pad, nPts)';
f = zeros(nPts, 1);
for i = 1:nPts
    f(i) = sum(exp(-0.5 * ((xi(i) - x) / bw).^2));
end
f = f / (n * bw * sqrt(2 * pi));

pkIdx = local_find_peaks(f, minPeakProm);

if numel(pkIdx) > maxModes
    proms = local_peak_prominence(f, pkIdx);
    relProms = proms ./ max(reshape(f(pkIdx), size(proms)), eps);
    [~, ord] = sort(relProms, 'descend');
    pkIdx = sort(pkIdx(ord(1:maxModes)));
end

nModes = numel(pkIdx);

if nModes <= 1
    nModes = 1;
    label = ones(n, 1);
    return
end

pkIdx = sort(pkIdx);
valleyX = zeros(nModes - 1, 1);
for m = 1:nModes - 1
    seg = f(pkIdx(m):pkIdx(m + 1));
    [~, relIdx] = min(seg);
    valleyIdx = pkIdx(m) + relIdx - 1;
    valleyX(m) = xi(valleyIdx);
end

label = ones(n, 1);
for m = 1:numel(valleyX)
    label = label + (x > valleyX(m));
end

% Defensive re-sort so label 1 always has the smallest mean(x), even in
% unusual/edge-case KDE shapes.
meanPerLabel = nan(nModes, 1);
for m = 1:nModes
    meanPerLabel(m) = mean(x(label == m));
end
[~, order] = sort(meanPerLabel);
remap = zeros(nModes, 1);
remap(order) = 1:nModes;
label = remap(label);
end


% ------------------------------------------------------------------
function pkIdx = local_find_peaks(y, promFrac)
%LOCAL_FIND_PEAKS  Local maxima of y whose prominence is at least promFrac
%   times their own height.
%   Simple, dependency-free re-implementation of the "prominence"
%   concept behind Signal Processing Toolbox's findpeaks: a peak's
%   prominence is its height above the higher of the two nearest
%   valley floors on either side. The threshold is taken relative to each
%   peak's own height rather than to the global range of y, so that a low,
%   broad peak is not held to a bar set by the tallest one.
n = numel(y);
pkIdx = [];
for i = 2:n - 1
    if y(i) > y(i - 1) && y(i) >= y(i + 1)
        pkIdx(end + 1) = i; %#ok<AGROW>
    end
end
if isempty(pkIdx)
    return
end
proms = local_peak_prominence(y, pkIdx);
pkHeight = reshape(y(pkIdx), size(proms));
pkIdx = pkIdx(proms >= promFrac * pkHeight);
end


% ------------------------------------------------------------------
function proms = local_peak_prominence(y, pkIdx)
%LOCAL_PEAK_PROMINENCE  Approximate topographic prominence of each peak.
n = numel(y);
proms = zeros(size(pkIdx));
for k = 1:numel(pkIdx)
    idx = pkIdx(k);

    leftMin = y(idx);
    j = idx;
    while j > 1
        j = j - 1;
        if y(j) > y(idx)
            break
        end
        leftMin = min(leftMin, y(j));
    end

    rightMin = y(idx);
    j = idx;
    while j < n
        j = j + 1;
        if y(j) > y(idx)
            break
        end
        rightMin = min(rightMin, y(j));
    end

    proms(k) = y(idx) - max(leftMin, rightMin);
end
end

% ------------------------------------------------------------------
function [H, isiAll, labAll] = estimate_survival_and_hazard( ...
    ISI, modeLabel, nModes, resolveClasses, opt)
%ESTIMATE_SURVIVAL_AND_HAZARD  Pooled and per-class survival/hazard, plus three
%   reconstructions of the pooled hazard from the class-wise pieces.
%
%   Weighting identity: with w_m(t) = p_m*S_m(t)/sum_j p_j*S_j(t), which is
%   the fraction of the risk set at age t belonging to class m,
%       h(t) = sum_m w_m(t) * h_m(t)
%   holds exactly at the level of Nelson-Aalen increments, because
%       sum_m [n_m/n] * [d_m/n_m] = sum_m d_m / n = d/n.
%   Survival mixes linearly, S(t) = sum_m p_m*S_m(t), and hazards do not add.
%
%   References: Kaplan & Meier (1958) J. Am. Stat. Assoc. 53:457-481;
%   Nelson (1972) Technometrics 14:945-966; Aalen (1978) Ann. Statist.
%   6:701-726; Ramlau-Hansen (1983) Ann. Statist. 11:453-466 (kernel
%   smoothing of counting-process intensities); Silverman (1986) Density
%   Estimation for Statistics and Data Analysis (bandwidth rule); Schuster
%   (1985) Commun. Stat. Theory Methods 14:1123-1136 (reflection at the
%   support boundary).

H = [];
isiAll = [];
labAll = [];
grid = opt.HazardGrid;
bw = opt.HazardBandwidth;

keep = isfinite(ISI) & ISI > 0;
if resolveClasses
    % Use the same sample for the pooled and class-specific curves.
    keep = keep & ismember(modeLabel, 1:nModes);
    labAll = modeLabel(keep);
    labAll = labAll(:);
end
isiAll = ISI(keep);
isiAll = isiAll(:);
if numel(isiAll) < 4
    return
end
nAll = numel(isiAll);

% --- shared evaluation grid -------------------------------------------
if isempty(grid)
    grid = linspace(0, local_quantile_simple(isiAll, opt.HazardGridQuantile), opt.HazardGridPoints)';
end
grid = grid(:);

% --- shared bandwidth: one value for pooled AND every class, so that the
% reconstructions are not contaminated by differing amounts of smoothing.
if isempty(bw)
    q = local_quantile_simple(isiAll, [0.25 0.75]);
    sig = min(std(isiAll), (q(2) - q(1)) / 1.349);
    if ~isfinite(sig) || sig <= 0
        sig = std(isiAll);
    end
    bw = 1.5 * 0.9 * sig * nAll^(-1/5);
    if ~isfinite(bw) || bw <= 0
        % A constant sample has zero spread. Use its time scale as a
        % fallback so the kernel still has a positive width.
        bw = median(isiAll) / 20;
    end
end

% --- adaptive bandwidth (a vector over the grid) if requested ----------
if strcmpi(opt.BandwidthMode, 'adaptive')
    bw = local_adaptive_bandwidth(isiAll, grid, bw, opt.BandwidthK);
end

% --- pooled ------------------------------------------------------------
[S_pooled, h_pooled_raw, tLastPooled] = local_surv_hazard_grid(isiAll, grid, bw, opt.HazardMinAtRisk);
h_pooled = h_pooled_raw;
if ~isnan(tLastPooled)
    h_pooled(grid > tLastPooled) = NaN;      % display / reporting copy
end

if ~resolveClasses
    H = struct('t', grid, 'S_pooled', S_pooled, 'h_pooled', h_pooled);
    return
end

% --- per class ----------------------------------------------------------
nG = numel(grid);
S_mode = nan(nG, nModes);
h_mode = nan(nG, nModes);
h_mode_raw = zeros(nG, nModes);
modeTLast = nan(nModes, 1);
p = zeros(nModes, 1);

for m = 1:nModes
    isiM = isiAll(labAll == m);
    p(m) = numel(isiM) / nAll;
    if isempty(isiM)
        continue
    end
    [S_mode(:, m), hRawM, modeTLast(m)] = ...
        local_surv_hazard_grid(isiM, grid, bw, opt.HazardMinAtRisk);
    h_mode_raw(:, m) = hRawM;
    h_mode(:, m) = hRawM;
    if ~isnan(modeTLast(m))
        h_mode(grid > modeTLast(m), m) = NaN;   % display copy only
    end
end

% --- weights: w_m(t) = p_m S_m(t) / sum_j p_j S_j(t) -------------------
pS = S_mode .* repmat(p(:)', nG, 1);
pS(isnan(pS)) = 0;
denom = sum(pS, 2);
w = pS ./ repmat(max(denom, eps), 1, nModes);
w(denom <= 0, :) = 0;

S_recon = denom;

% --- three reconstructions --------------------------------------------
% Contributions from a class whose hazard has been truncated (risk set
% exhausted) are dropped only where that class carries negligible weight.
% Contributions are obtained by smoothing the WEIGHTED Nelson-Aalen
% increments, not by weighting the already-smoothed per-class hazards.
% Since w_m(t_i) = n_m(t_i)/n(t_i), the weighted increment is
%     w_m(t_i) * d_m,i/n_m(t_i) = d_m,i / n(t_i),
% i.e. class m's events over the POOLED risk set. Summing over m restores
% d_i/n(t_i) exactly, so sum_m contrib(:,m) equals h_pooled at every grid
% point to floating-point precision. Weighting after smoothing instead
% evaluates w at the grid point rather than at each event inside the kernel
% window, which makes the components overshoot the total wherever the
% hazard changes quickly - most visibly near t = 0.
contrib = local_mode_contributions(isiAll, labAll, nModes, grid, bw);
h_recon = sum(contrib, 2);

h_constw = sum(h_mode_raw .* repmat(p(:)', nG, 1), 2);
h_sum    = sum(h_mode_raw, 2);

% the reconstruction is only meaningful while the POOLED risk set survives
if ~isnan(tLastPooled)
    outside = grid > tLastPooled;
    contrib(outside, :) = NaN;
    h_recon(outside)  = NaN;
    h_constw(outside) = NaN;
    h_sum(outside)    = NaN;
end

% --- residuals over a common comparison mask ---------------------------
floorH = opt.HazardRelativeFloor * max([0; h_pooled(isfinite(h_pooled))]);
if isnan(tLastPooled)
    tLastPooled = max(grid);
end
mask = ~isnan(h_pooled) & h_pooled > floorH & grid <= tLastPooled;

H = struct( ...
    't', grid, 'bandwidth', bw, 'bandwidthMedian', median(bw(:)), 'p', p, ...
    'S_pooled', S_pooled, 'h_pooled', h_pooled, ...
    'S_mode', S_mode, 'h_mode', h_mode, 'w', w, 'h_contrib', contrib, ...
    'S_recon', S_recon, 'h_recon', h_recon, ...
    'h_constw', h_constw, 'h_sum', h_sum, ...
    'resid_S', max(abs(S_recon - S_pooled)), ...
    'resid_recon',  local_reldev(h_recon,  h_pooled, mask), ...
    'resid_constw', local_reldev(h_constw, h_pooled, mask), ...
    'resid_sum',    local_reldev(h_sum,    h_pooled, mask), ...
    'mask', mask, 'modeIFR', false(nModes, 1), 'modeTLast', modeTLast);

for m = 1:nModes
    hv = h_mode(:, m);
    ok = ~isnan(hv);
    if sum(ok) > 2
        H.modeIFR(m) = all(diff(hv(ok)) >= -1e-9);
    end
end

end

% ------------------------------------------------------------------
function [pValue, reject] = validate_class_fit(H, isiAll, labAll, nModes, opt)
% Test whether V_min class labels separate ISI durations. The statistic is
% the RMS relative deviation of the constant-weight hazard from the pooled
% hazard. Permuting labels preserves class sizes, the pooled estimate, the
% bandwidth and the comparison mask. The null assumes exchangeable labels.
pValue = NaN;
reject = false;
nPerm = round(opt.NumPermutations);
if nPerm > 0 && nModes >= 2 && any(H.mask)
    if ~isempty(opt.PermSeed)
        try
            rng(opt.PermSeed);
        catch
            rand('state', opt.PermSeed); %#ok<RAND>
        end
    end
    permVals = nan(nPerm, 1);
    for it = 1:nPerm
        labPerm = labAll(randperm(numel(labAll)));
        dperm = local_constw_resid(isiAll, labPerm, nModes, H.t, ...
            H.bandwidth, H.h_pooled, H.mask, opt.HazardMinAtRisk);
        permVals(it) = dperm(2);
    end
    obs = H.resid_constw(2);
    if isfinite(obs) && all(isfinite(permVals))
        nExceed = sum(permVals >= obs);
        pValue = (1 + nExceed) / (1 + nPerm);
        reject = nExceed <= opt.MaxPermutationExceedances;
    end
end
end


% ------------------------------------------------------------------
function d = local_constw_resid(isiAll, labAll, nModes, grid, bw, h_pooled, mask, minAtRisk)
%LOCAL_CONSTW_RESID  [max, rms] relative deviation of the constant-weight
%   reconstruction sum_m p_m*h_m(t) from h_pooled, for an arbitrary
%   labelling of the same intervals. Used for the permutation null.
nG = numel(grid);
N = numel(isiAll);
hm = nan(nG, nModes);
p = zeros(nModes, 1);
for m = 1:nModes
    xm = isiAll(labAll == m);
    p(m) = numel(xm) / N;
    if isempty(xm)
        continue
    end
    [~, hm(:, m)] = local_surv_hazard_grid(xm, grid, bw, minAtRisk);
end
hc = sum(hm .* repmat(p(:)', nG, 1), 2);
d = local_reldev(hc, h_pooled, mask);
end


% ------------------------------------------------------------------
function [S, h, tLast] = local_surv_hazard_grid(x, g, b, minAtRisk)
%LOCAL_SURV_HAZARD_GRID  Empirical survival and Ramlau-Hansen kernel hazard
%   of the sample x, both evaluated on the grid g with bandwidth b.
%   S(t) = P(ISI > t) is the empirical survival (identical to Kaplan-Meier
%   for a complete sample). The hazard smooths the Nelson-Aalen increments
%   d_i/n_i with an Epanechnikov kernel, reflected about t = 0. h is
%   returned untruncated; tLast reports the last duration with at least minAtRisk
%   intervals at risk so the caller can truncate where appropriate.

x = sort(x(:));
N = numel(x);
nG = numel(g);

S = zeros(nG, 1);
for i = 1:nG
    S(i) = sum(x > g(i)) / N;
end

[tE, ~, ic] = unique(x);
tE = tE(:); ic = ic(:);
dE = accumarray(ic, 1, [numel(tE) 1]);
nE = N - [0; cumsum(dE(1:end-1))];
dH = dE ./ nE;

if isscalar(b)
    bv = repmat(b, nG, 1);
else
    bv = b(:);
end

h = zeros(nG, 1);
for i = 1:nG
    u = (g(i) - tE) / bv(i);
    w = 0.75 * (1 - u.^2) .* (abs(u) <= 1);
    u2 = (g(i) + tE) / bv(i);                   % reflection about t = 0
    w = w + 0.75 * (1 - u2.^2) .* (abs(u2) <= 1);
    h(i) = sum(w .* dH) / bv(i);
end

% h is returned UNtruncated. tLast is the last duration with at least minAtRisk
% intervals at risk; callers truncate for display, but the decomposition
% uses the raw values so that a class whose risk set is exhausted
% contributes ~0 rather than NaN.
k = find(nE >= minAtRisk, 1, 'last');
if isempty(k)
    tLast = NaN;
else
    tLast = tE(k);
end
end


% ------------------------------------------------------------------
function d = local_reldev(v, ref, mask)
%LOCAL_RELDEV  [max, rms] relative deviation of v from ref over mask.
sel = mask & ~isnan(v) & ~isnan(ref) & ref ~= 0;
if ~any(sel)
    d = [NaN NaN];
    return
end
r = (v(sel) - ref(sel)) ./ ref(sel);
d = [max(abs(r)), sqrt(mean(r.^2))];
end


% ------------------------------------------------------------------
function q = local_quantile_simple(v, pr)
%LOCAL_QUANTILE_SIMPLE  Linear-interpolation quantile, no toolboxes.
v = sort(v(:));
m = numel(v);
if m == 0
    q = nan(size(pr));
    return
end
pos = (0.5 : m - 0.5)' / m;
q = interp1(pos, v, pr(:), 'linear');
q(pr(:) <= pos(1)) = v(1);
q(pr(:) >= pos(end)) = v(end);
q = reshape(q, size(pr));
end


% ------------------------------------------------------------------
function [label, nModes] = local_merge_small_modes(label, x, nModes, minN)
%LOCAL_MERGE_SMALL_MODES  Merge under-populated classes into their nearest
%   neighbour in the fitted variable. Repeatedly takes the smallest class; if it holds
%   fewer than minN intervals it is absorbed by the surviving class whose
%   mean is closest, and the remaining classes are renumbered. Stops
%   when every class has at least minN members, or one class is left.
%   Labels are re-sorted at the end so class 1 has the lowest mean of the fitted variable.

while nModes > 1
    cnt = zeros(nModes, 1);
    mu  = nan(nModes, 1);
    for m = 1:nModes
        sel = (label == m);
        cnt(m) = sum(sel);
        if cnt(m) > 0
            mu(m) = mean(x(sel));
        end
    end

    [minCount, m0] = min(cnt);
    if minCount >= minN
        break
    end

    if minCount > 0
        d = abs(mu - mu(m0));
        d(m0) = Inf;
        d(isnan(d)) = Inf;
        [~, mt] = min(d);
        label(label == m0) = mt;
    end

    % renumber the surviving labels to 1..nModes-1, preserving their order
    keepIdx = setdiff(1:nModes, m0);
    remap = zeros(nModes, 1);
    remap(keepIdx) = 1:numel(keepIdx);
    ok = ~isnan(label);
    label(ok) = remap(label(ok));
    nModes = nModes - 1;
end

% Preserve ordering by the mean of the fitted variable.
if nModes > 1
    mu = nan(nModes, 1);
    for m = 1:nModes
        sel = (label == m);
        if any(sel)
            mu(m) = mean(x(sel));
        end
    end
    [~, order] = sort(mu);
    remap = zeros(nModes, 1);
    remap(order) = 1:nModes;
    ok = ~isnan(label);
    label(ok) = remap(label(ok));
end
end


% ------------------------------------------------------------------
function contrib = local_mode_contributions(x, lab, nModes, g, b)
%LOCAL_MODE_CONTRIBUTIONS  Per-class contributions to the pooled kernel hazard.
%   contrib(:,m) = (1/b) * sum_i K((t - t_i)/b) * d_m,i / n(t_i), with d_m,i
%   the number of class-m intervals ending at t_i and n(t_i) the POOLED risk
%   set. This is the Ramlau-Hansen estimator applied to the survival-weighted
%   increments w_m(t_i)*dH_m,i, and summing over m recovers the pooled hazard
%   exactly, because sum_m d_m,i = d_i. Intervals with lab == 0 (no class)
%   are counted in the risk set but contribute to no component.

x = x(:);
lab = lab(:);
N = numel(x);
nG = numel(g);
g = g(:);

[xs, ord] = sort(x);
ls = lab(ord);

[tE, ~, ic] = unique(xs);
tE = tE(:); ic = ic(:);
dAll = accumarray(ic, 1, [numel(tE) 1]);
nE = N - [0; cumsum(dAll(1:end-1))];

dMode = zeros(numel(tE), nModes);
for m = 1:nModes
    dMode(:, m) = accumarray(ic, double(ls == m), [numel(tE) 1]);
end

if isscalar(b)
    bv = repmat(b, nG, 1);
else
    bv = b(:);
end

contrib = zeros(nG, nModes);
for i = 1:nG
    u = (g(i) - tE) / bv(i);
    w = 0.75 * (1 - u.^2) .* (abs(u) <= 1);
    u2 = (g(i) + tE) / bv(i);                    % reflection about t = 0
    w = w + 0.75 * (1 - u2.^2) .* (abs(u2) <= 1);
    contrib(i, :) = (w ./ nE)' * dMode / bv(i);
end
end


% ------------------------------------------------------------------
function bv = local_adaptive_bandwidth(x, g, b0, k)
%LOCAL_ADAPTIVE_BANDWIDTH  Nearest-neighbour kernel half-width over the grid.
%   bv(i) is the distance from g(i) to the k-th nearest observed interval,
%   rescaled so that median(bv) = b0 and clamped to [0.3, 5]*b0. Rescaling
%   makes this a redistribution of smoothing rather than a global change:
%   the window narrows where intervals are dense and widens where they are
%   sparse, while the overall level stays at the fixed bandwidth.
%
%   Nearest-neighbour bandwidths for hazard estimation: Mueller & Wang
%   (1994) Biometrics 50:61-76; the variable-bandwidth idea more generally
%   is Abramson (1982) Ann. Statist. 10:1217-1223.

xs = sort(x(:));
n = numel(xs);
g = g(:);
nG = numel(g);
k = min(max(round(k), 2), n);

bv = zeros(nG, 1);
for i = 1:nG
    d = sort(abs(xs - g(i)));
    bv(i) = d(k);
end

med = median(bv);
if isfinite(med) && med > 0
    bv = bv * (b0 / med);
else
    bv = repmat(b0, nG, 1);
end
bv = min(max(bv, 0.3 * b0), 5 * b0);
end


% ------------------------------------------------------------------
function rho = compute_spearman_rho(a, b)
%COMPUTE_SPEARMAN_RHO  Spearman rank correlation, ties averaged, no toolboxes.
a = a(:); b = b(:);
ok = ~isnan(a) & ~isnan(b);
a = a(ok); b = b(ok);
if numel(a) < 3
    rho = NaN;
    return
end
ra = local_tiedrank(a);
rb = local_tiedrank(b);
ra = ra - mean(ra);
rb = rb - mean(rb);
d = sqrt(sum(ra.^2) * sum(rb.^2));
if d <= 0
    rho = NaN;
else
    rho = sum(ra .* rb) / d;
end
end


% ------------------------------------------------------------------
function r = local_tiedrank(v)
%LOCAL_TIEDRANK  Ranks of v, tied values sharing their average rank.
v = v(:);
n = numel(v);
[sv, ord] = sort(v);
r = zeros(n, 1);
i = 1;
while i <= n
    j = i;
    while j < n && sv(j+1) == sv(i)
        j = j + 1;
    end
    r(ord(i:j)) = (i + j) / 2;
    i = j + 1;
end
end


% ------------------------------------------------------------------
function [pValue, reject] = test_hazard_monotonicity(ISI, nBoot, testFraction, cutoffOverride, bandwidth, penaltyOverride, alpha, seed)
%TEST_HAZARD_MONOTONICITY  Test whether an uncensored hazard is nondecreasing.
%
%   The primary L1 statistic tests
%
%       H0: the hazard h(t) is nondecreasing on [0, cutoff)
%
%   using the isotonic projection test of Groeneboom and Jongbloed.  It
%   replaced an earlier total-time-on-test (TTT) test. The TTT settings
%   retain their original names.
%
%   The empirical cumulative hazard Hn(t) = -log(1-Fn(t)) is projected onto
%   the cone of convex functions.  The primary statistic is
%
%       L1 = integral {Fn(t-) - Fiso(t)} dFn(t),
%
%   where Fiso = 1-exp(-GCM(Hn)).  Large values reject H0.  Critical values
%   are obtained by sampling from a smoothed, penalized, nondecreasing
%   isotonic-hazard estimate.  This is important: resampling the observed
%   intervals or TTT spacings does not, in general, reproduce the null law.
%
%   The function assumes independent, identically distributed, uncensored
%   positive intervals.  Serially dependent ISIs invalidate this iid
%   bootstrap.  For a point-process sequence, use a block/bootstrap model or
%   a conditional-intensity model instead.
%
%   References:
%     Groeneboom, P. and Jongbloed, G. (2012), Isotonic L2-projection test
%       for local monotonicity of a hazard, arXiv:1101.3306.
%     Groeneboom, P. and Jongbloed, G. (2013), Testing monotonicity of a
%       hazard: asymptotic distribution theory, Bernoulli 19, 1965-1999.

pValue = NaN;
reject = false;

x = ISI(:);
x = x(isfinite(x) & x > 0);
n = numel(x);
nBoot = round(nBoot);
if n < 12 || nBoot <= 0
    return
end
xs = sort(x);
if isempty(cutoffOverride)
    cutoff = local_empirical_cutoff(xs, testFraction);
else
    cutoff = cutoffOverride;
end

% A finite empirical cumulative hazard requires at least one observation
% above the endpoint.  Four observations below it is a bare computational
% minimum, not a claim that the asymptotics are reliable at that size.
nTest = sum(x < cutoff);
if ~isfinite(cutoff) || cutoff <= 0 || nTest < 4 || nTest >= n
    return
end

% Normalize the test interval to [0,1].  This makes the smoothing rule scale
% equivariant; multiplying every ISI by a constant leaves all p-values fixed.
z = x / cutoff;
a = 1;

if isempty(bandwidth)
    bw = n^(-1/4);
else
    bw = bandwidth / cutoff;
end
if bw > a
    error([mfilename ':BandwidthTooLarge'], ...
        'Bandwidth must not exceed the test-interval length.');
end
if isempty(penaltyOverride)
    penalty = 2*n^(-2/3);
else
    penalty = penaltyOverride;
end

observed = local_projection_statistic(z, a);
N = local_smooth_null(z, a, penalty, bw);

bootL1 = nan(nBoot, 1);

if ~isempty(seed)
    oldRng = rng;
    cleanupRng = onCleanup(@() rng(oldRng));
    rng(seed, 'twister');
end

for b = 1:nBoot
    zb = local_draw_from_null(N, n);
    bootL1(b) = local_projection_statistic(zb, a);
end
pValue = (1 + sum(bootL1 >= observed)) / (nBoot + 1);
reject = pValue < alpha;

end

% -------------------------------------------------------------------------
function cutoff = local_empirical_cutoff(xs, q)
% Choose a point strictly between observed levels, leaving a nonempty tail.

n = numel(xs);
[levels, ~, group] = unique(xs, 'sorted');
if numel(levels) < 2
    cutoff = NaN;
    return
end

counts = accumarray(group, 1);
cumCounts = cumsum(counts);
target = min(n - 1, max(4, floor(q*n)));
j = find(cumCounts <= target, 1, 'last');
if isempty(j)
    j = 1;
end
if j >= numel(levels)
    j = numel(levels) - 1;
end
cutoff = levels(j) + (levels(j+1) - levels(j))/2;

end

% -------------------------------------------------------------------------
function statistic = local_projection_statistic(x, a)
% Integrated empirical-versus-isotonic CDF gap on [0,a).

x = sort(x(:));
n = numel(x);
xin = x(x < a);
[xu, ~, group] = unique(xin, 'sorted');
counts = accumarray(group, 1);
leftCounts = [0; cumsum(counts(1:end-1))];

Fleft = leftCounts / n;
Hleft = -log1p(-Fleft);
FaLeft = numel(xin) / n;
HaLeft = -log1p(-FaLeft);

px = [0; xu; a];
py = [0; Hleft; HaLeft];
[hx, hy] = local_gcm(px, py);
Hiso = interp1(hx, hy, xu, 'linear');
Fiso = -expm1(-Hiso);
gap = max(0, Fleft - Fiso);
mass = counts / n;

statistic = sum(mass .* gap);

end

% -------------------------------------------------------------------------
function N = local_smooth_null(x, a, penalty, bw)
% Penalized isotonic hazard, smoothed with a triweight kernel.

x = sort(x(:));
n = numel(x);
xin = x(x < a);
[xu, ~, group] = unique(xin, 'sorted');
counts = accumarray(group, 1);
leftCounts = [0; cumsum(counts(1:end-1))];
Hleft = -log1p(-leftCounts/n);
HaLeft = -log1p(-numel(xin)/n);

% Equation (2.12) of Groeneboom and Jongbloed: lift the interior
% cumulative-sum points before taking the greatest convex minorant.
px = [0; xu; a];
py = [0; Hleft + penalty; HaLeft];
[hx, hy] = local_gcm(px, py);
slopes = diff(hy) ./ diff(hx);

% Numerical positivity is needed for an invertible cumulative hazard.
meanHazard = max(HaLeft/a, realmin('double'));
slopes = max(slopes, sqrt(eps)*meanHazard);

nGrid = min(4097, max(513, 8*n + 1));
t = linspace(0, a, nGrid)';
hStep = zeros(nGrid, 1);
for j = 1:numel(slopes)
    if j < numel(slopes)
        take = t >= hx(j) & t < hx(j+1);
    else
        take = t >= hx(j) & t <= hx(j+1);
    end
    hStep(take) = slopes(j);
end

dt = t(2) - t(1);
radius = max(1, ceil(bw/dt));
offset = (-radius:radius)' * dt;
kw = local_triweight(offset / bw);
if sum(kw) <= 0
    kw(radius+1) = 1;
end
kw = kw / sum(kw);

% Constant extension prevents boundary attenuation.  Convolution with a
% nonnegative kernel preserves monotonicity of the isotonic step function.
padded = [repmat(hStep(1), radius, 1); hStep; ...
          repmat(hStep(end), radius, 1)];
smoothed = conv(padded, kw, 'same');
h = smoothed(radius + (1:nGrid));
h = h(:);
for j = 2:nGrid
    if h(j) < h(j-1)
        h(j) = h(j-1);
    end
end
h = max(h, sqrt(eps)*meanHazard);

H = cumtrapz(t, h);
if H(end) > 0
    scale = HaLeft / H(end);
    h = h * scale;
    H = H * scale;
end

% Guard against duplicate cumulative-hazard grid values in interp1.
tinySlope = eps(max(1, H(end))) / dt;
h = max(h, tinySlope);
H = cumtrapz(t, h);

N.t = t;
N.h = h;
N.H = H;
N.a = a;

end

% -------------------------------------------------------------------------
function x = local_draw_from_null(N, n)
% Inverse-transform draw; continue past a with the terminal hazard.

e = -log(max(rand(n, 1), realmin('double')));
x = zeros(n, 1);
inside = e <= N.H(end);
if any(inside)
    [Hu, ia] = unique(N.H, 'stable');
    tu = N.t(ia);
    x(inside) = interp1(Hu, tu, e(inside), 'linear');
end
if any(~inside)
    x(~inside) = N.a + (e(~inside) - N.H(end)) / N.h(end);
end

end

% -------------------------------------------------------------------------
function k = local_triweight(u)
k = zeros(size(u));
inside = abs(u) <= 1;
k(inside) = (35/32) * (1 - u(inside).^2).^3;
end

% -------------------------------------------------------------------------
function [hx, hy] = local_gcm(px, py)
% Greatest convex minorant of ordered points, returned at its hull knots.

px = px(:);
py = py(:);
if numel(px) ~= numel(py) || any(diff(px) <= 0)
    error([mfilename ':BadCusumDiagram'], ...
        'Cumulative-sum diagram abscissae must be strictly increasing.');
end

n = numel(px);
hx = zeros(n, 1);
hy = zeros(n, 1);
k = 0;
for i = 1:n
    while k >= 2
        lhs = (hy(k) - hy(k-1)) * (px(i) - hx(k));
        rhs = (py(i) - hy(k)) * (hx(k) - hx(k-1));
        if lhs >= rhs
            k = k - 1;
        else
            break
        end
    end
    k = k + 1;
    hx(k) = px(i);
    hy(k) = py(i);
end
hx = hx(1:k);
hy = hy(1:k);

end



% ------------------------------------------------------------------
function minN = local_resolve_min_mode_n(v, nSample)
%LOCAL_RESOLVE_MIN_MODE_N  Interpret 'MinModeN' as a count or a proportion.
%   v >= 1 is an absolute interval count; 0 < v < 1 is a proportion of
%   nSample, floored at 3 so that a short recording cannot end up with a
%   threshold of 0 or 1 that admits any sliver the KDE happens to find.
if isempty(v) || v <= 0
    minN = 0;
elseif v >= 1
    minN = round(v);
else
    minN = max(3, round(v * nSample));
end
end

% ------------------------------------------------------------------
function [label, nModes, cut, pVal] = local_asymmetry_split(x, k, alpha, minTail, nKdePts)
%LOCAL_ASYMMETRY_SPLIT  Shoulder criterion based on one-sided excess mass.
%   Fits a Gaussian to the dominant KDE peak and asks whether the points
%   beyond +/- k*sigma are distributed symmetrically about it. A shoulder
%   can put most of them on one side. Skewness can also yield an
%   asymmetric tail count. The split goes at mu -/+ k*sigma on the heavy side.
%
%   Two design points. The centre is the KDE peak and the scale is the
%   NARROWER half-width at half-maximum, i.e. the side a shoulder has not
%   inflated; a least-squares Gaussian fit to all the data is dragged toward
%   the very shoulder being looked for. And the decision is a symmetry test
%   BETWEEN the tails rather than a comparison against the Gaussian's own
%   tail mass, because an error in sigma inflates both tails equally and
%   cancels, whereas comparing to the Gaussian inherits the sigma error.
%
%   This test is deliberately liberal and is meant to run only when a
%   separate criterion has already indicated structure.

label = nan(size(x));
nModes = 1;
cut = NaN;
pVal = NaN;

v = x(~isnan(x));
n = numel(v);
if n < 12
    label(~isnan(x)) = 1;
    return
end

% --- KDE, same rule as local_partition_modes ---------------------------
bw = 1.06 * std(v) * n^(-1/5);
if ~isfinite(bw) || bw <= 0
    bw = (max(v) - min(v)) / 20;
end
if ~(bw > 0)
    label(~isnan(x)) = 1;
    return
end
xi = linspace(min(v) - 3*bw, max(v) + 3*bw, nKdePts)';
f = zeros(nKdePts, 1);
for i = 1:nKdePts
    f(i) = sum(exp(-0.5 * ((xi(i) - v) / bw).^2));
end
f = f / (n * bw * sqrt(2*pi));

[~, iPk] = max(f);
mu = xi(iPk);

% --- scale from the narrower half-width at half-maximum ----------------
half = f(iPk) / 2;
iL = iPk;
while iL > 1 && f(iL) > half
    iL = iL - 1;
end
iR = iPk;
while iR < nKdePts && f(iR) > half
    iR = iR + 1;
end
wL = mu - xi(iL);
wR = xi(iR) - mu;
w = min(wL, wR);
if ~(w > 0)
    w = max(wL, wR);
end
sigma = w / sqrt(2 * log(2));
if ~isfinite(sigma) || sigma <= 0
    label(~isnan(x)) = 1;
    return
end

% --- tail symmetry -----------------------------------------------------
nLo = sum(v < mu - k*sigma);
nHi = sum(v > mu + k*sigma);
tot = nLo + nHi;
if tot < minTail || tot > n/4
    label(~isnan(x)) = 1;
    return
end

hi = max(nLo, nHi);
pVal = 2 * local_binom_upper(hi, tot, 0.5);
pVal = min(pVal, 1);

if pVal >= alpha || hi < minTail
    label(~isnan(x)) = 1;
    return
end

if nLo > nHi
    cut = mu - k*sigma;
else
    cut = mu + k*sigma;
end

label(~isnan(x)) = 1 + double(x(~isnan(x)) > cut);
nModes = 2;

% Order the two classes by the mean of the fitted variable.
m1 = mean(x(label == 1));
m2 = mean(x(label == 2));
if m1 > m2
    sw = label;
    label(sw == 1) = 2;
    label(sw == 2) = 1;
end
end


% ------------------------------------------------------------------
function p = local_binom_upper(kObs, nTot, prob)
%LOCAL_BINOM_UPPER  P(X >= kObs) for X ~ Binomial(nTot, prob), no toolboxes.
if kObs <= 0
    p = 1;
    return
end
j = (kObs:nTot)';
logC = gammaln(nTot + 1) - gammaln(j + 1) - gammaln(nTot - j + 1);
p = sum(exp(logC + j*log(prob) + (nTot - j)*log(1 - prob)));
p = min(max(p, 0), 1);
end


% ------------------------------------------------------------------
function kappa = compute_cohen_kappa(modeLabel, nModes, isiModeLabel, nISI)
% Match class labels to maximize agreement before computing Cohen's kappa.
% Report NaN if class counts differ or expected agreement is one.
kappa = NaN;
both = ~isnan(modeLabel) & ~isnan(isiModeLabel);
n = nnz(both);
if n < 2 || nModes ~= nISI || nModes < 2
    return
end

conf = accumarray([modeLabel(both), isiModeLabel(both)], 1, [nModes, nISI]);
rowMarg = sum(conf, 2);
colMarg = sum(conf, 1);
P = perms(1:nISI);
best = -1;
matching = [];
for k = 1:size(P, 1)
    hit = sum(conf(sub2ind(size(conf), P(k, :), 1:nISI)));
    if hit > best
        best = hit;
        matching = P(k, :);
    end
end

expectedAgreement = (colMarg * rowMarg(matching(:))) / n^2;
if expectedAgreement < 1
    kappa = (best / n - expectedAgreement) / (1 - expectedAgreement);
end
end

%% Figures
function plot_spike_train_analysis(fig_label, ISI, Vmin, modeLabel, ...
    nModes, hazard, isiModeLabel, nISI, spk, Fs, VminT, monotonicityDecision, opt)
% Display the original per-neuron figures and the four-panel decomposition.
% Class colours follow the V_min labels across figures. Grey data markers
% and histogram bars denote intervals without a V_min class assignment.

if nModes <= 1
    cols = [0 0.6 0];
    nModes = 1;
else
    palette = [0.5 0 1; 0 0.5 1; 0.85 0.33 0.10; 0 0.6 0; 0.75 0.20 0.45];
    cols = palette(1:nModes,:);
end
nb = opt.HistogramBins;
visibility = 'off';
if opt.ShowFigures
    visibility = 'on';
end
figures = gobjects(1, 7);
fileNames = {'Vmin', 'ISI', 'scatter', 'trace', 'hazard', 'survival', 'hazard_decomp'};

% Each dashed line sits between adjacent classes in the fitted variable.
logISI = nan(size(ISI));
pos = isfinite(ISI) & ISI > 0;
logISI(pos) = log10(ISI(pos));
vminBoundaries = zeros(1, nModes - 1);
for m = 1:nModes - 1
    vminBoundaries(m) = (max(Vmin(modeLabel == m)) + min(Vmin(modeLabel == m + 1))) / 2;
end
isiBoundaries = [];
if nISI > 1
    isiBoundaries = zeros(1, nISI - 1);
    for m = 1:nISI - 1
        isiBoundaries(m) = (max(logISI(isiModeLabel == m)) + min(logISI(isiModeLabel == m + 1))) / 2;
    end
end

% --- 1. V_min histogram, stacked by class ------------------------------
v = Vmin(isfinite(Vmin));
if ~isempty(v)
    vPad = (min(v) == max(v)) * 0.01 * max(1, abs(v(1)));
    edges = linspace(min(v) - vPad, max(v) + vPad, nb + 1)';
    ctr = edges(1:end-1) + diff(edges)/2;
    C = zeros(nb, nModes);
    for m = 1:nModes
        C(:,m) = histcounts(Vmin(modeLabel == m), edges).';
    end
    unclassified = isfinite(Vmin) & isnan(modeLabel);
    barColors = cols;
    if any(unclassified)
        C(:,end+1) = histcounts(Vmin(unclassified), edges).';
        barColors(end+1,:) = [0.7 0.7 0.7];
    end

    figures(1) = figure('Color', 'w', 'Visible', visibility, 'Name', [fig_label ' | Vmin histogram']);
    hold on; box off
    hB = bar(ctr, C, 'stacked');
    for m = 1:numel(hB)
        set(hB(m), 'FaceColor', barColors(m,:), 'EdgeColor', [0.2 0.2 0.2]);
    end
    vminTitle = {sprintf('%s - V_{min} histogram', fig_label); ...
        'Bar colors: V_{min} classes'};
    if ~isempty(vminBoundaries)
        xline(vminBoundaries, 'k--', 'LineWidth', 1.8);
        vminTitle{2} = [vminTitle{2} '; dashed: V_{min} class boundaries'];
    end
    xlabel('V_{min} (mV)'); ylabel('Counts');
    title(vminTitle, 'Interpreter', 'tex');
end

% --- 2. log10(ISI) histogram, stacked by V_min class -------------------
d = logISI(pos);
if ~isempty(d)
    isiPad = (min(d) == max(d)) * 0.01 * max(1, abs(d(1)));
    edges = linspace(min(d) - isiPad, max(d) + isiPad, nb + 1)';
    ctr = edges(1:end-1) + diff(edges)/2;
    C = zeros(nb, nModes);
    for m = 1:nModes
        C(:,m) = histcounts(logISI(modeLabel == m), edges).';
    end
    unclassified = isfinite(logISI) & isnan(modeLabel);
    barColors = cols;
    if any(unclassified)
        C(:,end+1) = histcounts(logISI(unclassified), edges).';
        barColors(end+1,:) = [0.7 0.7 0.7];
    end

    figures(2) = figure('Color', 'w', 'Visible', visibility, 'Name', [fig_label ' | ISI histogram']);
    hold on; box off
    hB = bar(ctr, C, 'stacked');
    for m = 1:numel(hB)
        set(hB(m), 'FaceColor', barColors(m,:), 'EdgeColor', [0.2 0.2 0.2]);
    end
    xlabel('log_{10} ISI (s)'); ylabel('Counts');

    isiTitle = {sprintf('%s - log_{10}(ISI) histogram', fig_label); ...
        'Bar colors: V_{min} classes'};
    if ~isempty(isiBoundaries)
        xline(isiBoundaries, 'k--', 'LineWidth', 1.8);
        isiTitle{2} = [isiTitle{2} '; dashed: log_{10}(ISI) class boundaries'];
    end
    title(isiTitle, 'Interpreter', 'tex');
end

% --- 3. V_min vs log10(ISI) scatter, coloured by V_min class -----------
figures(3) = figure('Color', 'w', 'Visible', visibility, 'Name', [fig_label ' | Vmin vs ISI']);
hold on; box off
for m = 1:nModes
    sel = (modeLabel == m);
    if any(sel)
        plot(Vmin(sel), logISI(sel), 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', cols(m,:), 'MarkerEdgeColor', 'none');
    end
end
un = isnan(modeLabel) & ~isnan(Vmin) & ~isnan(logISI);
if any(un)
    plot(Vmin(un), logISI(un), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
end

% discrimination boundaries: vertical for the V_min partition, horizontal
% for the ISI partition. Taken midway between the adjacent classes' extreme
% values, which is where the labels actually change - after any MinModeN
% merging, that is not necessarily the original KDE valley.
xl = [min(Vmin(~isnan(Vmin))), max(Vmin(~isnan(Vmin)))];
yl = [min(logISI(~isnan(logISI))), max(logISI(~isnan(logISI)))];
if numel(xl) == 2 && diff(xl) > 0 && numel(yl) == 2 && diff(yl) > 0
    xlim(xl); ylim(yl);
end

scatterTitle = {sprintf('%s - V_{min} vs log_{10}(ISI)', fig_label); ...
    'Point colors: V_{min} classes'};
if ~isempty(vminBoundaries)
    xline(vminBoundaries, 'k--', 'LineWidth', 1.0);
    scatterTitle{end+1} = 'Vertical dashed: V_{min} class boundaries';
end
if ~isempty(isiBoundaries)
    yline(isiBoundaries, 'k--', 'LineWidth', 1.0);
    scatterTitle{end+1} = 'Horizontal dashed: log_{10}(ISI) class boundaries';
end
xlabel('V_{min} (mV)'); ylabel('log_{10} ISI (s)');
title(scatterTitle, 'Interpreter', 'tex');

% --- 4. filtered trace with V_min marked by class -----------------------
if isfield(spk, 'filtr') && ~isempty(spk.filtr)
    tr = spk.filtr(:);
    nS = numel(tr);
    tTr = (0:nS-1)' / Fs;

    % decimate the trace for drawing only; a 70 s recording at 20 kHz is
    % 1.4e6 points, which makes an unreasonably heavy PNG. The V_min
    % markers are drawn from the stored values, so they are unaffected.
    step = max(1, ceil(nS / opt.MaxTracePlotPoints));

    figures(4) = figure('Color', 'w', 'Visible', visibility, 'Name', [fig_label ' | Filtered trace'], 'Position', [100 100 1400 380]);
    hold on; box off
    plot(tTr(1:step:end), tr(1:step:end), '-', 'Color', [0.55 0.55 0.55], ...
        'LineWidth', 0.5);
    for m = 1:nModes
        sel = (modeLabel == m) & ~isnan(VminT) & ~isnan(Vmin);
        if any(sel)
            plot(VminT(sel), Vmin(sel), 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor', cols(m,:), 'MarkerEdgeColor', 'none');
        end
    end
    un = isnan(modeLabel) & ~isnan(VminT) & ~isnan(Vmin);
    if any(un)
        plot(VminT(un), Vmin(un), 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
    end
    if nS > 1
        xlim([tTr(1) tTr(end)]);
    end
    xlabel('time (s)'); ylabel('filtered V');
    title(sprintf('%s  -  filtered trace (marker colors: V_{min} classes)', fig_label), 'Interpreter', 'tex');
end

if ~isempty(hazard)
    hasClassCurves = isfield(hazard, 'S_mode');
    decisionText = strsplit(monotonicityDecision, ' - ');
    monotonicityStatus = sprintf('Nondecreasing hazard null: %s', decisionText{1});

    % --- 5. pooled hazard with the weighted components shaded -------------
    figures(5) = figure('Color', 'w', 'Visible', visibility, 'Name', [fig_label ' | Pooled hazard']);
    hold on; box off
    if hasClassCurves
        % Stacked bands: each component is drawn between the running total
        % below it and the running total including it, so band m occupies
        % [sum_{k<m} h_k, sum_{k<=m} h_k] and the top edge is the pooled
        % hazard. patch rather than area because FaceAlpha on area objects is
        % unsupported in older versions and in Octave.
        Cc = hazard.h_contrib;
        Cc(isnan(Cc)) = 0;
        Ccum = cumsum(Cc, 2);
        lower = [zeros(size(Cc,1), 1), Ccum(:, 1:end-1)];
        tt = hazard.t(:);
        for m = 1:size(Ccum, 2)
            patch([tt; flipud(tt)], [Ccum(:,m); flipud(lower(:,m))], cols(m,:), ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end
    end

    % Restrict the pooled curve to the support of the reconstruction.
    hp = hazard.h_pooled;
    if isfield(hazard, 'h_recon')
        hp(isnan(hazard.h_recon)) = NaN;
    end
    plot(hazard.t, hp, 'k', 'LineWidth', 1.8);
    xlabel('ISI t (s)'); ylabel('h(t)  (s^{-1})');
    title({sprintf('%s - hazard', fig_label); monotonicityStatus}, 'Interpreter', 'tex');

    % --- 6. pooled survival with the per-class survivals dashed ------------
    figures(6) = figure('Color', 'w', 'Visible', visibility, 'Name', [fig_label ' | Survival']);
    hold on; box off
    survivalTitle = {sprintf('%s - survival', fig_label); monotonicityStatus};
    if hasClassCurves
        for m = 1:size(hazard.S_mode, 2)
            plot(hazard.t, hazard.S_mode(:,m), '--', 'Color', cols(m,:), 'LineWidth', 1.4);
        end
        survivalTitle{end+1} = 'Colored dashed curves: V_{min} classes';
    end
    plot(hazard.t, hazard.S_pooled, 'k', 'LineWidth', 1.8);
    ylim([0 1]);
    xlabel('ISI t (s)'); ylabel('S(t)');
    title(survivalTitle, 'Interpreter', 'tex');
end

if ~isempty(hazard) && hasClassCurves
    % --- 7. four-panel hazard decomposition -------------------------------
    H = hazard;
    t = H.t;
    nModes = size(H.h_mode, 2);


    tmax = max(t(~isnan(H.h_pooled)));
    if isempty(tmax) || ~isfinite(tmax) || tmax <= 0
        tmax = max([t; eps]);
    end

    modeLab = cell(1, nModes);
    for m = 1:nModes
        modeLab{m} = sprintf('class %d (p=%.2f)', m, H.p(m));
    end

    figures(7) = figure('Color', 'w', 'Visible', visibility, 'Name', [fig_label ' | Hazard decomposition | ' monotonicityStatus]);

    % --- (1) survival ------------------------------------------------------
    subplot(2,2,1); hold on; box off
    hSc = cell(1, nModes);
    for m = 1:nModes
        hSc{m} = plot(t, H.S_mode(:,m), '-', 'Color', cols(m,:), 'LineWidth', 1.1);
    end
    hS = [hSc{:}];
    hP = plot(t, H.S_pooled, 'k', 'LineWidth', 1.8);
    hR = plot(t, H.S_recon, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);
    xlim([0 tmax]); ylim([0 1]);
    xlabel('ISI t (s)'); ylabel('S(t)');
    title({sprintf('Survival (max|S_{recon}-S_{pooled}| = %.1e)', H.resid_S); ...
        monotonicityStatus}, 'Interpreter', 'tex');
    legend([hS hP hR], [modeLab, {'pooled', '\Sigma p_m S_m'}], ...
        'Location', 'northeast', 'Box', 'off');

    % --- (2) weights -------------------------------------------------------
    subplot(2,2,2); hold on; box off
    for m = 1:nModes
        plot(t, H.w(:,m), '-', 'Color', cols(m,:), 'LineWidth', 1.5);
        plot([0 tmax], [H.p(m) H.p(m)], ':', 'Color', cols(m,:), 'LineWidth', 1.0);
    end
    xlim([0 tmax]); ylim([0 1]);
    xlabel('ISI t (s)'); ylabel('w_m(t)');
    title('Class weights (dotted = constant p_m)', 'Interpreter', 'tex');

    % --- (3) weighted decomposition of the pooled hazard -------------------
    subplot(2,2,3); hold on; box off
    C = H.h_contrib;
    C(isnan(C)) = 0;
    hA = area(t, C);
    for m = 1:nModes
        set(hA(m), 'FaceColor', cols(m,:), 'EdgeColor', 'none');
        try
            set(hA(m), 'FaceAlpha', 0.55);   % not available in all versions
        catch
        end
    end
    hPool = plot(t, H.h_pooled, 'k', 'LineWidth', 1.8);
    xlim([0 tmax]);
    xlabel('ISI t (s)'); ylabel('h(t)  (s^{-1})');
    title({'Pooled hazard and stacked class contributions'; monotonicityStatus}, 'Interpreter', 'tex');
    legend([hA(:)' hPool], [modeLab, {'pooled'}], 'Location', 'northwest', 'Box', 'off');

    % --- (4) the three reconstructions -------------------------------------
    subplot(2,2,4); hold on; box off
    h1 = plot(t, H.h_pooled, 'k', 'LineWidth', 1.8);
    h2 = plot(t, H.h_recon,  '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4);
    h3 = plot(t, H.h_constw, '-.', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.2);
    h4 = plot(t, H.h_sum,    ':',  'Color', [0.45 0.45 0.45], 'LineWidth', 1.6);
    set(gca, 'YScale', 'log');
    xlim([0 tmax]);
    xlabel('ISI t (s)'); ylabel('h(t)  (s^{-1})');
    title('Reconstructions (rms rel. dev.)', 'Interpreter', 'tex');
    legend([h1 h2 h3 h4], {'pooled', ...
        sprintf('Summed class contributions: %.2f', H.resid_recon(2)), ...
        sprintf('\\Sigma p_m h_m: %.1f', H.resid_constw(2)), ...
        sprintf('\\Sigma h_m: %.1f', H.resid_sum(2))}, ...
        'Location', 'southeast', 'Box', 'off');
end

if ~opt.ShowFigures
    close(figures(isgraphics(figures)));
end
end
