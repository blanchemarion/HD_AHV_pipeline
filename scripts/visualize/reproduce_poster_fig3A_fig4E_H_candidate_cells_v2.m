function results = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2(varargin)
% REPRODUCE_POSTER_FIG3A_FIG4E_H_CANDIDATE_CELLS
% Reproduce the candidate-cell analyses underlying poster Figure 3A and
% Figure 4E-H, using every candidate neuron from every fish.
%
% The function:
%   1. Recursively finds the newest *_candidate_neurons_clean.mat file for
%      each recording in surface, pachon and molino.
%   2. Loads all candidate cells, then identifies phase-tuned cells within
%      each fish using the ORI_V15 circular-shift information test followed
%      by Benjamini-Hochberg FDR across that fish's candidate cells.
%   3. Aligns phase within every fish by a rotation only (default: the peak
%      of the fish-average tuning curve is placed at 0 degrees).
%   4. Recomputes a two-dimensional cell PCA from the candidate traces and
%      resolves its arbitrary rotation/reflection against the aligned
%      preferred phases. Scores are circle-centred and normalized within
%      fish before pooling.
%   5. Plots all individual HD tuning curves as phase-aligned heat maps.
%   6. Reproduces angular resolution (FWHM), consistency (MVL), directional
%      information and the poster's reliability/SNR panel.
%   7. Additionally calculates a true alternating-block split-half tuning
%      reliability for every cell.
%   8. Uses fish, not neurons, as the unit of morph-level inference.
%
% IMPORTANT ABOUT PHASE ALIGNMENT
% --------------------------------
% With calcium traces and network phase alone there is no external absolute
% heading landmark. The default therefore produces a RELATIVE phase: each
% fish's population-tuning peak is set to 0 degrees. This is appropriate for
% comparing tuning shapes and for pooling PCA rings, and it cannot change
% FWHM, MVL, information or SNR. If network_phase_rad was already placed in
% a common anatomical/behavioural convention upstream, use:
%
%   'PhaseAlignment', 'preserve'
%
% A phase sign cannot be inferred safely without anatomy or behaviour.
% Recording-specific sign or offset corrections may be supplied explicitly
% with PhaseSigns and PhaseOffsetsDeg (see examples below).
%
% BASIC USAGE
% -----------
%   R = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2;
%
% This searches:
%   /path/to/data/surface
%   /path/to/data/pachon
%   /path/to/data/molino
%
% Use an explicit root or explicit files:
%   R = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
%       'DataRoot', '/path/to/data', ...
%       'OutputDir', '/path/to/data/poster_reproduction');
%
%   R = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
%       'Files', { ...
%       '/path/to/data/pachon/rec06/rec06_pachon_hdMetrics_20260513_162704_candidate_neurons_clean.mat', ...
%       '/path/to/data/surface/rec102/rec102_surface_hdMetrics_20260519_155841_candidate_neurons_clean.mat'});
%
% Manual phase convention corrections are keyed by recording name:
%   signs.rec06 = -1;
%   offsets.rec06 = 30;
%   R = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
%       'PhaseSigns', signs, 'PhaseOffsetsDeg', offsets);
%
% Poster Figure 4E-H was labelled "Properties of phase-tuned cells". The
% default keeps all candidates for PCA/heatmaps and restricts the metric
% panels to cells passing the ORI_V15 criterion:
%   - Skaggs information about network phase;
%   - 500 circular-shift shuffles per cell (10-90% of the recording);
%   - Benjamini-Hochberg q < 0.05 within each fish.
%
% To show all candidates in the metric panels instead, use:
%   R = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
%       'MetricCellSubset', 'all_candidates');
%
% For a quick code check only, PhaseTuningNShuffles can be reduced. Use the
% default 500 shuffles for the reported analysis.
%
% Restrict the run to the original poster cohort with explicit Files, or:
%   'IncludeFish', {'rec06','rec102', ...}
%
% KEY OUTPUT FILES
% ----------------
%   Figure3A_aligned_PC_projection.(png|svg)
%   Individual_cell_HD_tuning_heatmaps.(png|svg)
%   Figure4E_H_candidate_cell_metrics.(png|svg)
%   True_split_half_reliability.(png|svg)
%   Phase_alignment_QC.(png|svg)
%
% ORI FEATURE DEFINITIONS
% -----------------------
% The phase-rate curves are ORI_V15 Step 17 rate_bins: rectified activity,
% 36 phase bins and three-bin circular smoothing. Directional information
% and vector strength follow Step 17. FWHM and peak/FWHM selectivity follow
% Step 19/local_circ_fwhm, including the baseline-relative half height and
% interpolated circular crossings. ORI_V15 Step 19 has no SNR variable;
% Figure 4H uses the poster export's peak divided by the mean of the lowest
% 20% of those same rate bins. Stored feature values are QC references only,
% except exact ORI information is retained when phase_ok was not exported.
%
% STATISTICAL UNIT
% ----------------
% Cell values are shown descriptively. Each fish is reduced to its median
% cell metric for morph comparisons. The two a-priori rank comparisons are
% Surface-Pachon and Surface-Molino; fish labels are permuted and the raw,
% unadjusted p-value is reported for each planned contrast. Cell-bootstrap
% confidence intervals show that fish with different
% neuron counts have different within-fish precision, without allowing fish
% with many cells to dominate the inferential test.
%
% MATLAB requirements: R2020b or newer is recommended. No Circular
% Statistics Toolbox is required.

%% Parse inputs
defaultCfg = pipeline_config();
p = inputParser;
p.FunctionName = mfilename;

addParameter(p, 'DataRoot', defaultCfg.DataRoot, @(x) ischar(x) || isstring(x));
addParameter(p, 'Files', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'OutputDir', defaultCfg.PosterFigureDir, @(x) ischar(x) || isstring(x));
addParameter(p, 'MorphOrder', {'Surface','Molino','Pachon'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'IncludeFish', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'ExcludeFish', {}, @(x) ischar(x) || isstring(x) || iscell(x));

addParameter(p, 'PhaseAlignment', 'population_peak', @(x) ischar(x) || isstring(x));
addParameter(p, 'PhaseReferenceDeg', 0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PhaseSigns', struct(), @(x) isstruct(x) || isa(x, 'containers.Map') || isempty(x));
addParameter(p, 'PhaseOffsetsDeg', struct(), @(x) isstruct(x) || isa(x, 'containers.Map') || isempty(x));

addParameter(p, 'NPhaseBins', 36, @(x) isnumeric(x) && isscalar(x) && x >= 12);
addParameter(p, 'ZScoreTracesForPCA', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ReliabilityBlockSec', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MinReliabilityBins', 12, @(x) isnumeric(x) && isscalar(x) && x >= 3);

addParameter(p, 'MetricCellSubset', 'phase_tuned', @(x) ischar(x) || isstring(x));
addParameter(p, 'PhaseSelectivityThreshold', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PhaseTuningFDRAlpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'PhaseTuningNShuffles', 500, @(x) isnumeric(x) && isscalar(x) && x >= 20);
addParameter(p, 'PhaseTuningSmoothBins', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'PhaseTuningMinShiftFraction', 0.10, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'PhaseTuningMaxShiftFraction', 0.90, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'MinimumMetricCellsPerFish', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'PlannedContrasts', {'Surface','Pachon';'Surface','Molino'}, ...
    @(x) iscell(x) && size(x,2) == 2);

addParameter(p, 'NBootstrap', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
addParameter(p, 'NPermutations', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
addParameter(p, 'RandomSeed', 1701, @(x) isnumeric(x) && isscalar(x));
% Deprecated compatibility option. Figures are always exported as PNG and SVG.
addParameter(p, 'SaveFigFiles', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Visible', defaultCfg.FigureVisible, @(x) ischar(x) || isstring(x));

parse(p, varargin{:});
cfg = p.Results;
cfg.DataRoot = char(cfg.DataRoot);
cfg.OutputDir = char(cfg.OutputDir);
cfg.PhaseAlignment = lower(char(cfg.PhaseAlignment));
cfg.MetricCellSubset = lower(char(cfg.MetricCellSubset));
cfg.MorphOrder = cellstr(string(cfg.MorphOrder));
cfg.IncludeFish = cellstr(string(cfg.IncludeFish));
cfg.ExcludeFish = cellstr(string(cfg.ExcludeFish));
cfg.Visible = char(cfg.Visible);

validPhaseModes = {'population_peak', 'preserve'};
if ~ismember(cfg.PhaseAlignment, validPhaseModes)
    error('PhaseAlignment must be ''population_peak'' or ''preserve''.');
end

validMetricSubsets = {'all_candidates','phase_tuned','sel_phi_threshold'};
if ~ismember(cfg.MetricCellSubset, validMetricSubsets)
    error(['MetricCellSubset must be ''all_candidates'', ''phase_tuned'', ' ...
        'or ''sel_phi_threshold''.']);
end
if strcmp(cfg.MetricCellSubset,'sel_phi_threshold') && ...
        ~isfinite(cfg.PhaseSelectivityThreshold)
    error(['PhaseSelectivityThreshold must be supplied when ' ...
        'MetricCellSubset is ''sel_phi_threshold''.']);
end
if cfg.PhaseTuningMaxShiftFraction <= cfg.PhaseTuningMinShiftFraction
    error('PhaseTuningMaxShiftFraction must exceed PhaseTuningMinShiftFraction.');
end

rng(cfg.RandomSeed, 'twister');

if isempty(cfg.OutputDir)
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    cfg.OutputDir = fullfile(pwd, ['poster_candidate_cell_reproduction_' stamp]);
end
if exist(cfg.OutputDir, 'dir') ~= 7
    mkdir(cfg.OutputDir);
end

fprintf('\n=== Poster Figure 3A / Figure 4E-H candidate-cell analysis ===\n');
fprintf('Output directory:\n  %s\n', cfg.OutputDir);
fprintf('Phase alignment mode: %s\n', cfg.PhaseAlignment);
fprintf(['Phase-tuned criterion: ORI_V15 information shuffle, %d shifts/cell, ' ...
    'within-fish BH-FDR q < %.3g\n'], round(cfg.PhaseTuningNShuffles), ...
    cfg.PhaseTuningFDRAlpha);

%% Discover one newest candidate file per fish
fileTable = discoverCandidateFiles(cfg);
if isempty(fileTable)
    error(['No *_candidate_neurons_clean.mat files were found. ' ...
        'Check DataRoot or pass an explicit Files list.']);
end

fprintf('Found %d fish files.\n', height(fileTable));
disp(fileTable(:, {'Morph','Fish','File'}));

%% Shared aligned phase grid
nBins = round(cfg.NPhaseBins);
commonEdges = linspace(-pi, pi, nBins + 1);
commonCenters = (commonEdges(1:end-1) + commonEdges(2:end)) ./ 2;
commonCentersDeg = rad2deg(commonCenters);

%% Load and analyse each fish
fishData = repmat(emptyFishStruct(), height(fileTable), 1);
cellTables = cell(height(fileTable), 1);
alignmentTables = cell(height(fileTable), 1);

for iFish = 1:height(fileTable)
    fprintf('\n[%d/%d] %s %s\n  %s\n', iFish, height(fileTable), ...
        fileTable.Morph(iFish), fileTable.Fish(iFish), fileTable.File(iFish));

    S = load(char(fileTable.File(iFish)));
    F = analyseOneFish(S, char(fileTable.Morph(iFish)), ...
        char(fileTable.Fish(iFish)), char(fileTable.File(iFish)), ...
        commonCenters, commonEdges, cfg);

    fishData(iFish) = F;
    cellTables{iFish} = makeCellTable(F);
    alignmentTables{iFish} = makeAlignmentTable(F);

    fprintf(['  %d candidate cells | phase offset %+0.1f deg | ' ...
        '%d phase-tuned | PCA-phase fit R = %.3f\n'], F.nCells, ...
        rad2deg(F.phaseOffsetRad), nnz(F.PhaseTuned), F.pcaPhaseFitR);
end

cellTable = vertcat(cellTables{:});
alignmentTable = vertcat(alignmentTables{:});
featureMethodAudit = buildORIFeatureMethodAudit(cellTable);

% Figure 3A and tuning heatmaps always retain all candidates. Figure 4E-H
% can reproduce the poster's phase-tuned-cell subset explicitly.
[metricCellTable, subsetDescription] = selectMetricCellSubset(cellTable, cfg);
fprintf('\nFigure 4E-H cell subset: %s\n', subsetDescription);
fprintf('Retained %d / %d candidate cells for metric panels.\n', ...
    height(metricCellTable), height(cellTable));

%% Fish-level summaries with neuron-bootstrap uncertainty
metricDefs = posterMetricDefinitions();
fishSummary = buildFishSummary(metricCellTable, metricDefs, cfg.NBootstrap);

%% Fish-level permutation statistics
statsMetricDefs = [metricDefs; trueReliabilityDefinition()];
statsTable = runMorphStatistics(fishSummary, statsMetricDefs, cfg);
legacyPooledStats = runLegacyPooledCellDiagnostic(metricCellTable, metricDefs, cfg);

%% Figures
morphColors = [ ...
    0.88 0.25 0.28;  % Surface: red
    0.25 0.72 0.38;  % Molino: green
    0.95 0.78 0.10]; % Pachon: yellow

figPaths = struct();
figPaths.figure3A = plotAlignedPCA(fishData, cfg, morphColors);
figPaths.tuningHeatmaps = plotTuningHeatmaps(fishData, commonCentersDeg, cfg);
figPaths.figure4EH = plotMetricPanels(metricCellTable, fishSummary, statsTable, ...
    metricDefs, cfg, morphColors, 'poster');
figPaths.trueReliability = plotMetricPanels(metricCellTable, fishSummary, statsTable, ...
    trueReliabilityDefinition(), cfg, morphColors, 'reliability');
figPaths.phaseQC = plotPhaseAlignmentQC(fishData, commonCentersDeg, cfg, morphColors);

%% Return complete results in memory (no MAT/CSV files are written)
results = struct();
results.config = cfg;
results.inputFiles = fileTable;
results.cells = cellTable;
results.metricCells = metricCellTable;
results.metricCellSubsetDescription = subsetDescription;
results.fishSummary = fishSummary;
results.statistics = statsTable;
results.legacyPooledCellStatisticsDiagnosticOnly = legacyPooledStats;
results.phaseAlignment = alignmentTable;
results.ORIFeatureMethodAudit = featureMethodAudit;
results.phaseBinCentersDeg = commonCentersDeg;
results.fish = fishData;
results.figurePaths = figPaths;

fprintf('\nAnalysis complete.\n');
fprintf(['All candidates were retained for PCA/heatmaps; the ORI_V15 ' ...
    'phase-tuned subset was used for Figure 4E-H by default.\n']);
fprintf('Fish was the inferential unit; planned p-values are unadjusted.\n');
fprintf('PNG figures saved in:\n  %s\n', cfg.OutputDir);

end

%% ========================================================================
function fileTable = discoverCandidateFiles(cfg)
% Return one newest matching file per morph + recording.

if isempty(cfg.Files)
    allPaths = {};
    for i = 1:numel(cfg.MorphOrder)
        morphLower = lower(cfg.MorphOrder{i});
        root = fullfile(cfg.DataRoot, morphLower);
        if exist(root, 'dir') ~= 7
            warning('Morph directory not found: %s', root);
            continue;
        end
        d = dir(fullfile(root, '**', '*_candidate_neurons_clean.mat'));
        for k = 1:numel(d)
            allPaths{end+1,1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
        end
    end
else
    allPaths = cellstr(string(cfg.Files(:)));
end

if isempty(allPaths)
    fileTable = table();
    return;
end

n = numel(allPaths);
morph = strings(n,1);
fish = strings(n,1);
score = nan(n,1);

for i = 1:n
    if exist(allPaths{i}, 'file') ~= 2
        error('Input file does not exist: %s', allPaths{i});
    end
    [morph(i), fish(i)] = inferMorphAndFish(allPaths{i});
    score(i) = candidateFileTimeScore(allPaths{i});
end

% Keep newest file for duplicate analyses of the same recording.
key = lower(morph + "|" + fish);
uniqueKeys = unique(key, 'stable');
keep = false(n,1);
for i = 1:numel(uniqueKeys)
    idx = find(key == uniqueKeys(i));
    [~, j] = max(score(idx));
    keep(idx(j)) = true;
    if numel(idx) > 1
        fprintf('  Duplicate candidate files for %s; keeping newest:\n    %s\n', ...
            uniqueKeys(i), allPaths{idx(j)});
    end
end

allPaths = allPaths(keep);
morph = morph(keep);
fish = fish(keep);
score = score(keep);

% Optional cohort restriction. Entries may be "rec06" or
% "Pachon/rec06". Explicit Files remains the most reproducible option.
cohortKeep = true(numel(fish),1);
if ~isempty(cfg.IncludeFish)
    cohortKeep(:) = false;
    for i = 1:numel(fish)
        cohortKeep(i) = fishMatchesList(morph(i), fish(i), cfg.IncludeFish);
    end
end
if ~isempty(cfg.ExcludeFish)
    for i = 1:numel(fish)
        if fishMatchesList(morph(i), fish(i), cfg.ExcludeFish)
            cohortKeep(i) = false;
        end
    end
end
allPaths = allPaths(cohortKeep);
morph = morph(cohortKeep);
fish = fish(cohortKeep);
score = score(cohortKeep);

orderIndex = nan(numel(morph),1);
for i = 1:numel(cfg.MorphOrder)
    orderIndex(strcmpi(morph, cfg.MorphOrder{i})) = i;
end
orderIndex(~isfinite(orderIndex)) = numel(cfg.MorphOrder) + 1;

fileTable = table(morph, fish, string(allPaths), score, orderIndex, ...
    'VariableNames', {'Morph','Fish','File','FileTimeScore','MorphOrder'});
fileTable = sortrows(fileTable, {'MorphOrder','Fish'});
fileTable.MorphOrder = [];
end

function tf = fishMatchesList(morph, fish, list)
tf = false;
fishOnly = lower(strtrim(string(fish)));
combinedSlash = lower(strtrim(string(morph) + "/" + string(fish)));
combinedPipe = lower(strtrim(string(morph) + "|" + string(fish)));
entries = lower(strtrim(string(list(:))));
tf = any(entries == fishOnly | entries == combinedSlash | entries == combinedPipe);
end

function [morph, fish] = inferMorphAndFish(filePath)
[~, name] = fileparts(filePath);
low = lower(filePath);

if ~isempty(strfind(low, 'surface')) %#ok<STREMP>
    morph = "Surface";
elseif ~isempty(strfind(low, 'pachon')) %#ok<STREMP>
    morph = "Pachon";
elseif ~isempty(strfind(low, 'molino')) %#ok<STREMP>
    morph = "Molino";
else
    error('Cannot infer morph from path: %s', filePath);
end

tok = regexp(name, '^(rec[^_]*)_', 'tokens', 'once');
if isempty(tok)
    [folder, ~, ~] = fileparts(filePath);
    [~, fishName] = fileparts(folder);
    fish = string(fishName);
else
    fish = string(tok{1});
end
end

function score = candidateFileTimeScore(filePath)
[~, name, ~] = fileparts(filePath);
tok = regexp(name, '_(\d{8})_(\d{6})_candidate_neurons_clean$', 'tokens', 'once');
if ~isempty(tok)
    try
        score = datenum([tok{1} tok{2}], 'yyyymmddHHMMSS'); %#ok<DATNM>
        return;
    catch
    end
end
d = dir(filePath);
score = d.datenum;
end

%% ========================================================================
function F = analyseOneFish(S, morph, fish, filePath, commonCenters, commonEdges, cfg)

required = {'candidate_cell_ids','calcium_traces','network_phase_rad'};
for i = 1:numel(required)
    if ~isfield(S, required{i})
        error('%s is missing required variable "%s".', filePath, required{i});
    end
end

cellIds = double(S.candidate_cell_ids(:));
nCells = numel(cellIds);
traces = orientTraceMatrix(double(S.calcium_traces), nCells, filePath);
phase = double(S.network_phase_rad(:));

nTime = min(size(traces,1), numel(phase));
if size(traces,1) ~= numel(phase)
    warning('%s: traces and phase lengths differ; truncating both to %d.', fish, nTime);
end
traces = traces(1:nTime,:);
phase = wrapToPiLocal(phase(1:nTime));

if isfield(S, 'time_s') && numel(S.time_s) >= nTime
    time = double(S.time_s(1:nTime));
    time = time(:);
else
    fs = 3.41;
    if isfield(S, 'metadata') && isfield(S.metadata, 'default_sampling_rate_hz')
        fsCandidate = double(S.metadata.default_sampling_rate_hz);
        if isfinite(fsCandidate) && fsCandidate > 0
            fs = fsCandidate;
        end
    end
    time = (0:nTime-1)' ./ fs;
end

% Stored values are retained only for QC. Primary metrics are recomputed
% below from the exact ORI_V15 Step 17 rate curves.
features = extractStoredFeatures(S, cellIds);

% Recreate ORI_V15 Step 17 on ALL candidate cells. A cell is retained when
% its Skaggs information exceeds its circular-shift null after BH-FDR across
% the candidate cells from this fish. The selection is recomputed here and
% does not depend on an arbitrary sel_phi threshold or a stored flag.
[phaseOk, phaseValiditySource] = extractPhaseValidityMask(S, nTime);
if strcmp(phaseValiditySource, 'all finite phase frames (phase_ok unavailable)')
    warning(['%s: phase_ok was not exported. The local ORI shuffle selection ' ...
        'uses all finite phase frames; canonical stored ORI curves/information ' ...
        'are still used for feature calculation when available.'], fish);
end
phaseTest = oriPhaseTuningShuffleTest(traces, phase, phaseOk, cfg, fish);

% The candidate export contains ORI_V15 Step 17 rate_bins. Use those
% canonical curves when available; otherwise use the curves just recomputed
% by the identical Step 17 implementation above.
[curves, originalCenters, phaseCurveSource] = ...
    extractCanonicalORIPhaseCurves(S, nCells, phaseTest);

% Step 19 circular FWHM/selectivity and Step 17 vector strength are always
% recomputed from the same canonical rate curves.
[fwhmRad, peakPhi, phaseSelectivity] = ...
    oriCircularFWHMAndSelectivity(curves, originalCenters);
fwhm = rad2deg(fwhmRad);
[preferredOriginal, mvl] = oriPreferredPhaseAndMVL(curves, originalCenters);

% Directional information requires the exact Step 17 occupancy mask. The
% candidate export does not contain phase_ok, but does contain ORI's exact
% info_bits. Use that canonical value when present; otherwise use the full
% local Step 17 recomputation and expose both values in the output table.
infoBits = phaseTest.InformationBits;
hasStoredORIInformation = isfinite(features.Information);
infoBits(hasStoredORIInformation) = features.Information(hasStoredORIInformation);

% ORI_V15 Step 19 does not define an SNR variable. Figure 4H uses the
% poster-specific peak/baseline measure recovered from the candidate-file
% export: peak divided by the mean of the lowest 20% of phase bins. It is
% recomputed here on the ORI Step 17 rate curves so preprocessing and
% smoothing remain identical for all four displayed metrics.
snr = peakToBaseline(curves);

phaseSign = lookupRecordingValue(cfg.PhaseSigns, fish, 1);
if ~ismember(phaseSign, [-1 1])
    error('PhaseSigns.%s must be +1 or -1.', fish);
end
manualOffset = deg2rad(lookupRecordingValue(cfg.PhaseOffsetsDeg, fish, 0));

switch cfg.PhaseAlignment
    case 'population_peak'
        populationCurve = mean(curves, 1, 'omitnan');
        [~, iPeak] = max(populationCurve);
        populationPeak = originalCenters(iPeak);
        phaseOffset = deg2rad(cfg.PhaseReferenceDeg) - phaseSign .* populationPeak + manualOffset;
    case 'preserve'
        phaseOffset = manualOffset;
    otherwise
        error('Unsupported phase alignment mode.');
end
phaseOffset = wrapToPiLocal(phaseOffset);

preferredAligned = wrapToPiLocal(phaseSign .* preferredOriginal + phaseOffset);

% Evaluate each original curve on the shared aligned phase grid.
queryOriginal = phaseSign .* (commonCenters - phaseOffset);
alignedCurves = periodicInterpolateCurves(originalCenters, curves, queryOriginal);

% True reliability: alternating contiguous time blocks, not adjacent frames.
splitReliability = splitHalfReliability(traces, phase, time, originalCenters, ...
    cfg.ReliabilityBlockSec, cfg.MinReliabilityBins);

% Recompute PCA positions from candidate traces, then remove PCA's arbitrary
% reflection and rotation by matching cell angle to preferred network phase.
[pcaXY, alignedPcaXY, handedness, rotationRad, fitR] = ...
    alignedCellPCA(traces, preferredAligned, mvl, logical(cfg.ZScoreTracesForPCA));

F = emptyFishStruct();
F.morph = string(morph);
F.fish = string(fish);
F.file = string(filePath);
F.nCells = nCells;
F.cellIds = cellIds;
F.phaseSign = phaseSign;
F.phaseOffsetRad = phaseOffset;
F.originalPhaseCentersRad = originalCenters;
F.phaseCurveSource = string(phaseCurveSource);
F.phaseValiditySource = string(phaseValiditySource);
F.alignedPhaseCentersRad = commonCenters;
F.tuningCurvesOriginal = curves;
F.tuningCurvesAligned = alignedCurves;
F.preferredPhaseOriginalRad = preferredOriginal;
F.preferredPhaseAlignedRad = preferredAligned;
F.FWHMDeg = fwhm;
F.MVL = mvl;
F.InformationBits = infoBits;
F.PeakToBaselineSNR = snr;
F.PeakPhaseRate = peakPhi;
F.PhaseSelectivity = phaseSelectivity;
F.PhaseTuningInformationBits = phaseTest.InformationBits;
F.PhaseTuningPValue = phaseTest.PValue;
F.PhaseTuningQValue = phaseTest.QValue;
F.PhaseTuned = phaseTest.IsTuned;
F.PhaseTuningNValidFrames = phaseTest.NValidFrames;
F.StoredPhaseTuned = features.PhaseTuned;
F.StoredFWHMDeg = features.FWHM;
F.StoredMVL = features.MVL;
F.StoredInformationBits = features.Information;
F.StoredPeakToBaselineSNR = features.SNR;
F.StoredPhaseSelectivity = features.PhaseSelectivity;
F.InformationFromStoredORI = hasStoredORIInformation;
F.SplitHalfReliability = splitReliability;
F.pcaXY = pcaXY;
F.alignedPcaXY = alignedPcaXY;
F.pcaHandedness = handedness;
F.pcaRotationRad = rotationRad;
F.pcaPhaseFitR = fitR;
end

function F = emptyFishStruct()
F = struct( ...
    'morph', string.empty, 'fish', string.empty, 'file', string.empty, ...
    'nCells', 0, 'cellIds', [], 'phaseSign', 1, 'phaseOffsetRad', NaN, ...
    'phaseCurveSource', string.empty, 'phaseValiditySource', string.empty, ...
    'originalPhaseCentersRad', [], 'alignedPhaseCentersRad', [], ...
    'tuningCurvesOriginal', [], 'tuningCurvesAligned', [], ...
    'preferredPhaseOriginalRad', [], 'preferredPhaseAlignedRad', [], ...
    'FWHMDeg', [], 'MVL', [], 'InformationBits', [], ...
    'PeakToBaselineSNR', [], 'PeakPhaseRate', [], 'PhaseSelectivity', [], ...
    'PhaseTuningInformationBits', [], 'PhaseTuningPValue', [], ...
    'PhaseTuningQValue', [], 'PhaseTuned', [], 'StoredPhaseTuned', [], ...
    'PhaseTuningNValidFrames', 0, 'StoredFWHMDeg', [], 'StoredMVL', [], ...
    'StoredInformationBits', [], 'StoredPeakToBaselineSNR', [], ...
    'StoredPhaseSelectivity', [], ...
    'InformationFromStoredORI', [], ...
    'SplitHalfReliability', [], ...
    'pcaXY', [], 'alignedPcaXY', [], 'pcaHandedness', NaN, ...
    'pcaRotationRad', NaN, 'pcaPhaseFitR', NaN);
end

function traces = orientTraceMatrix(traces, nCells, filePath)
if size(traces,2) == nCells
    return;
elseif size(traces,1) == nCells
    traces = traces';
else
    error('%s: no calcium_traces dimension matches %d candidate IDs.', filePath, nCells);
end
end

function [phaseOk, source] = extractPhaseValidityMask(S, nTime)
% ORI_V15 uses phase_ok when available and otherwise accepts every finite
% phase sample. Candidate-neuron exports do not always carry this vector.
phaseOk = true(nTime,1);
source = 'all finite phase frames (phase_ok unavailable)';
aliases = {'phase_ok','network_phase_ok','phase_valid','network_phase_valid'};
for i = 1:numel(aliases)
    if ~isfield(S, aliases{i})
        continue;
    end
    candidate = S.(aliases{i});
    if (isnumeric(candidate) || islogical(candidate)) && numel(candidate) >= nTime
        phaseOk = logical(candidate(1:nTime));
        phaseOk = phaseOk(:);
        source = aliases{i};
        return;
    end
end
end

function out = oriPhaseTuningShuffleTest(traces, phase, phaseOk, cfg, fish)
% Faithful, self-contained implementation of ORI_V15 Step 17:
%  1) clip each frame at its 2nd/98th cell percentiles and frame-center;
%  2) rectify activity;
%  3) calculate 36-bin, 3-bin-smoothed Skaggs phase information;
%  4) compare each cell with circular time shifts spanning 10-90%;
%  5) apply Benjamini-Hochberg FDR across candidate cells within this fish.

nBins = round(cfg.NPhaseBins);
smoothBins = round(cfg.PhaseTuningSmoothBins);
nShuffles = round(cfg.PhaseTuningNShuffles);

Fz = clipCenterFramewiseLocal(double(traces), 2, 98);
valid = isfinite(phase(:)) & logical(phaseOk(:));
phi = wrapToPiLocal(double(phase(valid)));
activity = max(Fz(valid,:), 0);

if numel(phi) < 20
    error('%s: fewer than 20 valid phase frames for the ORI_V15 tuning test.', fish);
end

edges = linspace(-pi, pi, nBins+1);
centers = ((edges(1:end-1) + edges(2:end))./2)';
binIndex = discretize(phi, edges);
hasBin = isfinite(binIndex);
binIndex = binIndex(hasBin);
activity = activity(hasBin,:);

occupancyCounts = accumarray(binIndex, 1, [nBins,1], @sum, 0);
occupancy = occupancyCounts ./ max(sum(occupancyCounts),1);
nValid = size(activity,1);
binAverager = sparse(binIndex, (1:nValid)', ...
    1 ./ occupancyCounts(binIndex), nBins, nValid);

nCells = size(activity,2);
rateBins = circularSmoothPhaseCurves(binAverager*activity, nBins, smoothBins);
[observed, preferredPhase, meanVectorLength] = ...
    oriInformationAndVectorStrength(rateBins, occupancy, centers);
pValue = nan(nCells,1);

fprintf(['  ORI_V15 phase-tuning selection for %s: %d cells x %d ' ...
    'circular shifts...\n'], fish, nCells, nShuffles);

minShift = max(1, floor(cfg.PhaseTuningMinShiftFraction .* nValid));
maxShift = max(1, floor(cfg.PhaseTuningMaxShiftFraction .* nValid));
if maxShift <= minShift
    error('%s: invalid circular-shift range for the phase-tuning test.', fish);
end

timeIndex = (1:nValid)';
shuffleBatchSize = min(100,nShuffles);
for k = 1:nCells
    r0 = activity(:,k);

    nAtLeastObserved = 0;
    for firstShuffle = 1:shuffleBatchSize:nShuffles
        thisBatch = min(shuffleBatchSize, nShuffles-firstShuffle+1);
        shifts = randi([minShift maxShift], 1, thisBatch);
        % These indices reproduce circshift(r0,shift) column by column.
        shiftedIndex = mod(timeIndex - shifts - 1, nValid) + 1;
        shiftedActivity = r0(shiftedIndex);
        shuffledCurves = binAverager * shiftedActivity;
        shuffledInformation = phaseInformationFromBinnedCurves( ...
            shuffledCurves, occupancy, nBins, smoothBins);
        nAtLeastObserved = nAtLeastObserved + ...
            sum(shuffledInformation >= observed(k));
    end
    % ORI_V15 uses the direct empirical tail fraction (zero is therefore
    % possible with 500 shifts); preserve that definition exactly.
    pValue(k) = nAtLeastObserved ./ nShuffles;
end

qValue = benjaminiHochbergLocal(pValue);
isTuned = isfinite(qValue) & qValue < cfg.PhaseTuningFDRAlpha;

fprintf('  ORI_V15 phase tuned after within-fish BH-FDR: %d/%d cells.\n', ...
    nnz(isTuned), nCells);

out = struct();
out.InformationBits = observed;
out.PValue = pValue;
out.QValue = qValue;
out.IsTuned = isTuned;
out.NValidFrames = nValid;
out.Occupancy = occupancy;
out.PhaseCentersRad = centers;
out.RateBins = rateBins;
out.PreferredPhaseRad = preferredPhase;
out.MVL = meanVectorLength;
end

function Fz = clipCenterFramewiseLocal(F, pLow, pHigh)
% Same preprocessing as ORI_V15/clip_center_framewise: percentiles are
% taken across cells within each frame, then the clipped frame mean is
% subtracted from every cell.
lower = prctile(F, pLow, 2);
upper = prctile(F, pHigh, 2);
Fz = min(F, upper);
Fz = max(Fz, lower);
Fz = Fz - mean(Fz,2);
end

function information = phaseInformationFromBinnedCurves(curve, occupancy, nBins, smoothBins)
% curve is bins x one-or-more traces. Vectorising the shuffle batches gives
% the same result as ORI_V15's nested circshift/accumarray loop but is much
% faster for a full multi-fish dataset.
curve = circularSmoothPhaseCurves(curve, nBins, smoothBins);
if size(curve,1) ~= nBins
    error('Internal phase-tuning error: expected %d phase bins, obtained %d.', ...
        nBins, size(curve,1));
end

occupied = occupancy > 0;
meanRate = sum(occupancy(occupied) .* curve(occupied,:), 1);
ratio = curve(occupied,:) ./ max(meanRate,realmin('double'));
information = sum(occupancy(occupied) .* ratio .* log2(ratio + eps), 1);
information(~isfinite(meanRate) | meanRate <= 0) = 0;
end

function curve = circularSmoothPhaseCurves(curve, nBins, smoothBins)
% Exact matrix form of ORI_V15's three-copy circular convolution.
if smoothBins > 1
    extended = [curve; curve; curve];
    smoothed = conv2(extended, ones(smoothBins,1)./smoothBins, 'same');
    curve = smoothed((nBins+1):(2*nBins), :);
end
end

function [information, preferredPhase, meanVectorLength] = ...
        oriInformationAndVectorStrength(rateBins, occupancy, centers)
% ORI_V15 Step 17 definitions, evaluated for all cells at once.
occupied = occupancy > 0;
rates = rateBins(occupied,:);
p = occupancy(occupied);

meanRate = sum(p .* rates, 1);
ratio = rates ./ max(meanRate,realmin('double'));
information = sum(p .* ratio .* log2(ratio + eps), 1)';

vector = sum(rates .* exp(1i.*centers(occupied)), 1);
totalRate = sum(rates,1);
hasRate = totalRate > 0 & isfinite(totalRate);
vector(hasRate) = vector(hasRate) ./ totalRate(hasRate);
vector(~hasRate) = 0;
preferredPhase = angle(vector)';
meanVectorLength = abs(vector)';

noMeanRate = ~isfinite(meanRate) | meanRate <= 0;
information(noMeanRate) = 0;
preferredPhase(noMeanRate) = NaN;
meanVectorLength(noMeanRate) = 0;
end

function q = benjaminiHochbergLocal(p)
% Equivalent to ORI_V15's mafdr(...,'BHFDR',true) fallback.
p = p(:);
q = nan(size(p));
good = isfinite(p);
pv = p(good);
if isempty(pv)
    return;
end
[sorted, order] = sort(pv, 'ascend');
m = numel(sorted);
qSorted = sorted .* m ./ (1:m)';
for i = m-1:-1:1
    qSorted(i) = min(qSorted(i), qSorted(i+1));
end
qSorted = min(qSorted,1);
qUnsorted = nan(m,1);
qUnsorted(order) = qSorted;
q(good) = qUnsorted;
end

function [curves, centers, source] = ...
        extractCanonicalORIPhaseCurves(S, nCells, phaseTest)
% Candidate exports contain ORI_V15/hdMetrics.hd_single.phase.rate_bins.
% Preserve those exact curves when present because phase_ok is not always
% exported. The local Step 17 result is the principled fallback.
if isfield(S,'tuning_curves_phi_all_cells')
    curves = double(S.tuning_curves_phi_all_cells);
    if size(curves,1) == nCells
        % Cells x bins.
    elseif size(curves,2) == nCells
        curves = curves';
    else
        error('No tuning_curves_phi_all_cells dimension matches candidate count.');
    end

    centers = [];
    if isfield(S,'phi_bin_centers_rad')
        centers = double(S.phi_bin_centers_rad(:))';
    elseif isfield(S,'phi_bin_centers_deg')
        centers = deg2rad(double(S.phi_bin_centers_deg(:))');
    end
    if isempty(centers)
        nCurveBins = size(curves,2);
        edges = linspace(-pi, pi, nCurveBins + 1);
        centers = (edges(1:end-1) + edges(2:end)) ./ 2;
    end

    if numel(centers) ~= size(curves,2)
        error('Number of phase-bin centers does not match tuning-curve columns.');
    end
    source = 'stored ORI_V15 Step 17 rate_bins';
else
    curves = phaseTest.RateBins';
    centers = phaseTest.PhaseCentersRad';
    source = 'locally recomputed ORI_V15 Step 17 rate_bins';
end

centers = wrapToPiLocal(centers(:)');
[centers, order] = sort(centers);
curves = curves(:,order);
end

function features = extractStoredFeatures(S, cellIds)
n = numel(cellIds);
features = struct('FWHM', nan(n,1), 'MVL', nan(n,1), ...
    'Information', nan(n,1), 'SNR', nan(n,1), ...
    'PhaseSelectivity', nan(n,1), 'PhasePValue', nan(n,1), ...
    'PhaseTuned', nan(n,1));

if ~isfield(S, 'candidate_features')
    return;
end
T = S.candidate_features;

fwhm = getFeatureVector(T, {'FWHM_phi_deg','FWHM_deg','FWHM'}, n);
mvl = getFeatureVector(T, {'vecR','MVL','mvl','mean_vector_length'}, n);
info = getFeatureVector(T, {'info_bits','information_bits','Information'}, n);
snr = getFeatureVector(T, {'peak_to_baseline','SNR','snr'}, n);
phaseSelectivity = getFeatureVector(T, {'sel_phi','phase_selectivity','selectivity_phi'}, n);
phasePValue = getFeatureVector(T, {'p_phi','p_phase','phase_p_value','phase_pvalue', ...
    'pval_phi','phase_tuning_p'}, n);
phaseTuned = getFeatureVector(T, {'is_phase_tuned','phase_tuned','phaseTuned', ...
    'significant_phi','is_tuned_phi'}, n);
featureIds = getFeatureVector(T, {'candidate_cell_id','candidate_cell_ids','cell_id','roi_id'}, n);

if all(isfinite(featureIds)) && numel(unique(featureIds)) == n
    [tf, loc] = ismember(cellIds, featureIds);
    if all(tf)
        fwhm = fwhm(loc);
        mvl = mvl(loc);
        info = info(loc);
        snr = snr(loc);
        phaseSelectivity = phaseSelectivity(loc);
        phasePValue = phasePValue(loc);
        phaseTuned = phaseTuned(loc);
    end
end

features.FWHM = fwhm(:);
features.MVL = mvl(:);
features.Information = info(:);
features.SNR = snr(:);
features.PhaseSelectivity = phaseSelectivity(:);
features.PhasePValue = phasePValue(:);
features.PhaseTuned = phaseTuned(:);
end

function v = getFeatureVector(T, aliases, nExpected)
v = nan(nExpected,1);

if istable(T)
    names = T.Properties.VariableNames;
    for i = 1:numel(aliases)
        idx = find(strcmpi(names, aliases{i}), 1);
        if ~isempty(idx)
            candidate = T.(names{idx});
            if isnumeric(candidate) || islogical(candidate)
                candidate = double(candidate(:));
                if numel(candidate) == nExpected
                    v = candidate;
                end
            end
            return;
        end
    end
elseif isstruct(T)
    names = fieldnames(T);
    for i = 1:numel(aliases)
        idx = find(strcmpi(names, aliases{i}), 1);
        if ~isempty(idx)
            candidate = T.(names{idx});
            candidate = double(candidate(:));
            if numel(candidate) == nExpected
                v = candidate;
            end
            return;
        end
    end
end
end

%% ========================================================================
function [fwhmRad, peakRate, selectivity] = ...
        oriCircularFWHMAndSelectivity(curves, centers)
% Exact ORI_V15 Step 19 wrapper around local_circ_fwhm. Width is returned
% in radians; selectivity follows ORI: peak rate divided by FWHM radians.
nCells = size(curves,1);
if numel(centers) ~= size(curves,2)
    error('ORI FWHM: phase centres and curve columns do not match.');
end

fwhmRad = nan(nCells,1);
peakRate = nan(nCells,1);
selectivity = nan(nCells,1);
for k = 1:nCells
    [fwhmRad(k), peakRate(k)] = ...
        oriLocalCircularFWHM(centers(:), curves(k,:)');
    if isfinite(fwhmRad(k)) && fwhmRad(k) > 0
        selectivity(k) = peakRate(k) ./ fwhmRad(k);
    end
end
end

function [preferredPhase, meanVectorLength] = ...
        oriPreferredPhaseAndMVL(curves, centers)
% Exact ORI_V15 Step 17 vector definition. The already rectified and
% smoothed rate bins are used directly, without occupancy weighting.
rates = curves';
vector = sum(rates .* exp(1i.*centers(:)), 1);
totalRate = sum(rates,1);
hasRate = totalRate > 0 & isfinite(totalRate);
vector(hasRate) = vector(hasRate) ./ totalRate(hasRate);
vector(~hasRate) = 0;
preferredPhase = angle(vector)';
meanVectorLength = abs(vector)';
preferredPhase(~hasRate) = NaN;
meanVectorLength(~hasRate) = 0;
end

function [widthRad, peakValue] = oriLocalCircularFWHM(centers, rates)
% Verbatim methodology of ORI_V15/local_circ_fwhm: circular interpolation
% around the main peak, at half the baseline-to-peak amplitude.
theta = centers(:);
response = rates(:);
if numel(theta) ~= numel(response)
    error('ORI FWHM: centres and rates must have the same length.');
end
if all(~isfinite(response)) || ~any(isfinite(response))
    widthRad = NaN;
    peakValue = NaN;
    return;
end

response(~isfinite(response)) = NaN;
peakValue = max(response,[],'omitnan');
if ~isfinite(peakValue) || peakValue <= 0
    widthRad = NaN;
    return;
end

baseline = min(response,[],'omitnan');
halfHeight = baseline + (peakValue-baseline)./2;
n = numel(theta);
thetaExtended = [theta-2*pi; theta; theta+2*pi];
responseExtended = [response; response; response];

[~, centralPeakIndex] = max(response);
centralPeakIndex = centralPeakIndex + n;

leftIndex = centralPeakIndex;
while leftIndex > 1 && responseExtended(leftIndex) >= halfHeight
    leftIndex = leftIndex-1;
end
if leftIndex < 1 || leftIndex+1 > numel(responseExtended) || ...
        ~isfinite(responseExtended(leftIndex)) || ...
        ~isfinite(responseExtended(leftIndex+1)) || ...
        responseExtended(leftIndex) == responseExtended(leftIndex+1)
    widthRad = NaN;
    return;
end
leftCrossing = interp1(responseExtended(leftIndex:leftIndex+1), ...
    thetaExtended(leftIndex:leftIndex+1), halfHeight, 'linear');

rightIndex = centralPeakIndex;
nExtended = numel(responseExtended);
while rightIndex < nExtended && responseExtended(rightIndex) >= halfHeight
    rightIndex = rightIndex+1;
end
if rightIndex > nExtended || rightIndex-1 < 1 || ...
        ~isfinite(responseExtended(rightIndex-1)) || ...
        ~isfinite(responseExtended(rightIndex)) || ...
        responseExtended(rightIndex-1) == responseExtended(rightIndex)
    widthRad = NaN;
    return;
end
rightCrossing = interp1(responseExtended(rightIndex-1:rightIndex), ...
    thetaExtended(rightIndex-1:rightIndex), halfHeight, 'linear');

widthRad = rightCrossing-leftCrossing;
if ~isfinite(widthRad) || widthRad < 0
    widthRad = NaN;
elseif widthRad > 2*pi
    widthRad = mod(widthRad,2*pi);
end
end

function snr = peakToBaseline(curves)
% Poster definition recovered from the supplied candidate file: maximum
% curve value divided by the mean of the lowest 20% of phase bins.
nLow = max(1, round(0.20 .* size(curves,2)));
snr = nan(size(curves,1),1);
for i = 1:size(curves,1)
    y = curves(i,isfinite(curves(i,:)));
    if numel(y) < 3
        continue;
    end
    y = sort(y);
    baseline = mean(y(1:min(nLow,numel(y))));
    if baseline > 0
        snr(i) = max(y) ./ baseline;
    end
end
end

function edges = circularEdgesFromCenters(centers)
centers = unwrap(centers(:)');
if numel(centers) < 2
    error('At least two phase centres are needed.');
end
mid = (centers(1:end-1) + centers(2:end)) ./ 2;
first = centers(1) - (mid(1) - centers(1));
last = centers(end) + (centers(end) - mid(end));
edges = [first mid last];
end

function curves = tuningCurveFromSamples(traces, phase, edges)
nBins = numel(edges)-1;
nCells = size(traces,2);
curves = nan(nCells,nBins);
[~,~,bin] = histcounts(phase, edges);
for b = 1:nBins
    use = bin == b;
    if any(use)
        curves(:,b) = mean(traces(use,:), 1, 'omitnan')';
    end
end
end

function aligned = periodicInterpolateCurves(centers, curves, query)
centers = unwrap(centers(:)');
query = wrapToPiLocal(query(:)');
x = [centers-2*pi, centers, centers+2*pi];
y = [curves, curves, curves];
aligned = interp1(x', y', query', 'linear')';
end

%% ========================================================================
function reliability = splitHalfReliability(traces, phase, time, centers, blockSec, minBins)
edges = circularEdgesFromCenters(centers);
block = floor((time - time(1)) ./ blockSec);
useA = mod(block,2) == 0;
useB = ~useA;

curveA = tuningCurveFromSamples(traces(useA,:), phase(useA), edges);
curveB = tuningCurveFromSamples(traces(useB,:), phase(useB), edges);

nCells = size(traces,2);
reliability = nan(nCells,1);
for i = 1:nCells
    good = isfinite(curveA(i,:)) & isfinite(curveB(i,:));
    if sum(good) < minBins
        continue;
    end
    a = curveA(i,good);
    b = curveB(i,good);
    a = a - mean(a);
    b = b - mean(b);
    denom = sqrt(sum(a.^2) .* sum(b.^2));
    if denom > 0
        reliability(i) = sum(a .* b) ./ denom;
    end
end
end

%% ========================================================================
function [rawXY, alignedXY, handedness, rotationRad, fitR] = ...
    alignedCellPCA(traces, preferredPhase, weights, doZscore)

nCells = size(traces,2);
rawXY = nan(nCells,2);
alignedXY = nan(nCells,2);
handedness = NaN;
rotationRad = NaN;
fitR = NaN;

if nCells < 3 || size(traces,1) < 3
    return;
end

XtimeByCell = fillmissing(traces, 'linear', 1, 'EndValues', 'nearest');
XtimeByCell(~isfinite(XtimeByCell)) = 0;

if doZscore
    muCell = mean(XtimeByCell, 1);
    sdCell = std(XtimeByCell, 0, 1);
    sdCell(~isfinite(sdCell) | sdCell <= 0) = 1;
    XtimeByCell = bsxfun(@rdivide, bsxfun(@minus, XtimeByCell, muCell), sdCell);
end

% Neurons are observations; time samples are features, matching the poster
% pipeline. The Gram matrix gives the same cell scores as an economy SVD but
% avoids constructing thousands of temporal PC coefficients.
X = XtimeByCell';
X = bsxfun(@minus, X, mean(X,1));
G = X * X';
G = (G + G') ./ 2;
[V,D] = eig(G);
[lambda, order] = sort(real(diag(D)), 'descend');
if numel(lambda) < 2 || lambda(2) <= 0
    return;
end
score = bsxfun(@times, real(V(:,order(1:2))), sqrt(max(lambda(1:2),0))');

% Kasa circle centre/scale, as used in the upstream rPC script.
[center, radius] = fitCircleKasa(score(:,1), score(:,2));
if ~isfinite(radius) || radius <= 0
    center = mean(score,1,'omitnan');
    radius = median(sqrt(sum((score-center).^2,2)), 'omitnan');
end
rawXY = (score - center) ./ radius;

z = complex(rawXY(:,1), rawXY(:,2));
good = isfinite(preferredPhase) & isfinite(real(z)) & isfinite(imag(z));
w = double(weights(:));
w(~isfinite(w) | w < 0) = 0;
w = w + 0.05; % keep broad cells but favour sharply tuned cells
good = good & isfinite(w);
if sum(good) < 3
    alignedXY = rawXY;
    return;
end

bestR = -inf;
bestSign = 1;
bestDelta = 0;
for s = [1 -1]
    if s == 1
        zTest = z;
    else
        zTest = conj(z);
    end
    residual = preferredPhase(good) - angle(zTest(good));
    c = sum(w(good) .* exp(1i .* residual));
    thisR = abs(c) ./ sum(w(good));
    if thisR > bestR
        bestR = thisR;
        bestSign = s;
        bestDelta = angle(c);
    end
end

if bestSign == -1
    z = conj(z);
end
z = z .* exp(1i .* bestDelta);

alignedXY = [real(z), imag(z)];
handedness = bestSign;
rotationRad = bestDelta;
fitR = bestR;
end

function [center, radius] = fitCircleKasa(x, y)
x = double(x(:));
y = double(y(:));
good = isfinite(x) & isfinite(y);
x = x(good);
y = y(good);
if numel(x) < 3
    center = [NaN NaN];
    radius = NaN;
    return;
end
A = [2*x, 2*y, ones(size(x))];
b = x.^2 + y.^2;
p = A \ b;
center = p(1:2)';
radius = sqrt(max(p(3) + sum(center.^2), 0));
end

%% ========================================================================
function T = makeCellTable(F)
n = F.nCells;
T = table( ...
    repmat(F.morph,n,1), repmat(F.fish,n,1), repmat(F.file,n,1), ...
    repmat(F.phaseCurveSource,n,1), repmat(F.phaseValiditySource,n,1), ...
    F.cellIds(:), ...
    rad2deg(F.preferredPhaseOriginalRad(:)), ...
    rad2deg(F.preferredPhaseAlignedRad(:)), F.FWHMDeg(:), F.MVL(:), ...
    F.InformationBits(:), F.PeakToBaselineSNR(:), F.PeakPhaseRate(:), ...
    F.PhaseSelectivity(:), ...
    F.PhaseTuningInformationBits(:), F.PhaseTuningPValue(:), ...
    F.PhaseTuningQValue(:), F.PhaseTuned(:), F.StoredPhaseTuned(:), ...
    F.StoredFWHMDeg(:), F.StoredMVL(:), F.StoredInformationBits(:), ...
    F.StoredPeakToBaselineSNR(:), F.StoredPhaseSelectivity(:), ...
    F.InformationFromStoredORI(:), ...
    F.SplitHalfReliability(:), ...
    'VariableNames', {'Morph','Fish','SourceFile','PhaseCurveSource', ...
    'PhaseValiditySource', ...
    'CandidateCellID', ...
    'PreferredPhaseOriginalDeg','PreferredPhaseAlignedDeg','FWHMDeg','MVL', ...
    'InformationBits','PeakToBaselineSNR','PeakPhaseRate','PhaseSelectivity', ...
    'PhaseTuningInformationBits','PhaseTuningPValue','PhaseTuningQValue', ...
    'PhaseTuned','StoredPhaseTuned','StoredFWHMDeg','StoredMVL', ...
    'StoredInformationBits','StoredPeakToBaselineSNR','StoredPhaseSelectivity', ...
    'InformationFromStoredORI', ...
    'SplitHalfReliability'});
end

function T = makeAlignmentTable(F)
T = table(F.morph, F.fish, F.file, F.phaseCurveSource, ...
    F.phaseValiditySource, ...
    F.nCells, nnz(F.PhaseTuned), ...
    nnz(F.PhaseTuned)./max(F.nCells,1), F.PhaseTuningNValidFrames, F.phaseSign, ...
    rad2deg(F.phaseOffsetRad), F.pcaHandedness, rad2deg(F.pcaRotationRad), ...
    F.pcaPhaseFitR, ...
    'VariableNames', {'Morph','Fish','SourceFile','PhaseCurveSource', ...
    'PhaseValiditySource', ...
    'NCandidateCells', ...
    'NPhaseTuned','PhaseTunedFraction','PhaseTuningNValidFrames','PhaseSign', ...
    'PhaseOffsetDeg','PCAHandedness','PCARotationDeg','PCAPhaseFitR'});
end

function audit = buildORIFeatureMethodAudit(cellTable)
% Compare direct ORI-method calculations with the canonical values exported
% in candidate_features. This table is diagnostic and never drives tests.
definitions = table( ...
    string({'CircularFWHM';'MeanVectorLength';'DirectionalInformation'; ...
    'PeakToBaselineSNR';'PhaseSelectivity'}), ...
    string({'FWHMDeg';'MVL';'PhaseTuningInformationBits'; ...
    'PeakToBaselineSNR';'PhaseSelectivity'}), ...
    string({'StoredFWHMDeg';'StoredMVL';'StoredInformationBits'; ...
    'StoredPeakToBaselineSNR';'StoredPhaseSelectivity'}), ...
    'VariableNames', {'Metric','ComputedVariable','StoredVariable'});

fishKeys = unique(cellTable.Morph + "|" + cellTable.Fish, 'stable');
rows = cell(numel(fishKeys).*height(definitions),1);
counter = 0;
for f = 1:numel(fishKeys)
    useFish = (cellTable.Morph + "|" + cellTable.Fish) == fishKeys(f);
    morph = cellTable.Morph(find(useFish,1));
    fish = cellTable.Fish(find(useFish,1));
    for d = 1:height(definitions)
        x = double(cellTable.(char(definitions.ComputedVariable(d)))(useFish));
        y = double(cellTable.(char(definitions.StoredVariable(d)))(useFish));
        good = isfinite(x) & isfinite(y);
        difference = x(good)-y(good);
        if isempty(difference)
            medianAbsoluteDifference = NaN;
            maximumAbsoluteDifference = NaN;
            medianSignedDifference = NaN;
            correlation = NaN;
        else
            medianAbsoluteDifference = median(abs(difference));
            maximumAbsoluteDifference = max(abs(difference));
            medianSignedDifference = median(difference);
            if numel(difference) >= 3 && std(x(good)) > 0 && std(y(good)) > 0
                correlationMatrix = corrcoef(x(good),y(good));
                correlation = correlationMatrix(1,2);
            else
                correlation = NaN;
            end
        end
        counter = counter+1;
        rows{counter} = table(morph,fish,definitions.Metric(d),nnz(good), ...
            medianAbsoluteDifference,maximumAbsoluteDifference, ...
            medianSignedDifference,correlation, ...
            'VariableNames', {'Morph','Fish','Metric','NCompared', ...
            'MedianAbsoluteDifference','MaximumAbsoluteDifference', ...
            'MedianSignedDifference','Correlation'});
    end
end
audit = vertcat(rows{1:counter});
end

function [metricCells, description] = selectMetricCellSubset(cellTable, cfg)
switch cfg.MetricCellSubset
    case 'all_candidates'
        keep = true(height(cellTable),1);
        description = 'all candidate cells (not the poster phase-tuned subset)';

    case 'sel_phi_threshold'
        if ~any(isfinite(cellTable.PhaseSelectivity))
            error(['MetricCellSubset is ''sel_phi_threshold'', but ORI Step 19 ' ...
                'phase selectivity could not be computed.']);
        end
        keep = isfinite(cellTable.PhaseSelectivity) & ...
            cellTable.PhaseSelectivity >= cfg.PhaseSelectivityThreshold;
        description = sprintf('phase-selective cells: sel_phi >= %.6g', ...
            cfg.PhaseSelectivityThreshold);

    case 'phase_tuned'
        keep = logical(cellTable.PhaseTuned) & ...
            isfinite(cellTable.PhaseTuningQValue) & ...
            cellTable.PhaseTuningQValue < cfg.PhaseTuningFDRAlpha;
        description = sprintf(['ORI_V15 phase-tuned cells: circular-shift ' ...
            'information test, within-fish BH-FDR q < %.3g'], ...
            cfg.PhaseTuningFDRAlpha);

    otherwise
        error('Unsupported MetricCellSubset.');
end

metricCells = cellTable(keep,:);
if isempty(metricCells)
    error('The selected Figure 4E-H cell subset is empty. Check the tuning criterion.');
end

allFish = unique(cellTable.Morph + "|" + cellTable.Fish);
keptFish = unique(metricCells.Morph + "|" + metricCells.Fish);
missingFish = setdiff(allFish,keptFish);
if ~isempty(missingFish)
    warning('%d fish have no cells in the selected metric subset: %s', ...
        numel(missingFish), strjoin(cellstr(missingFish), ', '));
end
end

function defs = posterMetricDefinitions()
defs = table( ...
    string({'FWHMDeg';'MVL';'InformationBits';'PeakToBaselineSNR'}), ...
    string({'E Angular resolution';'F Consistency';'G Information';'H Reliability / SNR'}), ...
    string({'FWHM (degrees; lower = sharper)';'Mean vector length'; ...
    'Directional information (bits)';'Peak / baseline (lowest 20% of bins)'}), ...
    [false; false; false; true], ...
    'VariableNames', {'Variable','Title','YLabel','LogScale'});
end

function defs = trueReliabilityDefinition()
defs = table("SplitHalfReliability", "True split-half reliability", ...
    "Alternating-block tuning correlation (r)", false, ...
    'VariableNames', {'Variable','Title','YLabel','LogScale'});
end

function summary = buildFishSummary(cellTable, metricDefs, nBootstrap)
rows = {};
fishKeys = unique(cellTable.Morph + "|" + cellTable.Fish, 'stable');
counter = 0;

allDefs = [metricDefs; trueReliabilityDefinition()];
for i = 1:numel(fishKeys)
    idxFish = (cellTable.Morph + "|" + cellTable.Fish) == fishKeys(i);
    morph = cellTable.Morph(find(idxFish,1));
    fish = cellTable.Fish(find(idxFish,1));

    for m = 1:height(allDefs)
        variable = char(allDefs.Variable(m));
        x = double(cellTable.(variable)(idxFish));
        x = x(isfinite(x));
        if isempty(x)
            value = NaN;
            ci = [NaN NaN];
        else
            value = median(x);
            ci = bootstrapMedianCI(x, nBootstrap);
        end

        counter = counter + 1;
        rows{counter,1} = table(morph, fish, allDefs.Variable(m), numel(x), ...
            value, ci(1), ci(2), ...
            'VariableNames', {'Morph','Fish','Metric','NCells','Value','CILow','CIHigh'}); %#ok<AGROW>
    end
end
summary = vertcat(rows{:});
end

function ci = bootstrapMedianCI(x, nBootstrap)
x = x(:);
if numel(x) == 1
    ci = [x x];
    return;
end
idx = randi(numel(x), numel(x), nBootstrap);
boot = median(x(idx), 1);
ci = percentileLocal(boot, [2.5 97.5]);
end

%% ========================================================================
function statsTable = runMorphStatistics(fishSummary, metricDefs, cfg)
% Report each a-priori Surface-Pachon and Surface-Molino contrast directly.
% No across-contrast multiplicity adjustment is applied here.
rows = {};
counter = 0;

for m = 1:height(metricDefs)
    metric = metricDefs.Variable(m);
    use = fishSummary.Metric == metric & isfinite(fishSummary.Value) & ...
        fishSummary.NCells >= cfg.MinimumMetricCellsPerFish;
    values = double(fishSummary.Value(use));
    groups = fishSummary.Morph(use);

    % These are the only confirmatory comparisons. With two groups,
    % Kruskal-Wallis is equivalent to a Mann-Whitney/Wilcoxon rank-sum test.
    % We compute the KW rank statistic and obtain p by permuting FISH labels.
    for c = 1:size(cfg.PlannedContrasts,1)
        g1 = string(cfg.PlannedContrasts{c,1});
        g2 = string(cfg.PlannedContrasts{c,2});
        x1 = values(strcmpi(groups,g1));
        x2 = values(strcmpi(groups,g2));
        [pPair, hStatistic] = permutationTwoGroupKW(x1, x2, cfg.NPermutations);
        medianDifference = median(x2,'omitnan') - median(x1,'omitnan');
        cliffsDelta = cliffsDeltaLocal(x2, x1); % positive: group2 > group1

        counter = counter + 1;
        rows{counter,1} = table(metric, "planned_two_group_KW", g1, g2, ...
            numel(x1), numel(x2), hStatistic, medianDifference, cliffsDelta, ...
            pPair, ...
            'VariableNames', {'Metric','Test','Group1','Group2','NFish1', ...
            'NFish2','KWStatistic','MedianDifference','CliffsDelta', ...
            'PValue'}); %#ok<AGROW>
    end
end

statsTable = vertcat(rows{:});
end

function statsTable = runLegacyPooledCellDiagnostic(cellTable, metricDefs, cfg)
% Diagnostic only. This deliberately ignores fish nesting and therefore
% must not be reported as the confirmatory morph test. It is returned in
% memory so an old neuron-pooled poster analysis can be identified explicitly.
rows = {};
counter = 0;

for m = 1:height(metricDefs)
    metric = metricDefs.Variable(m);
    variable = char(metric);
    values = double(cellTable.(variable));
    groups = cellTable.Morph;

    for c = 1:size(cfg.PlannedContrasts,1)
        g1 = string(cfg.PlannedContrasts{c,1});
        g2 = string(cfg.PlannedContrasts{c,2});
        x1 = values(strcmpi(groups,g1) & isfinite(values));
        x2 = values(strcmpi(groups,g2) & isfinite(values));

        if isempty(x1) || isempty(x2)
            H = NaN;
            pValue = NaN;
            delta = NaN;
            medianDifference = NaN;
        else
            allValues = [x1(:); x2(:)];
            n1 = numel(x1);
            n2 = numel(x2);
            ranks = tiedRanksLocal(allValues);
            H = twoGroupKWFromRanks(ranks,n1);
            pValue = gammainc(H./2, 0.5, 'upper');
            u2 = sum(ranks(n1+1:end)) - n2.*(n2+1)./2;
            delta = 2.*u2./(n1.*n2)-1; % positive: group2 > group1
            medianDifference = median(x2)-median(x1);
        end

        counter = counter+1;
        rows{counter,1} = table(metric, ...
            "DIAGNOSTIC_ONLY_pooled_neurons_not_independent", g1, g2, ...
            numel(x1), numel(x2), H, medianDifference, delta, pValue, ...
            'VariableNames', {'Metric','Test','Group1','Group2','NCells1', ...
            'NCells2','KWStatistic','MedianDifference','CliffsDelta', ...
            'PValue'}); %#ok<AGROW>
    end
end

statsTable = vertcat(rows{:});
end

function [p, observedH] = permutationTwoGroupKW(x1, x2, nPerm)
x1 = x1(isfinite(x1));
x2 = x2(isfinite(x2));
if isempty(x1) || isempty(x2)
    p = NaN;
    observedH = NaN;
    return;
end
allValues = [x1(:); x2(:)];
n1 = numel(x1);
ranks = tiedRanksLocal(allValues);
observedH = twoGroupKWFromRanks(ranks, n1);
count = 0;
for i = 1:nPerm
    idx = randperm(numel(allValues));
    thisH = twoGroupKWFromRanks(ranks(idx), n1);
    count = count + (thisH >= observedH - eps);
end
p = (count + 1) ./ (nPerm + 1);
end

function H = twoGroupKWFromRanks(ranks, n1)
N = numel(ranks);
n2 = N-n1;
mean1 = mean(ranks(1:n1));
mean2 = mean(ranks(n1+1:end));
H = 12 ./ (N.*(N+1)) .* ...
    (n1.*(mean1-(N+1)./2).^2 + n2.*(mean2-(N+1)./2).^2);

% Tie correction. The tie structure is invariant under permutation.
[~,~,ic] = unique(ranks);
tieCounts = accumarray(ic,1);
tieCorrection = 1 - sum(tieCounts.^3-tieCounts) ./ (N.^3-N);
if tieCorrection > 0
    H = H ./ tieCorrection;
end
end

function ranks = tiedRanksLocal(x)
x = x(:);
[sorted, order] = sort(x);
ranksSorted = nan(size(sorted));
i = 1;
while i <= numel(sorted)
    j = i;
    while j < numel(sorted) && sorted(j+1) == sorted(i)
        j = j+1;
    end
    ranksSorted(i:j) = mean(i:j);
    i = j+1;
end
ranks = nan(size(x));
ranks(order) = ranksSorted;
end

function delta = cliffsDeltaLocal(x, y)
x = x(isfinite(x));
y = y(isfinite(y));
if isempty(x) || isempty(y)
    delta = NaN;
    return;
end
greater = 0;
less = 0;
for i = 1:numel(x)
    greater = greater + sum(x(i) > y);
    less = less + sum(x(i) < y);
end
delta = (greater-less) ./ (numel(x).*numel(y));
end

%% ========================================================================
function pathOut = plotAlignedPCA(fishData, cfg, morphColors)
fig = figure('Color','w','Visible',cfg.Visible,'Name','Figure 3A aligned PCA');
tl = tiledlayout(fig, 1, numel(cfg.MorphOrder), 'TileSpacing','compact','Padding','compact');
title(tl, {'Figure 3A - aligned candidate-cell PCA projection', ...
    'PCA ring normalized within fish; colour = aligned preferred network phase'});

for m = 1:numel(cfg.MorphOrder)
    ax = nexttile(tl);
    hold(ax,'on');
    nFish = 0;
    nCells = 0;
    for i = 1:numel(fishData)
        if ~strcmpi(fishData(i).morph, cfg.MorphOrder{m})
            continue;
        end
        xy = fishData(i).alignedPcaXY;
        c = rad2deg(fishData(i).preferredPhaseAlignedRad);
        good = all(isfinite(xy),2) & isfinite(c);
        if any(good)
            try
                scatter(ax, xy(good,1), xy(good,2), 13, c(good), 'filled', ...
                    'MarkerFaceAlpha',0.40, 'MarkerEdgeAlpha',0.15);
            catch
                scatter(ax, xy(good,1), xy(good,2), 13, c(good), 'filled');
            end
        end
        nFish = nFish + 1;
        nCells = nCells + sum(good);
    end
    plot(ax, cos(linspace(0,2*pi,361)), sin(linspace(0,2*pi,361)), ':', ...
        'Color',[0.25 0.25 0.25], 'LineWidth',1);
    xline(ax,0,':','Color',[0.65 0.65 0.65]);
    yline(ax,0,':','Color',[0.65 0.65 0.65]);
    axis(ax,'equal');
    xlim(ax,[-2.2 2.2]);
    ylim(ax,[-2.2 2.2]);
    xlabel(ax,'aligned PC1');
    ylabel(ax,'aligned PC2');
    title(ax, sprintf('%s (n=%d fish, %d cells)', cfg.MorphOrder{m}, nFish, nCells), ...
        'Color',morphColors(m,:));
    box(ax,'off');
    caxis(ax,[-180 180]);
end
colormap(fig, hsv(256));
cb = colorbar;
try
    cb.Layout.Tile = 'east';
catch
end
cb.Label.String = 'Aligned preferred phase (deg)';

pathOut = saveFigureBoth(fig, cfg.OutputDir, 'Figure3A_aligned_PC_projection');
end

function pathOut = plotTuningHeatmaps(fishData, centersDeg, cfg)
fig = figure('Color','w','Visible',cfg.Visible,'Name','Individual cell HD tuning');
tl = tiledlayout(fig, 1, numel(cfg.MorphOrder), 'TileSpacing','compact','Padding','compact');
title(tl, {'Individual candidate-cell HD tuning', ...
    'Activity as a function of aligned network phase; each row min-max normalized'});

for m = 1:numel(cfg.MorphOrder)
    curves = [];
    pref = [];
    nFish = 0;
    for i = 1:numel(fishData)
        if strcmpi(fishData(i).morph, cfg.MorphOrder{m})
            curves = [curves; fishData(i).tuningCurvesAligned]; %#ok<AGROW>
            pref = [pref; fishData(i).preferredPhaseAlignedRad]; %#ok<AGROW>
            nFish = nFish + 1;
        end
    end
    [~, order] = sort(pref);
    curves = curves(order,:);
    curves = normalizeRows01(curves);

    ax = nexttile(tl);
    imagesc(ax, centersDeg, 1:size(curves,1), curves);
    axis(ax,'xy');
    xline(ax,0,'w:','LineWidth',1);
    xlabel(ax,'Aligned network phase (deg)');
    ylabel(ax,'Candidate cells sorted by preferred phase');
    title(ax,sprintf('%s (n=%d fish, %d cells)',cfg.MorphOrder{m},nFish,size(curves,1)));
    caxis(ax,[0 1]);
end
colormap(fig, parula(256));
cb = colorbar;
try
    cb.Layout.Tile = 'east';
catch
end
cb.Label.String = 'Within-cell normalized activity';

pathOut = saveFigureBoth(fig, cfg.OutputDir, 'Individual_cell_HD_tuning_heatmaps');
end

function pathOut = plotMetricPanels(cellTable, fishSummary, statsTable, defs, cfg, morphColors, mode)
fig = figure('Color','w','Visible',cfg.Visible,'Name','Candidate cell metrics');
tl = tiledlayout(fig, 1, height(defs), 'TileSpacing','compact','Padding','compact');

if strcmp(mode,'poster')
    title(tl, {'Figure 4E-H - candidate-cell HD tuning properties', ...
        sprintf(['Bars: mean of fish medians +/- SEM; large points: fish medians; ' ...
        'brackets: planned fish-level KW tests; subset = %s'], ...
        strrep(cfg.MetricCellSubset,'_',' '))});
else
    title(tl, {'True tuning reliability (additional analysis)', ...
        'Bars: mean of fish medians +/- SEM; alternating 60-s blocks; brackets: planned fish-level KW tests'});
end

for d = 1:height(defs)
    ax = nexttile(tl);
    hold(ax,'on');
    variable = char(defs.Variable(d));

    for m = 1:numel(cfg.MorphOrder)
        morph = string(cfg.MorphOrder{m});

        useFish = fishSummary.Metric == defs.Variable(d) & strcmpi(fishSummary.Morph,morph) ...
            & isfinite(fishSummary.Value);
        yFish = double(fishSummary.Value(useFish));
        drawFishBar(ax, m, yFish, morphColors(m,:), defs.LogScale(d));

        useCell = strcmpi(cellTable.Morph,morph) & isfinite(cellTable.(variable));
        yCell = double(cellTable.(variable)(useCell));
        jitter = (rand(size(yCell))-0.5).*0.46;
        try
            scatter(ax, m+jitter, yCell, 9, morphColors(m,:), 'filled', ...
                'MarkerFaceAlpha',0.12, 'MarkerEdgeAlpha',0.05);
        catch
            scatter(ax, m+jitter, yCell, 9, morphColors(m,:), 'filled');
        end

        if ~isempty(yFish)
            xFish = m + linspace(-0.12,0.12,numel(yFish))';
            scatter(ax, xFish, yFish, 42, morphColors(m,:), 'filled', ...
                'MarkerEdgeColor','k','LineWidth',0.75);
        end
    end

    xlim(ax,[0.45 numel(cfg.MorphOrder)+0.55]);
    set(ax,'XTick',1:numel(cfg.MorphOrder),'XTickLabel',cfg.MorphOrder);
    ylabel(ax,defs.YLabel(d));
    box(ax,'off');
    grid(ax,'on');
    ax.GridAlpha = 0.12;
    if defs.LogScale(d)
        set(ax,'YScale','log');
    end

    addPlannedMorphBars(ax, statsTable, defs.Variable(d), ...
        cfg.MorphOrder, defs.LogScale(d));
    title(ax, char(defs.Title(d)));
end

if strcmp(mode,'poster')
    name = 'Figure4E_H_candidate_cell_metrics';
else
    name = 'True_split_half_reliability';
end
pathOut = saveFigureBoth(fig, cfg.OutputDir, name);
end

function addPlannedMorphBars(ax, statsTable, metric, morphOrder, logScale)
rows = find(statsTable.Metric == metric & ...
    statsTable.Test == "planned_two_group_KW" & isfinite(statsTable.PValue));
if isempty(rows)
    return;
end

% Draw the shorter comparison first and stack the wider comparison above it.
spans = nan(size(rows));
for i = 1:numel(rows)
    x1 = find(strcmpi(string(morphOrder), statsTable.Group1(rows(i))), 1);
    x2 = find(strcmpi(string(morphOrder), statsTable.Group2(rows(i))), 1);
    if ~isempty(x1) && ~isempty(x2)
        spans(i) = abs(x2-x1);
    end
end
[~, order] = sort(spans, 'ascend', 'MissingPlacement', 'last');
rows = rows(order);

limits = ylim(ax);
if logScale
    if limits(1) <= 0 || limits(2) <= 0
        return;
    end
    factor = 1.32;
    base = limits(2);
    for i = 1:numel(rows)
        r = rows(i);
        x1 = find(strcmpi(string(morphOrder), statsTable.Group1(r)), 1);
        x2 = find(strcmpi(string(morphOrder), statsTable.Group2(r)), 1);
        if isempty(x1) || isempty(x2)
            continue;
        end
        y = base .* factor.^i;
        tickBottom = y ./ 1.05;
        plot(ax, [x1 x1 x2 x2], [tickBottom y y tickBottom], 'k-', ...
            'LineWidth', 1.1, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(ax, mean([x1 x2]), y .* 1.025, significanceLabel(statsTable.PValue(r)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', 'Clipping', 'off');
    end
    ylim(ax, [limits(1), base .* factor.^(numel(rows) + 0.55)]);
else
    span = limits(2) - limits(1);
    if ~isfinite(span) || span <= 0
        span = max(abs(limits));
        if ~isfinite(span) || span <= 0
            span = 1;
        end
    end
    base = limits(2) + 0.08 .* span;
    step = 0.13 .* span;
    tick = 0.035 .* span;
    for i = 1:numel(rows)
        r = rows(i);
        x1 = find(strcmpi(string(morphOrder), statsTable.Group1(r)), 1);
        x2 = find(strcmpi(string(morphOrder), statsTable.Group2(r)), 1);
        if isempty(x1) || isempty(x2)
            continue;
        end
        y = base + (i-1) .* step;
        plot(ax, [x1 x1 x2 x2], [y-tick y y y-tick], 'k-', ...
            'LineWidth', 1.1, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(ax, mean([x1 x2]), y + 0.015 .* span, ...
            significanceLabel(statsTable.PValue(r)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', 'Clipping', 'off');
    end
    ylim(ax, [limits(1), base + (numel(rows)-1) .* step + 0.14 .* span]);
end
addPlannedPValueText(ax, statsTable(rows,:));
end

function addPlannedPValueText(ax, rows)
pieces = strings(height(rows),1);
for i = 1:height(rows)
    pieces(i) = string(rows.Group1(i)) + "-" + string(rows.Group2(i)) + ...
        ": p=" + string(sprintf('%.3g',rows.PValue(i)));
end
text(ax,0.02,0.98,strjoin(cellstr(pieces),newline), ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','top','FontSize',8,'Interpreter','none', ...
    'BackgroundColor','w','Margin',2,'Clipping','off');
end

function label = significanceLabel(p)
if p < 0.001
    label = '***';
elseif p < 0.01
    label = '**';
elseif p < 0.05
    label = '*';
else
    label = 'ns';
end
end

function drawFishBar(ax, x, y, color, logScale)
y = y(isfinite(y));
if isempty(y)
    return;
end
barMean = mean(y,'omitnan');
if numel(y) > 1
    sem = std(y,0,'omitnan') ./ sqrt(numel(y));
else
    sem = 0;
end

b = bar(ax, x, barMean, 0.62, 'FaceColor',color, ...
    'EdgeColor',color,'LineWidth',1.3);
if logScale
    % The SNR panel is logarithmic; its natural neutral value is one.
    % A positive base keeps the bar drawable after the axis becomes log.
    b.BaseValue = 1;
end
try
    b.FaceAlpha = 0.25;
catch
end

lowerError = sem;
if logScale
    positiveFloor = max(realmin('double'), barMean .* 1e-6);
    lowerError = min(sem, max(barMean-positiveFloor,0));
end
errorbar(ax, x, barMean, lowerError, sem, 'LineStyle','none', ...
    'Color',[0.15 0.15 0.15], 'LineWidth',1.3, 'CapSize',8);
end

function pathOut = plotPhaseAlignmentQC(fishData, centersDeg, cfg, morphColors)
fig = figure('Color','w','Visible',cfg.Visible,'Name','Phase alignment QC');
ax = axes(fig);
hold(ax,'on');

lineHandles = gobjects(numel(cfg.MorphOrder),1);
for m = 1:numel(cfg.MorphOrder)
    morphCurves = [];
    for i = 1:numel(fishData)
        if ~strcmpi(fishData(i).morph,cfg.MorphOrder{m})
            continue;
        end
        c = mean(fishData(i).tuningCurvesAligned,1,'omitnan');
        c = normalizeRows01(c);
        paleColor = 0.78.*[1 1 1] + 0.22.*morphColors(m,:);
        plot(ax,centersDeg,c,'Color',paleColor,'LineWidth',0.8);
        morphCurves = [morphCurves; c]; %#ok<AGROW>
    end
    if ~isempty(morphCurves)
        lineHandles(m) = plot(ax,centersDeg,mean(morphCurves,1,'omitnan'), ...
            'Color',morphColors(m,:),'LineWidth',2.5);
    end
end
xline(ax,cfg.PhaseReferenceDeg,'k:','Reference','LineWidth',1.2);
xlabel(ax,'Aligned network phase (deg)');
ylabel(ax,'Normalized fish-average tuning');
title(ax,{'Phase-alignment quality control', ...
    'Thin lines are individual fish; thick lines are morph means'});
legend(ax,lineHandles,cfg.MorphOrder,'Location','best');
xlim(ax,[-180 180]);
box(ax,'off');
grid(ax,'on');
ax.GridAlpha = 0.12;

pathOut = saveFigureBoth(fig, cfg.OutputDir, 'Phase_alignment_QC');
end

function pathOut = saveFigureBoth(fig, outputDir, baseName)
pathOut = fullfile(outputDir, [baseName '.png']);
try
    exportgraphics(fig, pathOut, 'Resolution', 300);
catch
    print(fig, pathOut, '-dpng', '-r300');
end
svgPath = fullfile(outputDir, [baseName '.svg']);
try
    exportgraphics(fig, svgPath, 'ContentType', 'vector');
catch
    print(fig, svgPath, '-dsvg');
end
end

%% ========================================================================
function out = normalizeRows01(x)
rowMin = min(x,[],2,'omitnan');
rowMax = max(x,[],2,'omitnan');
range = rowMax-rowMin;
range(~isfinite(range) | range <= 0) = 1;
out = bsxfun(@rdivide, bsxfun(@minus,x,rowMin), range);
end

function value = lookupRecordingValue(container, fish, defaultValue)
value = defaultValue;
if isempty(container)
    return;
end
key = matlab.lang.makeValidName(char(fish));
if isstruct(container)
    if isfield(container,key)
        value = double(container.(key));
    end
elseif isa(container,'containers.Map')
    if isKey(container,char(fish))
        value = double(container(char(fish)));
    elseif isKey(container,key)
        value = double(container(key));
    end
end
if ~isscalar(value) || ~isfinite(value)
    error('Manual phase value for %s must be one finite scalar.', fish);
end
end

function q = percentileLocal(x, percent)
x = sort(x(isfinite(x)));
percent = percent(:)';
q = nan(size(percent));
if isempty(x)
    return;
elseif numel(x) == 1
    q(:) = x;
    return;
end
for i = 1:numel(percent)
    pos = 1 + (numel(x)-1) .* percent(i) ./ 100;
    lo = floor(pos);
    hi = ceil(pos);
    frac = pos-lo;
    q(i) = x(lo).*(1-frac) + x(hi).*frac;
end
end

function angle = wrapToPiLocal(angle)
angle = mod(angle + pi, 2*pi) - pi;
end
