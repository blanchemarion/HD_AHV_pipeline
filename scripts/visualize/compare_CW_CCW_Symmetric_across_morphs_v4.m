function results = compare_CW_CCW_Symmetric_across_morphs_v4(varargin)
% COMPARE_CW_CCW_SYMMETRIC_ACROSS_MORPHS
% Class-focused cross-fish analysis using ONLY files named
% *_candidate_neurons_clean_classified.mat. The clean, unclassified sibling
% file is never searched for or loaded.
%
% Main figures (rows are morphs whenever possible):
%   1) Preferred firing direction: 3 x 4 binned-probability panel.
%   2) Cell projection in the common candidate-cell PC space: 3 x 3.
%   3) |pairwise activity correlation| versus circular phase difference:
%      3 x 6 for CW-CW, CCW-CCW, Sym-Sym, CW-CCW, CW-Sym and CCW-Sym.
%   4) Class-specific rPC phase versus estimated heading:
%      fish-level |Pearson r| summaries plus representative z-score overlays.
%   5) Figure 4E-H properties in a class x metric 3 x 4 panel, with
%      Surface-Molino and Surface-Pachon planned-test brackets.
%
% Every retained figure is produced twice: once from its original input
% subset and once from the highly tuned cells within each fish.
%
% Fish is the independent unit. Direction histograms are normalized within
% fish before morph averaging, and pairwise-correlation curves are binned
% within fish before morph averaging. A fish with many neurons or pairs
% therefore does not dominate the displayed morph mean.
%
% BASIC USAGE
% -----------
%   R = compare_CW_CCW_Symmetric_across_morphs_v4;
%
%   R = compare_CW_CCW_Symmetric_across_morphs_v4( ...
%       'DataRoot','/path/to/data', ...
%       'OutputDir','/path/to/data/class_comparison_functional');
%
% Explicit inputs must also be classified files:
%   R = compare_CW_CCW_Symmetric_across_morphs_v4( ...
%       'Files',{'/path/recXX_..._candidate_neurons_clean_classified.mat'});
%
% PHASE CONVENTION
% ----------------
% By default the upstream network_phase_rad convention is preserved because
% ORI anatomy-aligns the rPC space. This is essential if absolute preferred
% directions are to be compared across fish. If files do not share a common
% upstream convention, use 'PhaseAlignment','population_peak' to rotate each
% fish's all-candidate population-tuning peak to 0 deg; interpret the result
% as relative phase only.
%
% CLASS-SPECIFIC rPC PHASE
% ------------------------
% All candidate cells define one aligned PC ring per fish. For each class,
% the rPC population vector is then recomputed using only that class's
% traces and its positions on this common ring. Thus class trajectories are
% comparable and do not acquire separate arbitrary PCA rotations.
%
% HEADING INPUT
% -------------
% Heading is first sought in the classified file (including ORI heading
% overlay exports). If absent, the newest matching *_swimResults_pass2.mat
% in the recording folder is used automatically. Missing heading causes
% only the heading panels for that fish to be skipped.
%
% MATLAB R2020b or newer is recommended. No Circular Statistics Toolbox is
% required.

%% Parse inputs
defaultCfg = pipeline_config();
p = inputParser;
p.FunctionName = mfilename;

addParameter(p, 'DataRoot', defaultCfg.ClassifiedCandidateDir, @(x) ischar(x) || isstring(x));
addParameter(p, 'Files', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'OutputDir', defaultCfg.ClassComparisonDir, @(x) ischar(x) || isstring(x));
addParameter(p, 'MorphOrder', {'Surface','Molino','Pachon'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'ClassOrder', {'CW','CCW','Symmetric'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'IncludeFish', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'ExcludeFish', {}, @(x) ischar(x) || isstring(x) || iscell(x));

addParameter(p, 'PhaseAlignment', 'preserve', @(x) ischar(x) || isstring(x));
addParameter(p, 'PhaseReferenceDeg', 0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PhaseSigns', struct(), @(x) isstruct(x) || isa(x, 'containers.Map') || isempty(x));
addParameter(p, 'PhaseOffsetsDeg', struct(), @(x) isstruct(x) || isa(x, 'containers.Map') || isempty(x));

addParameter(p, 'NPhaseBins', 36, @(x) isnumeric(x) && isscalar(x) && x >= 12);
addParameter(p, 'PreferredDirectionNBins', 18, @(x) isnumeric(x) && isscalar(x) && x >= 6);
addParameter(p, 'PreferredPhaseKDEBandwidthDeg', 20, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 180);
addParameter(p, 'PreferredPhaseKDEGridPoints', 181, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 73);
addParameter(p, 'MinimumPreferredPhaseCellsPerFish', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'PairPhaseDifferenceNBins', 12, @(x) isnumeric(x) && isscalar(x) && x >= 6);
addParameter(p, 'ZScoreTracesForPCA', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ReliabilityBlockSec', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MinReliabilityBins', 12, @(x) isnumeric(x) && isscalar(x) && x >= 3);

addParameter(p, 'BehaviorFilePattern', '*_swimResults_pass2.mat', @(x) ischar(x) || isstring(x));
addParameter(p, 'HeadingUnits', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'HeadingOverlayWindowSec', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PhaseConfidencePercentile', 10, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 100);
addParameter(p, 'MinimumClassCellsForPhase', 3, @(x) isnumeric(x) && isscalar(x) && x >= 2);

addParameter(p, 'MetricCellSubset', 'phase_tuned', @(x) ischar(x) || isstring(x));
addParameter(p, 'PhaseSelectivityThreshold', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PhaseTuningFDRAlpha', 0.01, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'PhaseTuningNShuffles', 500, @(x) isnumeric(x) && isscalar(x) && x >= 20);
addParameter(p, 'PhaseTuningSmoothBins', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'PhaseTuningMinShiftFraction', 0.10, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'PhaseTuningMaxShiftFraction', 0.90, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'MinimumMetricCellsPerFish', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'PlannedContrasts', {'Surface','Molino';'Surface','Pachon'}, ...
    @(x) iscell(x) && size(x,2) == 2);

addParameter(p, 'NBootstrap', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
addParameter(p, 'NPermutations', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
addParameter(p, 'RandomSeed', 1701, @(x) isnumeric(x) && isscalar(x));
% Deprecated compatibility option; output is always PNG and SVG.
addParameter(p, 'SaveFigFiles', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Visible', defaultCfg.FigureVisible, @(x) ischar(x) || isstring(x));

parse(p, varargin{:});
cfg = p.Results;
cfg.DataRoot = char(cfg.DataRoot);
cfg.OutputDir = char(cfg.OutputDir);
cfg.PhaseAlignment = lower(char(cfg.PhaseAlignment));
cfg.MetricCellSubset = lower(char(cfg.MetricCellSubset));
cfg.HeadingUnits = lower(char(cfg.HeadingUnits));
cfg.BehaviorFilePattern = char(cfg.BehaviorFilePattern);
cfg.MorphOrder = cellstr(string(cfg.MorphOrder(:))');
cfg.ClassOrder = cellstr(string(cfg.ClassOrder(:))');
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
if ~ismember(cfg.HeadingUnits, {'auto','radians','degrees'})
    error('HeadingUnits must be ''auto'', ''radians'' or ''degrees''.');
end
if numel(cfg.ClassOrder) ~= 3 || ...
        ~all(ismember(["CW","CCW","Symmetric"],string(cfg.ClassOrder)))
    error('ClassOrder must contain CW, CCW and Symmetric exactly once.');
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
    cfg.OutputDir = fullfile(pwd, ['CW_CCW_Symmetric_comparison_' stamp]);
end
if exist(cfg.OutputDir, 'dir') ~= 7
    mkdir(cfg.OutputDir);
end

fprintf('\n=== CW / CCW / Symmetric cross-morph analysis ===\n');
fprintf('Output directory:\n  %s\n', cfg.OutputDir);
fprintf('Phase alignment mode: %s\n', cfg.PhaseAlignment);
fprintf(['Phase-tuned criterion: ORI_V15 information shuffle, %d shifts/cell, ' ...
    'within-fish BH-FDR q < %.3g\n'], round(cfg.PhaseTuningNShuffles), ...
    cfg.PhaseTuningFDRAlpha);

%% Discover one newest CLASSIFIED candidate file per fish
fileTable = discoverClassifiedCandidateFiles(cfg);
if isempty(fileTable)
    error(['No *_candidate_neurons_clean_classified.mat files were found. ' ...
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

    sourceFile = char(fileTable.File(iFish));
    S = load(sourceFile);
    classLabels = extractClassificationFromStruct(S, sourceFile, cfg.ClassOrder);
    F = analyseOneFish(S, char(fileTable.Morph(iFish)), ...
        char(fileTable.Fish(iFish)), sourceFile, classLabels, ...
        commonCenters, commonEdges, cfg);

    fishData(iFish) = F;
    cellTables{iFish} = makeCellTable(F);
    alignmentTables{iFish} = makeAlignmentTable(F);

    fprintf(['  %d candidate cells | phase offset %+0.1f deg | ' ...
        '%d phase-tuned | PCA-phase fit R = %.3f | heading: %s\n'], F.nCells, ...
        rad2deg(F.phaseOffsetRad), nnz(F.PhaseTuned), F.pcaPhaseFitR, ...
        char(F.headingSource));
end

cellTable = vertcat(cellTables{:});
alignmentTable = vertcat(alignmentTables{:});
featureMethodAudit = buildORIFeatureMethodAudit(cellTable);

% Figure 4E-H can reproduce the poster phase-tuned-cell subset explicitly.
[metricCellTable, subsetDescription] = selectMetricCellSubset(cellTable, cfg);
fprintf('\nClass Figure 4E-H cell subset: %s\n', subsetDescription);
fprintf('Retained %d / %d candidate cells for metric panels.\n', ...
    height(metricCellTable), height(cellTable));

%% Class-focused summaries
metricDefs = posterMetricDefinitions();
classFishSummary = buildClassFishSummary(metricCellTable, ...
    metricDefs, cfg.NBootstrap, cfg.ClassOrder);
classStatsTable = runClassMorphStatistics(classFishSummary, metricDefs, cfg);

highlyTunedCellTable = cellTable(logical(cellTable.PhaseTuned) & ...
    isfinite(cellTable.PhaseTuningQValue) & ...
    cellTable.PhaseTuningQValue < cfg.PhaseTuningFDRAlpha,:);
highlyTunedSubsetDescription = sprintf(['Highly tuned cells within each fish: ' ...
    'ORI_V15 circular-shift information test, BH-FDR q < %.3g'], ...
    cfg.PhaseTuningFDRAlpha);
highlyTunedClassFishSummary = buildClassFishSummary( ...
    highlyTunedCellTable, metricDefs, cfg.NBootstrap, cfg.ClassOrder);
highlyTunedClassStatsTable = runClassMorphStatistics( ...
    highlyTunedClassFishSummary, metricDefs, cfg);

preferredDirectionBins = buildPreferredDirectionBinTable(fishData, cfg, false);
highlyTunedPreferredDirectionBins = ...
    buildPreferredDirectionBinTable(fishData, cfg, true);
pairwiseCorrelationBins = vertcat(fishData.pairwiseCorrelationBins);
highlyTunedPairwiseCorrelationBins = ...
    vertcat(fishData.pairwiseCorrelationBinsHighlyTuned);
phaseHeadingFish = buildPhaseHeadingFishTable(fishData, cfg, false);
highlyTunedPhaseHeadingFish = ...
    buildPhaseHeadingFishTable(fishData, cfg, true);

%% Figures
morphColors = [ ...
    0.88 0.25 0.28;  % Surface: red
    0.25 0.72 0.38;  % Molino: green
    0.95 0.78 0.10]; % Pachon: yellow

figPaths = struct();
figPaths.preferredDirections = plotPreferredDirectionsByClass( ...
    preferredDirectionBins, cfg, false);
figPaths.preferredDirectionsHighlyTuned = plotPreferredDirectionsByClass( ...
    highlyTunedPreferredDirectionBins, cfg, true);
figPaths.classPCProjection = plotClassPCProjection( ...
    fishData, cfg, morphColors, false);
figPaths.classPCProjectionHighlyTuned = plotClassPCProjection( ...
    fishData, cfg, morphColors, true);
figPaths.pairwiseCorrelationByPhase = plotPairwiseCorrelationByPhase( ...
    pairwiseCorrelationBins, cfg, false);
figPaths.pairwiseCorrelationByPhaseHighlyTuned = ...
    plotPairwiseCorrelationByPhase( ...
    highlyTunedPairwiseCorrelationBins, cfg, true);
figPaths.phaseHeadingCorrelation = plotClassPhaseHeadingCorrelation( ...
    phaseHeadingFish, cfg, false);
figPaths.phaseHeadingCorrelationHighlyTuned = ...
    plotClassPhaseHeadingCorrelation( ...
    highlyTunedPhaseHeadingFish, cfg, true);
figPaths.phaseHeadingOverlays = plotClassPhaseHeadingOverlays( ...
    fishData, phaseHeadingFish, cfg, false);
figPaths.phaseHeadingOverlaysHighlyTuned = plotClassPhaseHeadingOverlays( ...
    fishData, highlyTunedPhaseHeadingFish, cfg, true);
figPaths.classMetrics = plotClassMetricPanels(metricCellTable, ...
    classFishSummary, classStatsTable, metricDefs, cfg, morphColors, false);
figPaths.classMetricsHighlyTuned = plotClassMetricPanels( ...
    highlyTunedCellTable, highlyTunedClassFishSummary, ...
    highlyTunedClassStatsTable, metricDefs, cfg, morphColors, true);

%% Return numeric results in memory; only PNG figures are written.
results = struct();
results.config = cfg;
results.inputFiles = fileTable;
results.cells = cellTable;
results.metricCells = metricCellTable;
results.metricCellSubsetDescription = subsetDescription;
results.highlyTunedCells = highlyTunedCellTable;
results.highlyTunedSubsetDescription = highlyTunedSubsetDescription;
results.phaseAlignment = alignmentTable;
results.ORIFeatureMethodAudit = featureMethodAudit;
results.classFishSummary = classFishSummary;
results.classStatistics = classStatsTable;
results.highlyTunedClassFishSummary = highlyTunedClassFishSummary;
results.highlyTunedClassStatistics = highlyTunedClassStatsTable;
results.preferredDirectionBins = preferredDirectionBins;
results.highlyTunedPreferredDirectionBins = highlyTunedPreferredDirectionBins;
results.pairwiseCorrelationBins = pairwiseCorrelationBins;
results.highlyTunedPairwiseCorrelationBins = ...
    highlyTunedPairwiseCorrelationBins;
results.classPhaseHeading = phaseHeadingFish;
results.highlyTunedClassPhaseHeading = highlyTunedPhaseHeadingFish;
results.phaseBinCentersDeg = commonCentersDeg;
results.fish = fishData;
results.figurePaths = figPaths;

fprintf('\nAnalysis complete.\n');
fprintf(['Only classified candidate files were used. All candidates define ' ...
    'the common PC ring; class trajectories use class-specific activity.\n']);
fprintf(['Fish was the summary/inferential unit; each class x metric has one ' ...
    'three-morph KW test and two unadjusted planned KW comparisons.\n']);
fprintf('PNG figures saved in:\n  %s\n', cfg.OutputDir);

end

%% ========================================================================
function fileTable = discoverClassifiedCandidateFiles(cfg)
% Return one newest classified file per morph + recording.

if isempty(cfg.Files)
    allPaths = {};
    for i = 1:numel(cfg.MorphOrder)
        morphLower = lower(cfg.MorphOrder{i});
        root = fullfile(cfg.DataRoot, morphLower);
        if exist(root, 'dir') ~= 7
            warning('Morph directory not found: %s', root);
            continue;
        end
        d = dir(fullfile(root, '**', '*_candidate_neurons_clean_classified.mat'));
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
    if isempty(regexp(allPaths{i}, '_candidate_neurons_clean_classified\.mat$', 'once'))
        error(['Every explicit input must be a ' ...
            '*_candidate_neurons_clean_classified.mat file:\n  %s'], allPaths{i});
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
tok = regexp(name, ['_(\d{8})_(\d{6})_candidate_neurons_clean' ...
    '_classified$'], 'tokens', 'once');
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
function F = analyseOneFish(S, morph, fish, filePath, classLabels, ...
    commonCenters, commonEdges, cfg)

required = {'candidate_cell_ids','calcium_traces','network_phase_rad'};
for i = 1:numel(required)
    if ~isfield(S, required{i})
        error('%s is missing required variable "%s".', filePath, required{i});
    end
end

cellIds = double(S.candidate_cell_ids(:));
nCells = numel(cellIds);
if numel(classLabels) ~= nCells
    error('%s: classification length does not match candidate_cell_ids.', filePath);
end
traces = orientTraceMatrix(double(S.calcium_traces), nCells, filePath);
classifiedKeep = ismember(classLabels, string(cfg.ClassOrder));
if ~any(classifiedKeep)
    error('%s has no fitted CW, CCW or Symmetric cells.', filePath);
end
if any(~classifiedKeep)
    fprintf('  Restricting class analysis to %d/%d fitted phase-tuned cells.', ...
        nnz(classifiedKeep),numel(classifiedKeep));
    fprintf('%s', newline);
end
cellIds = cellIds(classifiedKeep);
classLabels = classLabels(classifiedKeep);
traces = traces(:,classifiedKeep);
nCells = numel(cellIds);
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
    extractCanonicalORIPhaseCurves(S, nCells, phaseTest, classifiedKeep);

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

% Recompute one rPC population-vector trajectory per class using only the
% cells of that class, while keeping the common all-candidate PC ring.
[headingAtImaging, headingSource] = extractHeadingAtImagingTimes( ...
    S, filePath, time, cfg);
classPhase = computeAllClassPhases(traces, alignedPcaXY, classLabels, ...
    time, headingAtImaging, cfg);

% Correlations and phase differences are calculated within fish. Only the
% binned fish summaries are retained to avoid an enormous pair-level file.
pairwiseCorrelationBins = computePairwiseCorrelationBins(traces, ...
    preferredAligned, classLabels, morph, fish, cfg);

% Recompute trajectory and pairwise summaries from the highly tuned cells
% within this fish. The common PCA coordinates remain those established by
% all classified candidates so the two figure versions share one geometry.
highlyTuned = logical(phaseTest.IsTuned(:)) & ...
    isfinite(phaseTest.QValue(:)) & ...
    phaseTest.QValue(:) < cfg.PhaseTuningFDRAlpha;
classPhaseHighlyTuned = computeAllClassPhases( ...
    traces(:,highlyTuned), alignedPcaXY(highlyTuned,:), ...
    classLabels(highlyTuned), time, headingAtImaging, cfg);
pairwiseCorrelationBinsHighlyTuned = computePairwiseCorrelationBins( ...
    traces(:,highlyTuned), preferredAligned(highlyTuned), ...
    classLabels(highlyTuned), morph, fish, cfg);

F = emptyFishStruct();
F.morph = string(morph);
F.fish = string(fish);
F.file = string(filePath);
F.nCells = nCells;
F.cellIds = cellIds;
F.ClassLabel = classLabels(:);
F.HasClassification = true(nCells,1);
F.classifiedFile = string(filePath);
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
F.headingSource = string(headingSource);
F.classPhase = classPhase;
F.pairwiseCorrelationBins = pairwiseCorrelationBins;
F.classPhaseHighlyTuned = classPhaseHighlyTuned;
F.pairwiseCorrelationBinsHighlyTuned = pairwiseCorrelationBinsHighlyTuned;
end

function labels = extractClassificationFromStruct(S, filePath, classOrder)
% Read class labels directly from the classified candidate file. Both a
% label vector and the index-vector format produced by the classifier are
% accepted. No clean sibling file is consulted.

if ~isfield(S,'candidate_cell_ids')
    error('%s has no candidate_cell_ids.', filePath);
end
cellIds = double(S.candidate_cell_ids(:));
n = numel(cellIds);
if numel(unique(cellIds)) ~= n
    error('%s contains duplicate candidate_cell_ids.', filePath);
end

containers = {S};
if isfield(S,'neuron_classification') && isstruct(S.neuron_classification)
    containers{end+1} = S.neuron_classification; %#ok<AGROW>
end

rawLabels = [];
sourceIds = [];
for k = 1:numel(containers)
    C = containers{k};
    rawLabels = firstFieldCaseInsensitive(C, ...
        {'neuron_class_label','neuron_class_labels','labels','class_label','class_labels'});
    if ~isempty(rawLabels)
        sourceIds = firstFieldCaseInsensitive(C, ...
            {'candidate_cell_ids','cell_ids','candidate_ids'});
        break;
    end
end

if ~isempty(rawLabels)
    rawLabels = string(rawLabels(:));
    if ~isempty(sourceIds)
        sourceIds = double(sourceIds(:));
        if numel(sourceIds) ~= numel(rawLabels) || ...
                numel(unique(sourceIds)) ~= numel(sourceIds)
            error('%s has inconsistent classification labels and IDs.', filePath);
        end
        [matched,pos] = ismember(cellIds,sourceIds);
        if ~all(matched)
            error('%s does not classify every candidate_cell_id.', filePath);
        end
        rawLabels = rawLabels(pos);
    elseif numel(rawLabels) ~= n
        error('%s: class label count does not match candidate count.', filePath);
    end
    labels = normalizeClassLabels(rawLabels);
else
    labels = repmat("",n,1);
    aliases = { ...
        {'cw_candidate_indices','cw_indices','cw_idx','idx_cw','cw_mask','cw_cell_ids','cw_ids'}, ...
        {'ccw_candidate_indices','ccw_indices','ccw_idx','idx_ccw','ccw_mask','ccw_cell_ids','ccw_ids', ...
         'cww_candidate_indices','cww_indices','cww_idx','idx_cww'}, ...
        {'symmetric_candidate_indices','symmetric_indices','symmetric_idx', ...
         'idx_symmetric','symmetric_mask','symmetric_cell_ids','symmetric_ids', ...
         'sym_indices','sym_idx','idx_sym'}};
    canonical = ["CW","CCW","Symmetric"];
    for c = 1:numel(canonical)
        [values,fieldName] = firstFieldAcrossContainers(containers,aliases{c});
        if isempty(values)
            continue;
        end
        idx = classificationValuesToCandidateIndices(values,fieldName,cellIds);
        if any(strlength(labels(idx)) > 0)
            error('%s has overlapping class indices.', filePath);
        end
        labels(idx) = canonical(c);
    end
end

valid = ismember(labels,string(classOrder));
ignored = ismember(labels,["NotPhaseTuned","NotFit"]);
if any(~valid & ~ignored)
    bad = unique(labels(~valid & ~ignored));
    error(['%s has %d unrecognized candidates. Labels: %s.'], ...
        filePath, nnz(~valid & ~ignored), strjoin(cellstr(bad),', '));
end
fprintf(['  Classification: %d CW, %d CCW, %d Symmetric; ' ...
    '%d cells not fitted/phase-tuned and ignored.'], ...
    nnz(labels == "CW"),nnz(labels == "CCW"), ...
    nnz(labels == "Symmetric"),nnz(ignored));
fprintf('%s', newline);
end

function value = firstFieldCaseInsensitive(S, aliases)
value = [];
if ~isstruct(S) || ~isscalar(S)
    return;
end
names = fieldnames(S);
for i = 1:numel(aliases)
    hit = find(strcmpi(names,aliases{i}),1);
    if ~isempty(hit)
        value = S.(names{hit});
        return;
    end
end
end

function [value,fieldName] = firstFieldAcrossContainers(containers,aliases)
value = [];
fieldName = '';
for k = 1:numel(containers)
    C = containers{k};
    names = fieldnames(C);
    for i = 1:numel(aliases)
        hit = find(strcmpi(names,aliases{i}),1);
        if ~isempty(hit)
            value = C.(names{hit});
            fieldName = names{hit};
            return;
        end
    end
end
end

function labels = normalizeClassLabels(rawLabels)
rawLabels = string(rawLabels(:));
labels = repmat("",size(rawLabels));
for i = 1:numel(rawLabels)
    key = lower(strtrim(rawLabels(i)));
    key = regexprep(key,'[^a-z]','');
    if startsWith(key,"ccw") || startsWith(key,"cww") || ...
            startsWith(key,"counterclockwise") || startsWith(key,"leftshifter")
        labels(i) = "CCW";
    elseif startsWith(key,"cw") || startsWith(key,"clockwise") || ...
            startsWith(key,"rightshifter")
        labels(i) = "CW";
    elseif startsWith(key,"symmetric") || startsWith(key,"symmetry") || ...
            key == "sym"
        labels(i) = "Symmetric";
    elseif startsWith(key,"notphasetuned")
        labels(i) = "NotPhaseTuned";
    elseif startsWith(key,"notfit")
        labels(i) = "NotFit";
    end
end
end

function idx = classificationValuesToCandidateIndices(values,fieldName,cellIds)
n = numel(cellIds);
if islogical(values)
    if numel(values) ~= n
        error('Classification mask %s must have %d elements.',fieldName,n);
    end
    idx = find(values(:));
    return;
end
if ~isnumeric(values)
    error('Classification field %s must be numeric or logical.',fieldName);
end
v = double(values(:));
v = v(isfinite(v));
if numel(v) == n && all(v == 0 | v == 1)
    idx = find(logical(v));
    return;
end
if contains(lower(fieldName),'id')
    [matched,idx] = ismember(v,cellIds);
    if ~all(matched)
        error('Classification field %s contains unknown candidate IDs.',fieldName);
    end
else
    if any(v < 1 | v > n | v ~= round(v))
        error(['Classification field %s is neither a valid 1-based ' ...
            'candidate-index vector nor an ID field.'],fieldName);
    end
    idx = v;
end
idx = unique(idx(:),'stable');
end

function F = emptyFishStruct()
F = struct( ...
    'morph', string.empty, 'fish', string.empty, 'file', string.empty, ...
    'nCells', 0, 'cellIds', [], 'ClassLabel', string.empty, ...
    'HasClassification', [], 'classifiedFile', string.empty, ...
    'phaseSign', 1, 'phaseOffsetRad', NaN, ...
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
    'pcaRotationRad', NaN, 'pcaPhaseFitR', NaN, ...
    'headingSource', string.empty, 'classPhase', struct([]), ...
    'pairwiseCorrelationBins', table(), ...
    'classPhaseHighlyTuned', struct([]), ...
    'pairwiseCorrelationBinsHighlyTuned', table());
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

function [headingAtImaging,source] = extractHeadingAtImagingTimes( ...
    S,filePath,timeImaging,cfg)
% Recover heading either from the classified file or, if necessary, from
% the newest pass2 behavior file in the same recording folder.

headingAtImaging = nan(size(timeImaging));
[heading,timeHeading,fieldName,found] = findHeadingSeries(S);
sourcePrefix = 'classified file';
headingStruct = S;

if ~found
    behaviorPath = '';
    if isfield(S,'source_behavior_path') && ...
            isfile(char(string(S.source_behavior_path)))
        behaviorPath = char(string(S.source_behavior_path));
    else
        folder = fileparts(filePath);
        d = dir(fullfile(folder,cfg.BehaviorFilePattern));
        d = d(~[d.isdir]);
        if ~isempty(d)
            [~,order] = sort([d.datenum],'descend');
            d = d(order);
            behaviorPath = fullfile(d(1).folder,d(1).name);
            if numel(d) > 1
                warning('%s: multiple behavior files; using newest: %s', ...
                    filePath,behaviorPath);
            end
        end
    end
    if ~isempty(behaviorPath)
        B = load(behaviorPath);
        [heading,timeHeading,fieldName,found] = findHeadingSeries(B);
        sourcePrefix = behaviorPath;
        headingStruct = B;
    end
end

if ~found
    source = 'unavailable';
    warning('%s: no heading estimate found; heading analyses will be skipped.',filePath);
    return;
end

heading = convertHeadingToRadians(double(heading(:)),fieldName,cfg.HeadingUnits);
heading = unwrapFiniteSegments(heading);

if isempty(timeHeading)
    if numel(heading) == numel(timeImaging)
        timeHeading = timeImaging;
    else
        fps = findBehaviorFPS(headingStruct);
        if ~isfinite(fps)
            fps = 250;
        end
        timeHeading = (0:numel(heading)-1)'./fps;
    end
else
    timeHeading = double(timeHeading(:));
end

n = min(numel(heading),numel(timeHeading));
heading = heading(1:n);
timeHeading = timeHeading(1:n);
good = isfinite(timeHeading) & isfinite(heading);
if nnz(good) < 5
    source = 'unavailable: fewer than 5 finite heading samples';
    return;
end
timeHeading = timeHeading(good);
heading = heading(good);
[timeHeading,uniqueIdx] = unique(timeHeading,'stable');
heading = heading(uniqueIdx);

% If both time bases have different absolute origins but compatible
% durations, rebase them to zero. This handles already-cropped pass2 files.
overlap = min(max(timeImaging),max(timeHeading)) - ...
    max(min(timeImaging),min(timeHeading));
if overlap <= 0
    tIm = timeImaging - timeImaging(1);
    tHd = timeHeading - timeHeading(1);
else
    tIm = timeImaging;
    tHd = timeHeading;
end
headingAtImaging = interp1(tHd,heading,tIm,'linear',NaN);
if nnz(isfinite(headingAtImaging)) < 5
    source = 'unavailable: no imaging/heading time overlap';
    headingAtImaging(:) = NaN;
    return;
end
source = sprintf('%s [%s]',sourcePrefix,fieldName);
end

function [heading,timeHeading,fieldName,found] = findHeadingSeries(S)
heading = [];
timeHeading = [];
fieldName = '';
found = false;

% Highest-priority source: ORI's exported overlap already has an explicit
% time base and is exactly the series used by its phase-heading analysis.
overlayContainers = {};
if isfield(S,'heading_phase_global_export') && ...
        isstruct(S.heading_phase_global_export)
    overlayContainers{end+1} = S.heading_phase_global_export; %#ok<AGROW>
end
if isfield(S,'heading_overlay') && isstruct(S.heading_overlay)
    overlayContainers{end+1} = S.heading_overlay; %#ok<AGROW>
end
if isfield(S,'hdMetrics') && isstruct(S.hdMetrics) && ...
        isfield(S.hdMetrics,'r1pi') && isfield(S.hdMetrics.r1pi,'phase')
    tags = fieldnames(S.hdMetrics.r1pi.phase);
    for t = 1:numel(tags)
        P = S.hdMetrics.r1pi.phase.(tags{t});
        if isstruct(P) && isfield(P,'heading_overlay') && isstruct(P.heading_overlay)
            overlayContainers{end+1} = P.heading_overlay; %#ok<AGROW>
        end
    end
end
for i = 1:numel(overlayContainers)
    H = overlayContainers{i};
    heading = firstFieldCaseInsensitive(H,{'heading_overlap','heading_est','heading'});
    timeHeading = firstFieldCaseInsensitive(H,{'T_overlap','time_overlap','t_beh','time_s'});
    if ~isempty(heading)
        fieldName = 'heading_overlap_rad';
        found = true;
        return;
    end
end

containers = {S};
nestedNames = {'pass2FileResult','swimResults_pass2','swimResults','result', ...
    'behavior','beh','headingContinuous'};
for i = 1:numel(nestedNames)
    if isfield(S,nestedNames{i}) && isstruct(S.(nestedNames{i}))
        containers{end+1} = S.(nestedNames{i}); %#ok<AGROW>
    end
end

headingAliases = {'headingRadContinuous','heading_est_imaging_rad', ...
    'heading_imaging_rad','heading_estimate_rad','heading_est_imaging', ...
    'continuousBoutHeadingRad','continuousBoutHeading', ...
    'continuous_heading','heading_est','heading_estimate','heading'};
timeAliases = {'heading_time_s','t_beh','behavior_time_s','time_sec','time_s','Tsec'};
for k = 1:numel(containers)
    C = containers{k};
    names = fieldnames(C);
    for a = 1:numel(headingAliases)
        hit = find(strcmpi(names,headingAliases{a}),1);
        if isempty(hit)
            continue;
        end
        candidate = C.(names{hit});
        if ~isnumeric(candidate) || ~isvector(candidate) || numel(candidate) < 5
            continue;
        end
        heading = candidate;
        fieldName = names{hit};
        timeHeading = firstFieldCaseInsensitive(C,timeAliases);
        if ~isempty(timeHeading) && numel(timeHeading) ~= numel(heading)
            timeHeading = [];
        end
        found = true;
        return;
    end
end
end

function fps = findBehaviorFPS(S)
fps = NaN;
containers = {S};
nestedNames = {'pass2FileResult','swimResults_pass2','swimResults','result','behavior','beh'};
for i = 1:numel(nestedNames)
    if isfield(S,nestedNames{i}) && isstruct(S.(nestedNames{i}))
        containers{end+1} = S.(nestedNames{i}); %#ok<AGROW>
    end
end
for k = 1:numel(containers)
    value = firstFieldCaseInsensitive(containers{k}, ...
        {'fps_beh','behavior_fps','fps','frame_rate_hz'});
    if isnumeric(value) && isscalar(value) && isfinite(value) && value > 0
        fps = double(value);
        return;
    end
end
end

function headingRad = convertHeadingToRadians(heading,fieldName,mode)
switch mode
    case 'degrees'
        headingRad = deg2rad(heading);
    case 'radians'
        headingRad = heading;
    otherwise
        low = lower(fieldName);
        if contains(low,'deg')
            headingRad = deg2rad(heading);
        elseif contains(low,'rad') || contains(low,'heading_est') || ...
                contains(low,'heading_overlap')
            headingRad = heading;
        else
            finiteHeading = heading(isfinite(heading));
            if ~isempty(finiteHeading) && ...
                    percentileLocal(abs(finiteHeading),95) > 3*pi
                headingRad = deg2rad(heading);
            else
                headingRad = heading;
            end
        end
end
end

function xUnwrapped = unwrapFiniteSegments(x)
x = double(x(:));
xUnwrapped = nan(size(x));
finiteMask = isfinite(x);
starts = find(finiteMask & [true; ~finiteMask(1:end-1)]);
stops = find(finiteMask & [~finiteMask(2:end); true]);
for k = 1:numel(starts)
    idx = starts(k):stops(k);
    xUnwrapped(idx) = unwrap(x(idx));
end
end

function classPhase = computeAllClassPhases(traces,alignedPcaXY,classLabels, ...
    time,headingAtImaging,cfg)
template = struct('Class',string.empty,'NCells',0,'NValidFrames',0, ...
    'SignedR',NaN,'AbsoluteR',NaN,'Time',time(:), ...
    'PhaseZ',nan(size(time(:))),'HeadingZ',nan(size(time(:))), ...
    'VectorMagnitude',nan(size(time(:))),'ConfidenceThreshold',NaN);
classPhase = repmat(template,numel(cfg.ClassOrder),1);

for c = 1:numel(cfg.ClassOrder)
    className = string(cfg.ClassOrder{c});
    use = classLabels(:) == className & all(isfinite(alignedPcaXY),2);
    classPhase(c).Class = className;
    classPhase(c).NCells = nnz(use);
    if nnz(use) < cfg.MinimumClassCellsForPhase
        continue;
    end

    X = fillmissing(double(traces(:,use)),'linear',1,'EndValues','nearest');
    X(~isfinite(X)) = 0;
    mu = mean(X,1);
    sd = std(X,0,1);
    sd(~isfinite(sd) | sd <= 0) = 1;
    X = bsxfun(@rdivide,bsxfun(@minus,X,mu),sd);
    X = clipCenterFramewiseLocal(X,2,98);

    alpha = atan2(alignedPcaXY(use,2),alignedPcaXY(use,1));
    Ravg = [cos(alpha(:)),sin(alpha(:))];
    V = (X*Ravg)./nnz(use);
    phaseWrapped = atan2(V(:,2),V(:,1));
    phaseUnwrapped = unwrapFiniteSegments(phaseWrapped);
    vmag = sqrt(sum(V.^2,2));
    threshold = percentileLocal(vmag(isfinite(vmag)),cfg.PhaseConfidencePercentile);
    valid = isfinite(phaseUnwrapped) & isfinite(headingAtImaging) & ...
        isfinite(vmag) & vmag > threshold;

    phaseZ = zscoreUsingMask(phaseUnwrapped,valid);
    headingZ = zscoreUsingMask(headingAtImaging,valid);
    phaseZ(~valid) = NaN;
    headingZ(~valid) = NaN;

    classPhase(c).NValidFrames = nnz(valid);
    classPhase(c).Time = time(:);
    classPhase(c).PhaseZ = phaseZ;
    classPhase(c).HeadingZ = headingZ;
    classPhase(c).VectorMagnitude = vmag;
    classPhase(c).ConfidenceThreshold = threshold;
    if nnz(valid) >= 10 && std(phaseZ(valid)) > 0 && std(headingZ(valid)) > 0
        r = corr(phaseZ(valid),headingZ(valid),'Rows','complete','Type','Pearson');
        classPhase(c).SignedR = r;
        classPhase(c).AbsoluteR = abs(r);
    end
end
end

function z = zscoreUsingMask(x,mask)
x = double(x(:));
mask = logical(mask(:)) & isfinite(x);
z = nan(size(x));
if nnz(mask) < 2
    return;
end
mu = mean(x(mask));
sd = std(x(mask));
if ~isfinite(sd) || sd <= 0
    return;
end
z = (x-mu)./sd;
end

function T = computePairwiseCorrelationBins(traces,preferredPhase,classLabels, ...
    morph,fish,cfg)
pairA = ["CW","CCW","Symmetric","CW","CW","CCW"];
pairB = ["CW","CCW","Symmetric","CCW","Symmetric","Symmetric"];
pairLabel = ["CW-CW","CCW-CCW","Symmetric-Symmetric", ...
    "CW-CCW","CW-Symmetric","CCW-Symmetric"];
edges = linspace(0,pi,round(cfg.PairPhaseDifferenceNBins)+1);
centersDeg = rad2deg((edges(1:end-1)+edges(2:end))./2);

C = corr(double(traces),'Rows','pairwise','Type','Pearson');
rows = cell(numel(pairLabel),1);
for p = 1:numel(pairLabel)
    idxA = find(classLabels(:) == pairA(p));
    idxB = find(classLabels(:) == pairB(p));
    if pairA(p) == pairB(p)
        [ia,ib] = find(triu(true(numel(idxA)),1));
        cellA = idxA(ia);
        cellB = idxA(ib);
    else
        [ia,ib] = ndgrid(idxA,idxB);
        cellA = ia(:);
        cellB = ib(:);
    end
    if isempty(cellA)
        dPhase = [];
        absCorrelation = [];
    else
        dPhase = abs(wrapToPiLocal(preferredPhase(cellA)-preferredPhase(cellB)));
        absCorrelation = abs(C(sub2ind(size(C),cellA,cellB)));
    end
    bin = discretize(dPhase,edges);
    meanAbs = nan(numel(centersDeg),1);
    nPairs = zeros(numel(centersDeg),1);
    for b = 1:numel(centersDeg)
        good = bin == b & isfinite(absCorrelation);
        nPairs(b) = nnz(good);
        if any(good)
            meanAbs(b) = mean(absCorrelation(good));
        end
    end
    rows{p} = table(repmat(string(morph),numel(centersDeg),1), ...
        repmat(string(fish),numel(centersDeg),1), ...
        repmat(pairLabel(p),numel(centersDeg),1),centersDeg(:), ...
        meanAbs,nPairs,repmat(numel(idxA),numel(centersDeg),1), ...
        repmat(numel(idxB),numel(centersDeg),1), ...
        'VariableNames',{'Morph','Fish','PairType','PhaseDifferenceCenterDeg', ...
        'MeanAbsCorrelation','NPairs','NCellsClass1','NCellsClass2'});
end
T = vertcat(rows{:});
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
        extractCanonicalORIPhaseCurves(S, nCells, phaseTest, classifiedKeep)
% Candidate exports contain ORI_V15/hdMetrics.hd_single.phase.rate_bins.
% Preserve those exact curves when present because phase_ok is not always
% exported. The local Step 17 result is the principled fallback.
if isfield(S,'tuning_curves_phi_all_cells')
    curves = double(S.tuning_curves_phi_all_cells);
    nAllCandidates = numel(classifiedKeep);
    if size(curves,1) == nAllCandidates
        curves = curves(classifiedKeep,:);
    elseif size(curves,2) == nAllCandidates
        curves = curves(:,classifiedKeep)';
    elseif size(curves,1) == nCells
        % Already restricted to fitted cells.
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
    F.cellIds(:), F.ClassLabel(:), logical(F.HasClassification(:)), ...
    repmat(string(F.classifiedFile),n,1), ...
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
    'CandidateCellID','ClassLabel','HasClassification','ClassifiedSourceFile', ...
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

function T = buildPreferredDirectionBinTable(fishData,cfg,highlyTunedOnly)
edges = linspace(-180,180,round(cfg.PreferredDirectionNBins)+1);
centers = (edges(1:end-1)+edges(2:end))./2;
groups = ["All",string(cfg.ClassOrder)];
rows = cell(numel(fishData)*numel(groups),1);
counter = 0;
for i = 1:numel(fishData)
    phaseDeg = rad2deg(fishData(i).preferredPhaseAlignedRad(:));
    if highlyTunedOnly
        subsetMask = logical(fishData(i).PhaseTuned(:));
    else
        subsetMask = true(size(phaseDeg));
    end
    for g = 1:numel(groups)
        if groups(g) == "All"
            use = subsetMask & isfinite(phaseDeg);
        else
            use = subsetMask & fishData(i).ClassLabel(:) == groups(g) & ...
                isfinite(phaseDeg);
        end
        if any(use)
            probability = histcounts(phaseDeg(use),edges,'Normalization','probability')';
        else
            probability = nan(numel(centers),1);
        end
        counter = counter+1;
        rows{counter} = table(repmat(fishData(i).morph,numel(centers),1), ...
            repmat(fishData(i).fish,numel(centers),1), ...
            repmat(groups(g),numel(centers),1),centers(:),probability, ...
            repmat(nnz(use),numel(centers),1), ...
            'VariableNames',{'Morph','Fish','ClassGroup','DirectionCenterDeg', ...
            'Probability','NCells'});
    end
end
T = vertcat(rows{1:counter});
end

function [densityTable, summaryTable] = buildPreferredPhaseKDETables(fishData,cfg)
% Circular Gaussian KDE and circular summary statistics. Each fish is kept
% separate so that animals, rather than neurons, remain the replicate unit.
nGrid = round(cfg.PreferredPhaseKDEGridPoints);
if mod(nGrid,2) == 0
    nGrid = nGrid + 1;
end
gridDeg = linspace(-180,180,nGrid)';
gridRad = deg2rad(gridDeg);
bandwidthRad = deg2rad(cfg.PreferredPhaseKDEBandwidthDeg);

densityRows = cell(numel(fishData)*numel(cfg.ClassOrder),1);
summaryRows = cell(numel(fishData)*numel(cfg.ClassOrder),1);
counter = 0;
for i = 1:numel(fishData)
    allPhase = fishData(i).preferredPhaseAlignedRad(:);
    labels = fishData(i).ClassLabel(:);
    for c = 1:numel(cfg.ClassOrder)
        className = string(cfg.ClassOrder{c});
        use = labels == className & isfinite(allPhase);
        phase = wrapToPiLocal(allPhase(use));
        if numel(phase) < cfg.MinimumPreferredPhaseCellsPerFish
            continue;
        end

        density = circularGaussianKDE(phase,gridRad,bandwidthRad);
        z = mean(exp(1i.*phase));
        meanPhaseDeg = rad2deg(angle(z));
        absoluteCentroidDisplacementDeg = abs(meanPhaseDeg);
        resultantLength = abs(z);
        circularVariance = 1-resultantLength;

        % Fixed-bandwidth KDE entropy: 1 indicates uniform phase coverage;
        % lower values indicate concentration into a smaller ring sector.
        mass = density(:);
        mass = mass ./ sum(mass);
        positive = mass > 0 & isfinite(mass);
        normalizedEntropy = -sum(mass(positive).*log(mass(positive))) ./ ...
            log(numel(mass));

        counter = counter+1;
        densityRows{counter} = table( ...
            repmat(fishData(i).morph,nGrid,1), ...
            repmat(fishData(i).fish,nGrid,1), ...
            repmat(className,nGrid,1), gridDeg, density, ...
            repmat(numel(phase),nGrid,1), ...
            'VariableNames',{'Morph','Fish','Class','PreferredPhaseDeg', ...
            'CircularDensityPerRad','NCells'});
        summaryRows{counter} = table(fishData(i).morph,fishData(i).fish, ...
            className,numel(phase),meanPhaseDeg, ...
            absoluteCentroidDisplacementDeg,resultantLength, ...
            circularVariance,normalizedEntropy, ...
            'VariableNames',{'Morph','Fish','Class','NCells', ...
            'CircularMeanPhaseDeg','AbsoluteCircularCentroidDisplacementDeg', ...
            'ResultantLength','CircularVariance', ...
            'NormalizedKDEEntropy'});
    end
end

if counter == 0
    densityTable = table('Size',[0 6], ...
        'VariableTypes',{'string','string','string','double','double','double'}, ...
        'VariableNames',{'Morph','Fish','Class','PreferredPhaseDeg', ...
        'CircularDensityPerRad','NCells'});
    summaryTable = table('Size',[0 9], ...
        'VariableTypes',{'string','string','string','double','double', ...
        'double','double','double','double'}, ...
        'VariableNames',{'Morph','Fish','Class','NCells', ...
        'CircularMeanPhaseDeg','AbsoluteCircularCentroidDisplacementDeg', ...
        'ResultantLength','CircularVariance', ...
        'NormalizedKDEEntropy'});
else
    densityTable = vertcat(densityRows{1:counter});
    summaryTable = vertcat(summaryRows{1:counter});
end
end

function density = circularGaussianKDE(phaseRad,gridRad,bandwidthRad)
% Gaussian kernel applied to the shortest signed angular distance. This is
% circular at +/-pi and therefore has no edge artefact at +/-180 degrees.
phaseRad = phaseRad(:)';
gridRad = gridRad(:);
delta = wrapToPiLocal(gridRad-phaseRad);
density = mean(exp(-0.5.*(delta./bandwidthRad).^2),2) ./ ...
    (sqrt(2*pi).*bandwidthRad);
area = trapz(gridRad,density);
if isfinite(area) && area > 0
    density = density./area;
end
end

function T = buildPhaseHeadingFishTable(fishData,cfg,highlyTunedOnly)
rows = cell(numel(fishData)*numel(cfg.ClassOrder),1);
counter = 0;
for i = 1:numel(fishData)
    if highlyTunedOnly
        classPhase = fishData(i).classPhaseHighlyTuned;
    else
        classPhase = fishData(i).classPhase;
    end
    for c = 1:numel(classPhase)
        P = classPhase(c);
        counter = counter+1;
        rows{counter} = table(fishData(i).morph,fishData(i).fish,P.Class, ...
            P.NCells,P.NValidFrames,P.SignedR,P.AbsoluteR, ...
            fishData(i).headingSource, ...
            'VariableNames',{'Morph','Fish','Class','NClassCells', ...
            'NValidFrames','SignedR','AbsoluteR','HeadingSource'});
    end
end
if counter == 0
    T = table();
else
    T = vertcat(rows{1:counter});
end
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

function summary = buildClassFishSummary(cellTable, metricDefs, nBootstrap, classOrder)
% One value per fish, class and metric. Neurons never act as independent
% replicates for the cross-morph inference.
if isempty(cellTable)
    summary = table('Size',[0 9], ...
        'VariableTypes',{'string','string','string','string','double','double','double','double','string'}, ...
        'VariableNames',{'Morph','Fish','Class','Metric','NCells','Value','CILow','CIHigh','ClassifiedSourceFile'});
    return;
end

rows = {};
counter = 0;
fishKeys = unique(cellTable.Morph + "|" + cellTable.Fish, 'stable');
for i = 1:numel(fishKeys)
    idxFish = (cellTable.Morph + "|" + cellTable.Fish) == fishKeys(i);
    morph = cellTable.Morph(find(idxFish,1));
    fish = cellTable.Fish(find(idxFish,1));
    source = cellTable.ClassifiedSourceFile(find(idxFish,1));

    for c = 1:numel(classOrder)
        className = string(classOrder{c});
        idxClass = idxFish & cellTable.ClassLabel == className;
        for m = 1:height(metricDefs)
            variable = char(metricDefs.Variable(m));
            x = double(cellTable.(variable)(idxClass));
            x = x(isfinite(x));
            if isempty(x)
                value = NaN;
                ci = [NaN NaN];
            else
                value = median(x);
                ci = bootstrapMedianCI(x, nBootstrap);
            end
            counter = counter + 1;
            rows{counter,1} = table(morph, fish, className, ...
                metricDefs.Variable(m), numel(x), value, ci(1), ci(2), source, ...
                'VariableNames', {'Morph','Fish','Class','Metric','NCells', ...
                'Value','CILow','CIHigh','ClassifiedSourceFile'}); %#ok<AGROW>
        end
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

function statsTable = runClassMorphStatistics(fishSummary, metricDefs, cfg)
if isempty(fishSummary)
    statsTable = table('Size',[0 11], ...
        'VariableTypes',{'string','string','string','string','string','double','double','double','double','double','double'}, ...
        'VariableNames',{'Class','Metric','Test','Group1','Group2','NFish1', ...
        'NFish2','KWStatistic','MedianDifference','CliffsDelta','PValue'});
    return;
end

rows = {};
counter = 0;
for c = 1:numel(cfg.ClassOrder)
    className = string(cfg.ClassOrder{c});
    for m = 1:height(metricDefs)
        metric = metricDefs.Variable(m);
        use = fishSummary.Class == className & fishSummary.Metric == metric & ...
            isfinite(fishSummary.Value) & ...
            fishSummary.NCells >= cfg.MinimumMetricCellsPerFish;
        values = double(fishSummary.Value(use));
        groups = fishSummary.Morph(use);

        % One omnibus KW test across Surface, Molino and Pachon. Requiring
        % every requested morph prevents a two-morph subset from being
        % mislabeled as the three-morph omnibus test.
        [pOmnibus, hOmnibus, nByMorph] = kruskalWallisForGroups( ...
            values, groups, string(cfg.MorphOrder));
        counter = counter + 1;
        rows{counter,1} = table(className, metric, ...
            "omnibus_three_group_KW", "All morphs", "", ...
            sum(nByMorph), NaN, hOmnibus, NaN, NaN, pOmnibus, ...
            'VariableNames', {'Class','Metric','Test','Group1','Group2', ...
            'NFish1','NFish2','KWStatistic','MedianDifference', ...
            'CliffsDelta','PValue'}); %#ok<AGROW>

        % Prespecified two-group KW tests. These raw p-values are reported
        % regardless of the omnibus result and receive no Holm correction.
        for k = 1:size(cfg.PlannedContrasts,1)
            g1 = string(cfg.PlannedContrasts{k,1});
            g2 = string(cfg.PlannedContrasts{k,2});
            x1 = values(strcmpi(groups,g1));
            x2 = values(strcmpi(groups,g2));
            pairValues = [x1(:); x2(:)];
            pairGroups = [repmat(g1,numel(x1),1); repmat(g2,numel(x2),1)];
            [pPair, hStatistic] = kruskalWallisForGroups( ...
                pairValues, pairGroups, [g1 g2]);
            medianDifference = median(x2,'omitnan') - median(x1,'omitnan');
            cliffsDelta = cliffsDeltaLocal(x2, x1);
            counter = counter + 1;
            rows{counter,1} = table(className, metric, ...
                "planned_two_group_KW_unadjusted", ...
                g1, g2, numel(x1), numel(x2), hStatistic, medianDifference, ...
                cliffsDelta, pPair, 'VariableNames', {'Class','Metric','Test', ...
                'Group1','Group2','NFish1','NFish2','KWStatistic', ...
                'MedianDifference','CliffsDelta','PValue'}); %#ok<AGROW>
        end
    end
end
statsTable = vertcat(rows{:});
end

function [pValue, hStatistic, nByGroup] = kruskalWallisForGroups(values, groups, requiredGroups)
% Run MATLAB's standard asymptotic Kruskal-Wallis test without its figure.
% The function returns NaN unless every requested group has data.
values = double(values(:));
groups = string(groups(:));
requiredGroups = string(requiredGroups(:));

keep = isfinite(values) & ~ismissing(groups);
values = values(keep);
groups = groups(keep);
nByGroup = zeros(numel(requiredGroups),1);
selected = false(size(values));
for g = 1:numel(requiredGroups)
    thisGroup = strcmpi(groups, requiredGroups(g));
    nByGroup(g) = nnz(thisGroup);
    selected = selected | thisGroup;
end

pValue = NaN;
hStatistic = NaN;
if any(nByGroup == 0)
    return;
end

values = values(selected);
groups = groups(selected);
try
    [pValue, anovaTable] = kruskalwallis(values, cellstr(groups), 'off');
    hStatistic = anovaTable{2,4};
catch ME
    warning('Kruskal-Wallis test failed: %s', ME.message);
    pValue = NaN;
    hStatistic = NaN;
end
end

function statsTable = runLegacyPooledCellDiagnostic(cellTable, metricDefs, cfg)
% Diagnostic only. This deliberately ignores fish nesting and therefore
% must not be reported as the confirmatory morph test. It is saved so that
% an old neuron-pooled poster analysis can be identified explicitly.
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
function pathOut = plotPreferredDirectionsByClass(T,cfg,highlyTunedOnly)
[subsetTitle,fileSuffix] = figureSubsetMetadata(highlyTunedOnly);
groups = ["All",string(cfg.ClassOrder)];
nRows = numel(cfg.MorphOrder);
nCols = numel(groups);
fig = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Preferred firing direction by class','Position',[20 20 1800 1100]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
if strcmp(cfg.PhaseAlignment,'preserve')
    phaseDescription = 'Upstream anatomy-aligned phase convention preserved';
else
    phaseDescription = 'Relative phase: each fish all-candidate peak aligned to 0 deg';
end
title(tl,{['Preferred firing direction by functional class - ' subsetTitle], ...
    sprintf('Each fish normalized before morph averaging; %s',phaseDescription)});
finiteProbability = T.Probability(isfinite(T.Probability));
if isempty(finiteProbability)
    yMaximum = 1;
else
    yMaximum = min(1,max(0.20,1.08.*max(finiteProbability)));
end

for m = 1:nRows
    for g = 1:nCols
        ax = nexttile(tl,(m-1)*nCols+g);
        hold(ax,'on');
        use = strcmpi(T.Morph,cfg.MorphOrder{m}) & T.ClassGroup == groups(g);
        fish = unique(T.Fish(use),'stable');
        x = unique(T.DirectionCenterDeg(use));
        x = sort(double(x(:)))';
        Y = nan(numel(fish),numel(x));
        for f = 1:numel(fish)
            rows = use & T.Fish == fish(f);
            [~,order] = sort(T.DirectionCenterDeg(rows));
            values = T.Probability(rows);
            Y(f,:) = values(order);
        end
        color = classDisplayColor(groups(g));
        pale = 0.82.*[1 1 1]+0.18.*color;
        for f = 1:size(Y,1)
            plot(ax,x,Y(f,:),'Color',pale,'LineWidth',0.8);
        end
        plotMeanSemBand(ax,x,Y,color);
        xlim(ax,[-180 180]);
        ylim(ax,[0 yMaximum]);
        xticks(ax,[-180 -90 0 90 180]);
        grid(ax,'on'); ax.GridAlpha = 0.10; box(ax,'off');
        if m == 1
            title(ax,sprintf('%s (n=%d fish)',groups(g),numel(fish)));
        else
            title(ax,sprintf('n=%d fish',numel(fish)));
        end
        if g == 1
            ylabel(ax,sprintf('%s\nProbability',cfg.MorphOrder{m}));
        else
            set(ax,'YTickLabel',[]);
        end
        if m == nRows
            xlabel(ax,'Preferred phase (deg)');
        else
            set(ax,'XTickLabel',[]);
        end
    end
end
pathOut = saveFigureBoth(fig,cfg.OutputDir, ...
    ['Preferred_firing_direction_All_CW_CCW_Symmetric_3x4' fileSuffix], ...
    cfg.SaveFigFiles);
end

function pathOut = plotClassPCProjection(fishData,cfg,morphColors,highlyTunedOnly)
[subsetTitle,fileSuffix] = figureSubsetMetadata(highlyTunedOnly);
nRows = numel(cfg.MorphOrder);
nCols = numel(cfg.ClassOrder);
fig = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Class-specific PC projection','Position',[40 30 1450 1250]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
title(tl,{['Candidate-cell PC projection separated by class - ' subsetTitle], ...
    'One aligned, circle-normalized all-candidate PC space per fish; colour = preferred phase'});

for m = 1:nRows
    for c = 1:nCols
        ax = nexttile(tl,(m-1)*nCols+c);
        hold(ax,'on');
        nFish = 0;
        nCells = 0;
        for i = 1:numel(fishData)
            if ~strcmpi(fishData(i).morph,cfg.MorphOrder{m})
                continue;
            end
            use = fishData(i).ClassLabel == string(cfg.ClassOrder{c});
            if highlyTunedOnly
                use = use & logical(fishData(i).PhaseTuned(:));
            end
            xy = fishData(i).alignedPcaXY(use,:);
            phaseDeg = rad2deg(fishData(i).preferredPhaseAlignedRad(use));
            good = all(isfinite(xy),2) & isfinite(phaseDeg);
            if any(good)
                try
                    scatter(ax,xy(good,1),xy(good,2),14,phaseDeg(good),'filled', ...
                        'MarkerFaceAlpha',0.42,'MarkerEdgeAlpha',0.10);
                catch
                    scatter(ax,xy(good,1),xy(good,2),14,phaseDeg(good),'filled');
                end
                nCells = nCells+nnz(good);
            end
            nFish = nFish+1;
        end
        theta = linspace(0,2*pi,361);
        plot(ax,cos(theta),sin(theta),':','Color',[0.25 0.25 0.25]);
        xline(ax,0,':','Color',[0.75 0.75 0.75]);
        yline(ax,0,':','Color',[0.75 0.75 0.75]);
        axis(ax,'equal'); xlim(ax,[-2.2 2.2]); ylim(ax,[-2.2 2.2]);
        caxis(ax,[-180 180]); box(ax,'off');
        if m == 1
            title(ax,sprintf('%s | %d cells',cfg.ClassOrder{c},nCells));
        else
            title(ax,sprintf('%d cells',nCells));
        end
        if c == 1
            ylabel(ax,sprintf('%s\naligned PC2',cfg.MorphOrder{m}), ...
                'Color',morphColors(min(m,size(morphColors,1)),:));
        else
            set(ax,'YTickLabel',[]);
        end
        if m == nRows
            xlabel(ax,'aligned PC1');
        else
            set(ax,'XTickLabel',[]);
        end
    end
end
colormap(fig,hsv(256));
cb = colorbar;
try, cb.Layout.Tile = 'east'; catch, end
cb.Label.String = 'Aligned preferred phase (deg)';
pathOut = saveFigureBoth(fig,cfg.OutputDir, ...
    ['PC_projection_CW_CCW_Symmetric_by_morph_3x3' fileSuffix], ...
    cfg.SaveFigFiles);
end

function pathOut = plotPairwiseCorrelationByPhase(T,cfg,highlyTunedOnly)
[subsetTitle,fileSuffix] = figureSubsetMetadata(highlyTunedOnly);
pairs = ["CW-CW","CCW-CCW","Symmetric-Symmetric", ...
    "CW-CCW","CW-Symmetric","CCW-Symmetric"];
pairColors = [0.12 0.42 0.85; 0.90 0.25 0.18; 0.45 0.45 0.45; ...
    0.58 0.20 0.72; 0.10 0.64 0.62; 0.93 0.55 0.12];
nRows = numel(cfg.MorphOrder);
nCols = numel(pairs);
fig = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Pairwise correlation by phase and class pair', ...
    'Position',[10 30 2200 1050]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
title(tl,{['Pairwise activity-correlation magnitude versus phase difference - ' subsetTitle], ...
    'Thin lines = fish; thick line and band = fish mean +/- SEM (pairs are never pooled across fish)'});

for m = 1:nRows
    for p = 1:nCols
        ax = nexttile(tl,(m-1)*nCols+p);
        hold(ax,'on');
        use = strcmpi(T.Morph,cfg.MorphOrder{m}) & T.PairType == pairs(p);
        fish = unique(T.Fish(use),'stable');
        x = unique(T.PhaseDifferenceCenterDeg(use));
        x = sort(double(x(:)))';
        Y = nan(numel(fish),numel(x));
        for f = 1:numel(fish)
            rows = use & T.Fish == fish(f);
            [~,order] = sort(T.PhaseDifferenceCenterDeg(rows));
            values = T.MeanAbsCorrelation(rows);
            Y(f,:) = values(order);
        end
        pale = 0.82.*[1 1 1]+0.18.*pairColors(p,:);
        for f = 1:size(Y,1)
            plot(ax,x,Y(f,:),'Color',pale,'LineWidth',0.75);
        end
        plotMeanSemBand(ax,x,Y,pairColors(p,:));
        xlim(ax,[0 180]); ylim(ax,[0 1]); xticks(ax,[0 90 180]);
        grid(ax,'on'); ax.GridAlpha = 0.10; box(ax,'off');
        if m == 1
            title(ax,sprintf('%s (n=%d fish)',strrep(pairs(p),'-','-'),numel(fish)));
        else
            title(ax,sprintf('n=%d fish',numel(fish)));
        end
        if p == 1
            ylabel(ax,sprintf('%s\nMean |Pearson r|',cfg.MorphOrder{m}));
        else
            set(ax,'YTickLabel',[]);
        end
        if m == nRows
            xlabel(ax,'|Phase difference| (deg)');
        else
            set(ax,'XTickLabel',[]);
        end
    end
end
pathOut = saveFigureBoth(fig,cfg.OutputDir, ...
    ['Pairwise_abs_correlation_vs_phase_difference_class_pairs_3x6' fileSuffix], ...
    cfg.SaveFigFiles);
end

function pathOut = plotClassPhaseHeadingCorrelation(T,cfg,highlyTunedOnly)
[subsetTitle,fileSuffix] = figureSubsetMetadata(highlyTunedOnly);
fig = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Class rPC phase-heading correlation','Position',[100 40 1050 1200]);
tl = tiledlayout(fig,numel(cfg.MorphOrder),1, ...
    'TileSpacing','compact','Padding','compact');
title(tl,{['Class-specific rPC phase versus estimated heading - ' subsetTitle], ...
    'Absolute Pearson r from unwrapped trajectories; points are fish'});
for m = 1:numel(cfg.MorphOrder)
    ax = nexttile(tl,m); hold(ax,'on');
    for c = 1:numel(cfg.ClassOrder)
        use = strcmpi(T.Morph,cfg.MorphOrder{m}) & ...
            T.Class == string(cfg.ClassOrder{c}) & isfinite(T.AbsoluteR);
        y = double(T.AbsoluteR(use));
        color = classDisplayColor(string(cfg.ClassOrder{c}));
        if ~isempty(y)
            b = bar(ax,c,mean(y),0.62,'FaceColor',color,'EdgeColor',color);
            try, b.FaceAlpha = 0.24; catch, end
            sem = std(y)./sqrt(numel(y));
            errorbar(ax,c,mean(y),sem,'LineStyle','none','Color',[.15 .15 .15], ...
                'LineWidth',1.2,'CapSize',8);
            jitter = linspace(-0.12,0.12,numel(y))';
            scatter(ax,c+jitter,y,48,color,'filled','MarkerEdgeColor','k');
        end
    end
    xlim(ax,[0.45 numel(cfg.ClassOrder)+0.55]); ylim(ax,[0 1]);
    set(ax,'XTick',1:numel(cfg.ClassOrder),'XTickLabel',cfg.ClassOrder);
    ylabel(ax,'|r|'); title(ax,cfg.MorphOrder{m});
    grid(ax,'on'); ax.GridAlpha = 0.10; box(ax,'off');
end
pathOut = saveFigureBoth(fig,cfg.OutputDir, ...
    ['Class_rPC_phase_heading_absolute_correlation_by_morph_3rows' fileSuffix], ...
    cfg.SaveFigFiles);
end

function pathOut = plotClassPhaseHeadingOverlays(fishData,T,cfg,highlyTunedOnly)
[subsetTitle,fileSuffix] = figureSubsetMetadata(highlyTunedOnly);
nRows = numel(cfg.MorphOrder);
nCols = numel(cfg.ClassOrder);
fig = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Representative class phase-heading overlays', ...
    'Position',[30 20 1700 1200]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
title(tl,{sprintf(['Representative z-scored class rPC phase and heading - %s ' ...
    '(first %.0f s of valid overlap)'],subsetTitle,cfg.HeadingOverlayWindowSec), ...
    'Representative fish is closest to the morph x class median |r|; title reports full-overlap |r|'});

for m = 1:nRows
    for c = 1:nCols
        ax = nexttile(tl,(m-1)*nCols+c); hold(ax,'on');
        className = string(cfg.ClassOrder{c});
        use = strcmpi(T.Morph,cfg.MorphOrder{m}) & T.Class == className & ...
            isfinite(T.AbsoluteR);
        rows = find(use);
        if isempty(rows)
            text(ax,0.5,0.5,'Heading unavailable','Units','normalized', ...
                'HorizontalAlignment','center','Color',[.45 .45 .45]);
            axis(ax,'off');
            continue;
        end
        values = T.AbsoluteR(rows);
        target = median(values,'omitnan');
        [~,pick] = min(abs(values-target));
        row = rows(pick);
        fishName = T.Fish(row);
        fishIdx = find(strcmpi([fishData.morph],cfg.MorphOrder{m}) & ...
            [fishData.fish] == fishName,1);
        if isempty(fishIdx)
            axis(ax,'off');
            continue;
        end
        classIdx = find([fishData(fishIdx).classPhase.Class] == className,1);
        if highlyTunedOnly
            P = fishData(fishIdx).classPhaseHighlyTuned(classIdx);
        else
            P = fishData(fishIdx).classPhase(classIdx);
        end
        good = isfinite(P.Time) & isfinite(P.PhaseZ) & isfinite(P.HeadingZ);
        if nnz(good) < 5
            text(ax,0.5,0.5,'Insufficient overlap','Units','normalized', ...
                'HorizontalAlignment','center'); axis(ax,'off'); continue;
        end
        firstTime = P.Time(find(good,1));
        show = good & P.Time >= firstTime & ...
            P.Time <= firstTime+cfg.HeadingOverlayWindowSec;
        t = P.Time(show)-firstTime;
        plot(ax,t,P.PhaseZ(show),'Color',classDisplayColor(className),'LineWidth',1.15);
        plot(ax,t,P.HeadingZ(show),'Color',[0.15 0.15 0.15],'LineWidth',1.05);
        yline(ax,0,':','Color',[.75 .75 .75]);
        grid(ax,'on'); ax.GridAlpha = 0.10; box(ax,'off');
        title(ax,sprintf('%s | %s | |r|=%.3f',className,fishName,T.AbsoluteR(row)));
        if c == 1
            ylabel(ax,sprintf('%s\nz-score',cfg.MorphOrder{m}));
        else
            set(ax,'YTickLabel',[]);
        end
        if m == nRows
            xlabel(ax,'Time in displayed window (s)');
        else
            set(ax,'XTickLabel',[]);
        end
        if m == 1 && c == 1
            legend(ax,{'class rPC phase','heading'},'Location','best');
        end
    end
end
pathOut = saveFigureBoth(fig,cfg.OutputDir, ...
    ['Class_rPC_phase_heading_zscore_representative_overlays_3x3' fileSuffix], ...
    cfg.SaveFigFiles);
end

function [titleText,fileSuffix] = figureSubsetMetadata(highlyTunedOnly)
if highlyTunedOnly
    titleText = 'highly tuned cells within each fish';
    fileSuffix = '_highly_tuned_cells_per_fish';
else
    titleText = 'all classified candidate cells';
    fileSuffix = '';
end
end

function color = classDisplayColor(className)
className = string(className);
if className == "CW"
    color = [0.12 0.42 0.85];
elseif className == "CCW"
    color = [0.90 0.25 0.18];
elseif className == "Symmetric"
    color = [0.48 0.38 0.62];
else
    color = [0.18 0.18 0.18];
end
end

function plotMeanSemBand(ax,x,Y,color)
if isempty(x) || isempty(Y)
    return;
end
mu = mean(Y,1,'omitnan');
n = sum(isfinite(Y),1);
sem = std(Y,0,1,'omitnan')./sqrt(max(n,1));
good = isfinite(x) & isfinite(mu);
if nnz(good) >= 2
    xx = x(good);
    lo = mu(good)-sem(good);
    hi = mu(good)+sem(good);
    fill(ax,[xx fliplr(xx)],[lo fliplr(hi)],color, ...
        'FaceAlpha',0.16,'EdgeColor','none');
end
plot(ax,x,mu,'Color',color,'LineWidth',2.3);
end

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

pathOut = saveFigureBoth(fig, cfg.OutputDir, 'Figure3A_aligned_PC_projection', cfg.SaveFigFiles);
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

pathOut = saveFigureBoth(fig, cfg.OutputDir, 'Individual_cell_HD_tuning_heatmaps', cfg.SaveFigFiles);
end

function pathOut = plotMetricPanels(cellTable, fishSummary, statsTable, defs, cfg, morphColors, mode)
fig = figure('Color','w','Visible',cfg.Visible,'Name','Candidate cell metrics');
tl = tiledlayout(fig, 1, height(defs), 'TileSpacing','compact','Padding','compact');

if strcmp(mode,'poster')
    title(tl, {'Figure 4E-H - candidate-cell HD tuning properties', ...
        sprintf(['Bars: mean of fish medians +/- SEM; large points: fish medians; ' ...
        'raw planned fish-level KW p-values; subset = %s'], ...
        strrep(cfg.MetricCellSubset,'_',' '))});
else
    title(tl, {'True tuning reliability (additional analysis)', ...
        'Bars: mean of fish medians +/- SEM; alternating 60-s blocks'});
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

    if strcmp(mode,'poster')
        statText = plannedStatsText(statsTable, defs.Variable(d));
    else
        statText = 'descriptive; included to separate reliability from SNR';
    end
    title(ax,char(defs.Title(d)));
    text(ax,0.02,0.98,statText,'Units','normalized', ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',8,'Interpreter','none','BackgroundColor','w', ...
        'Margin',2,'Clipping','off');
end

if strcmp(mode,'poster')
    name = 'Figure4E_H_candidate_cell_metrics';
else
    name = 'True_split_half_reliability';
end
pathOut = saveFigureBoth(fig, cfg.OutputDir, name, cfg.SaveFigFiles);
end

function pathOut = plotClassMetricPanels(cellTable, fishSummary, statsTable, defs, cfg, morphColors, highlyTunedOnly)
[subsetTitle,fileSuffix] = figureSubsetMetadata(highlyTunedOnly);
if ~highlyTunedOnly
    subsetTitle = [strrep(cfg.MetricCellSubset,'_',' ') ' metric subset'];
end
% Rows are CW, CCW and Symmetric; columns are the four poster metrics.
fig = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Candidate metrics by functional class', ...
    'Position',[40 40 1900 1250]);
tl = tiledlayout(fig, numel(cfg.ClassOrder), height(defs), ...
    'TileSpacing','compact','Padding','compact');
title(tl, {['Candidate-cell HD properties by functional class - ' subsetTitle], ...
    'Bars = mean of fish medians +/- SEM; brackets show raw planned fish-level KW tests'});

for c = 1:numel(cfg.ClassOrder)
    className = string(cfg.ClassOrder{c});
    for d = 1:height(defs)
        tileNumber = (c-1).*height(defs) + d;
        ax = nexttile(tl, tileNumber);
        hold(ax,'on');
        variable = char(defs.Variable(d));

        for m = 1:numel(cfg.MorphOrder)
            morph = string(cfg.MorphOrder{m});
            useFish = fishSummary.Class == className & ...
                fishSummary.Metric == defs.Variable(d) & ...
                strcmpi(fishSummary.Morph,morph) & isfinite(fishSummary.Value);
            yFish = double(fishSummary.Value(useFish));
            drawFishBar(ax, m, yFish, morphColors(m,:), defs.LogScale(d));

            useCell = cellTable.ClassLabel == className & ...
                strcmpi(cellTable.Morph,morph) & isfinite(cellTable.(variable));
            yCell = double(cellTable.(variable)(useCell));
            jitter = (rand(size(yCell))-0.5).*0.46;
            try
                scatter(ax, m+jitter, yCell, 8, morphColors(m,:), 'filled', ...
                    'MarkerFaceAlpha',0.10,'MarkerEdgeAlpha',0.04);
            catch
                scatter(ax, m+jitter, yCell, 8, morphColors(m,:), 'filled');
            end

            if ~isempty(yFish)
                xFish = m + linspace(-0.12,0.12,numel(yFish))';
                scatter(ax, xFish, yFish, 38, morphColors(m,:), 'filled', ...
                    'MarkerEdgeColor','k','LineWidth',0.7);
            end
        end

        xlim(ax,[0.45 numel(cfg.MorphOrder)+0.55]);
        set(ax,'XTick',1:numel(cfg.MorphOrder), ...
            'XTickLabel',cfg.MorphOrder);
        if c < numel(cfg.ClassOrder)
            set(ax,'XTickLabel',[]);
        end
        ylabel(ax,defs.YLabel(d));
        box(ax,'off');
        grid(ax,'on');
        ax.GridAlpha = 0.12;
        if defs.LogScale(d)
            set(ax,'YScale','log');
        end

        addPlannedClassStatBars(ax,statsTable,className, ...
            defs.Variable(d),cfg.MorphOrder,defs.LogScale(d));
        title(ax,sprintf('%s - %s',className,defs.Title(d)));
    end
end

pathOut = saveFigureBoth(fig, cfg.OutputDir, ...
    ['Figure4E_H_metrics_by_CW_CCW_Symmetric_class_3x4' fileSuffix], ...
    cfg.SaveFigFiles);
end

function addPlannedClassStatBars(ax,statsTable,className,metric,morphOrder,isLogScale)
rows = find(statsTable.Class == className & ...
    statsTable.Metric == metric & ...
    statsTable.Test == "planned_two_group_KW_unadjusted");
if isempty(rows)
    return;
end

limits = ylim(ax);
if isLogScale
    dataTop = limits(2);
    if ~isfinite(dataTop) || dataTop <= 0
        return;
    end
    levelFactor = 1.32;
    for i = 1:numel(rows)
        row = rows(i);
        x1 = find(strcmpi(string(morphOrder),statsTable.Group1(row)),1);
        x2 = find(strcmpi(string(morphOrder),statsTable.Group2(row)),1);
        if isempty(x1) || isempty(x2); continue; end
        y = dataTop .* levelFactor.^i;
        yTick = y ./ 1.045;
        plot(ax,[x1 x1 x2 x2],[yTick y y yTick],'k-', ...
            'LineWidth',1.15,'Clipping','off');
        text(ax,mean([x1 x2]),y.*1.025, ...
            plannedSignificanceLabel(statsTable.PValue(row)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontWeight','bold');
    end
    ylim(ax,[limits(1) dataTop.*levelFactor.^(numel(rows)+0.55)]);
else
    dataRange = limits(2)-limits(1);
    if ~isfinite(dataRange) || dataRange <= 0
        dataRange = max(1,abs(limits(2)));
    end
    base = limits(2) + 0.08.*dataRange;
    step = 0.13.*dataRange;
    tick = 0.035.*dataRange;
    for i = 1:numel(rows)
        row = rows(i);
        x1 = find(strcmpi(string(morphOrder),statsTable.Group1(row)),1);
        x2 = find(strcmpi(string(morphOrder),statsTable.Group2(row)),1);
        if isempty(x1) || isempty(x2); continue; end
        y = base + (i-1).*step;
        plot(ax,[x1 x1 x2 x2],[y-tick y y y-tick],'k-', ...
            'LineWidth',1.15,'Clipping','off');
        text(ax,mean([x1 x2]),y+0.02.*dataRange, ...
            plannedSignificanceLabel(statsTable.PValue(row)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontWeight','bold');
    end
    ylim(ax,[limits(1) base + max(numel(rows)-1,0).*step + 0.16.*dataRange]);
end
addPlannedClassPValueText(ax,statsTable(rows,:));
end

function addPlannedClassPValueText(ax,rows)
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

function label = plannedSignificanceLabel(p)
if ~isfinite(p)
    label = 'n/a';
elseif p < 0.001
    label = '***';
elseif p < 0.01
    label = '**';
elseif p < 0.05
    label = '*';
else
    label = 'ns';
end
end

function textOut = plannedClassStatsText(statsTable, className, metric)
omnibusRow = find(statsTable.Class == className & ...
    statsTable.Metric == metric & ...
    statsTable.Test == "omnibus_three_group_KW", 1);
pairRows = find(statsTable.Class == className & ...
    statsTable.Metric == metric & ...
    statsTable.Test == "planned_two_group_KW_unadjusted");
if isempty(omnibusRow) && isempty(pairRows)
    textOut = 'fish-level Kruskal-Wallis tests unavailable';
    return;
end

pieces = strings(1 + numel(pairRows),1);
if isempty(omnibusRow)
    pieces(1) = "KW all p=NA";
else
    pOmnibus = statsTable.PValue(omnibusRow);
    pieces(1) = "KW all p=" + string(formatP(pOmnibus)) + ...
        significanceStars(pOmnibus);
end
for i = 1:numel(pairRows)
    r = pairRows(i);
    g1 = abbreviateMorph(statsTable.Group1(r));
    g2 = abbreviateMorph(statsTable.Group2(r));
    pValue = statsTable.PValue(r);
    pieces(i+1) = g1 + "-" + g2 + " p=" + string(formatP(pValue)) + ...
        significanceStars(pValue);
end
textOut = strjoin(cellstr(pieces), '; ');
end

function textOut = plannedStatsText(statsTable, metric)
rows = find(statsTable.Metric == metric & ...
    statsTable.Test == "planned_two_group_KW");
if isempty(rows)
    textOut = 'planned fish-level rank comparisons unavailable';
    return;
end
pieces = strings(numel(rows),1);
for i = 1:numel(rows)
    r = rows(i);
    g1 = abbreviateMorph(statsTable.Group1(r));
    g2 = abbreviateMorph(statsTable.Group2(r));
    p = statsTable.PValue(r);
    pieces(i) = g1 + "-" + g2 + " p=" + string(formatP(p)) + ...
        significanceStars(p);
end
textOut = strjoin(cellstr(pieces), '; ');
end

function a = abbreviateMorph(morph)
morph = lower(string(morph));
if morph == "surface"
    a = "S";
elseif morph == "pachon"
    a = "P";
elseif morph == "molino"
    a = "M";
else
    a = extractBefore(string(morph),2);
end
end

function stars = significanceStars(p)
if ~isfinite(p)
    stars = "";
elseif p < 0.001
    stars = "***";
elseif p < 0.01
    stars = "**";
elseif p < 0.05
    stars = "*";
else
    stars = " n.s.";
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

pathOut = saveFigureBoth(fig, cfg.OutputDir, 'Phase_alignment_QC', cfg.SaveFigFiles);
end

function pathOut = saveFigureBoth(fig, outputDir, baseName, ~)
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

function s = formatP(p)
if ~isfinite(p)
    s = 'NA';
elseif p < 0.001
    s = '<0.001';
else
    s = sprintf('%.3f',p);
end
end

function angle = wrapToPiLocal(angle)
angle = mod(angle + pi, 2*pi) - pi;
end
