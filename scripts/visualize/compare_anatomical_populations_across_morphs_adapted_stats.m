function results = compare_anatomical_populations_across_morphs_adapted_stats(cfg)
%COMPARE_ANATOMICAL_POPULATIONS_ACROSS_MORPHS_ADAPTED_STATS
% Anatomical-only comparison of CCW, CW and Symmetric neuron populations
% across Surface, Pachon and Molino Astyanax morphs.
%
% IMPORTANT:
%   - Functional model coefficients, AHV tuning, phase tuning, PCA/rPC
%     geometry and behavior are NOT compared in this script.
%   - HD_AHV_model_results.mat is used only to obtain:
%         1) recording/session identity,
%         2) neuron identity,
%         3) the final CCW / CW / Symmetric class label.
%   - Anatomical coordinates come directly from xCenteredUm/yCenteredUm
%     stored in each morph's HD_AHV_model_results.mat. These are the exact
%     original-ROI, per-recording FOV-corrected coordinates used by the
%     morph-specific anatomical plots.
%
% Main fish-level anatomical metrics, computed separately for each class:
%   1) Number of cells.
%   2) Percentage of all classified cells in that fish.
%   3) Spatial density within the sampled anatomical region.
%   4) Local nearest-neighbour density.
%   5) Compactness: hull area per cell.
%   6) Position along normalized x.
%   8) Central x-region enrichment and fixed-bandwidth KDE peak separation.
% Additional overall plots compare all cells with the pooled highly tuned
% population (CCW + CW + Symmetric) for cell count and compactness.
%
% Statistics (simple distribution-aware version):
%   - The statistical unit is the fish/recording; neurons are never treated
%     as independent biological replicates.
%   - Counts are compared with a quasi-Poisson GLM and an offset equal to the
%     number of all imaged neurons (read from the original *_RASTER.mat).
%   - Proportions are compared with a quasi-binomial GLM using the exact
%     numerator and denominator for every fish.
%   - Spatial summaries are compared with a fish-level weighted permutation
%     test; precision weights depend on the number of contributing neurons.
%   - KDE peak separation uses the same permutation test, but its precision
%     comes from the existing within-fish bootstrap confidence interval.
%   - Morph-pair p-values are reported without multiplicity adjustment.
%
% Required per morph:
%   - one HD_AHV_model_results.mat containing xCenteredUm, yCenteredUm,
%     session, neuronID and class for every modeled neuron.
%
% The script saves figures only. Numeric results are returned in memory.

if nargin < 1 || isempty(cfg)
    cfg = pipeline_config();
end
assert(isstruct(cfg) && isfield(cfg, 'DataRoot'), ...
    'Pass the struct returned by pipeline_config.');

%% ============================ USER PARAMETERS ============================
P = struct();

% The defaults below use the three morphs from the previous script.
% Change only these paths if your folders move.
P.morphs = struct([]);

P.morphs(1).name      = 'surface';
P.morphs(1).display   = 'Surface';
P.morphs(1).rootDir   = fullfile(cfg.DataRoot, 'surface');
P.morphs(1).modelFile = fullfile(P.morphs(1).rootDir, ...
    'HD_AHV_neuron_model_outputs_candidate_clean', 'HD_AHV_model_results.mat');
P.morphs(1).swapXY    = false;
P.morphs(1).flipX     = false;
P.morphs(1).flipY     = false;

P.morphs(2).name      = 'pachon';
P.morphs(2).display   = 'Pachon';
P.morphs(2).rootDir   = fullfile(cfg.DataRoot, 'pachon');
P.morphs(2).modelFile = fullfile(P.morphs(2).rootDir, ...
    'HD_AHV_neuron_model_outputs_candidate_clean', 'HD_AHV_model_results.mat');
P.morphs(2).swapXY    = false;
P.morphs(2).flipX     = false;
P.morphs(2).flipY     = false;

P.morphs(3).name      = 'molino';
P.morphs(3).display   = 'Molino';
P.morphs(3).rootDir   = fullfile(cfg.DataRoot, 'molino');
P.morphs(3).modelFile = fullfile(P.morphs(3).rootDir, ...
    'HD_AHV_neuron_model_outputs_candidate_clean', 'HD_AHV_model_results.mat');
P.morphs(3).swapXY    = false;
P.morphs(3).flipX     = false;
P.morphs(3).flipY     = false;

P.outputDir = cfg.AnatomyComparisonDir;
P.figureDir = fullfile(P.outputDir, 'figures');

% Original RASTER files used to count every imaged neuron in each fish.
% Generated population-specific RASTER files are explicitly excluded.
P.totalRaster.pattern = '*_RASTER.mat';
P.totalRaster.excludeTokens = { ...
    'SELECTED_ROIS_RASTER.mat', ...
    'NON_SELECTED_ROIS_RASTER.mat', ...
    'ANATOMICAL_ZONE_ROIS_RASTER.mat', ...
    'NOT_IN_ANATOMICAL_ZONE_ROIS_RASTER.mat'};
P.totalRaster.neuronVariable = 'raster';
P.totalRaster.neuronDimension = 2; % project convention: time x neurons
P.totalRaster.failOnMissing = false;

% Legacy selected-file settings are retained only because older helper
% functions remain at the end of the script; the main analysis does not use them.
P.selectedAllCellsPattern = '*SELECTED_ROIS_ALL_CELLS*.mat';
P.excludeTokens = {'NON_SELECTED_ROIS', 'NON_SELECTED'};

% Population order requested for all outputs.
P.classes = {'CCW', 'CW', 'Symmetric'};
P.classColors = [ ...
    0.90 0.20 0.20; ... % CCW
    0.10 0.65 0.90; ... % CW
    0.55 0.55 0.55];    % Symmetric

% Morph colors used in boxplots and density curves.
P.morphColors = [ ...
    0.20 0.70 0.30; ... % Surface
    0.95 0.85 0.20; ... % Pachon
    0.85 0.18 0.22];    % Molino

% Anatomy and matching.
% 'row_order' is useful when the model table was created in selected-ROI
% order. Every fallback is recorded in the returned anatomyMatchReport so it
% can be checked. Set to 'none' to forbid this fallback.
P.anatomy.matchingFallback = 'none';  % 'row_order' or 'none'
P.anatomy.minCellsForSpread = 2;
P.anatomy.minCellsForHull = 3;
P.anatomy.minCellsForNearestNeighbour = 2;
P.anatomy.axisEdges = linspace(-1.05, 1.05, 43);
P.anatomy.densitySmoothingSigmaBins = 1.0;

% High-resolution x-density visualization in physical micrometre units.
% A standard logarithmic x-axis cannot represent negative positions or zero,
% so only a signed-log display is saved; it expands the region around x = 0.
P.anatomy.xDensity.fullBinWidthUm = 5;
P.anatomy.xDensity.smoothingSigmaUm = 3;
P.anatomy.xDensity.signedLogScaleUm = 10;

% Fish-level central x-axis metrics in physical micrometre coordinates.
% Fractions are normalized by the total number of neurons in the relevant
% class in that fish, not by the number of neurons inside the central window.
P.anatomy.centralX.ccwRangeUm = [-20 0];        % -20 <= x < 0
P.anatomy.centralX.cwRangeUm = [0 20];          % 0 < x <= 20
P.anatomy.centralX.symmetricRangeUm = [-10 10];% -10 <= x <= 10
P.anatomy.centralX.minCellsForFraction = 5;

% Directional peak positions are estimated separately for CW and CCW with
% a Gaussian KDE using the same fixed bandwidth for every fish and morph.
% Peaks are searched only inside the central range below, so distant minor
% modes cannot replace the central anatomical peak of interest.
P.anatomy.centralX.kdeBandwidthUm = 8;
P.anatomy.centralX.kdeSearchRangeUm = [-50 50];
P.anatomy.centralX.kdeGridStepUm = 0.25;
P.anatomy.centralX.minCellsForPeak = 5;

% Neuron-resampling bootstrap uncertainty for each fish's peak positions.
P.anatomy.centralX.bootstrapIterations = 1000;
P.anatomy.centralX.bootstrapCI = 95;
P.anatomy.centralX.bootstrapRandomSeed = 1701;

% CW-versus-CCW left-right organization.
% A fish is included only when both CW and CCW contain at least this many
% anatomically matched cells.
P.anatomy.minCellsPerDirectionForLROpposition = 5;

% Independent anatomy-matching audit.
% The main analysis uses the model-table coordinates directly. This audit
% independently reloads the ORIGINAL ALL_CELLS file saved in each Result,
% extracts each ROI using neuronID, reconstructs the FOV-corrected position,
% and compares it with x/y and xCenteredUm/yCenteredUm in the model table.
P.matchAudit.enabled = true;
P.matchAudit.minNeuronsForCorrelation = 5;
P.matchAudit.minValidOriginalRoiFraction = 0.99;
P.matchAudit.goodCorrelation = 0.999;
P.matchAudit.maxMedianPixelErrorPx = 1e-6;
P.matchAudit.maxMedianFovErrorUm = 1e-6;
P.matchAudit.addSessionLabels = false;
P.matchAudit.makeMorphOverlayPanels = true;

% Statistical inference. These defaults keep the analysis deliberately
% simple. Increase permutationIterations for the final frozen analysis if
% desired; 20,000 already gives stable pairwise p-values for most datasets.
P.stats.alpha = 0.05;
P.stats.minFishPerMorph = 2;
P.stats.permutationIterations = 20000;
P.stats.maxExactPermutations = 50000;
P.stats.maxRelativeWeight = 5;
P.stats.peakSEFloorUm = P.anatomy.centralX.kdeGridStepUm / 2;
P.stats.permutationRandomSeed = 4317;
P.stats.showNonsignificant = true;
P.stats.labelMode = 'stars'; % 'stars' or 'exact'

% Figures.
P.figureVisible = cfg.FigureVisible;
P.figureFormat = 'png';  % Primary format; SVG is also always saved
P.figureResolution = 300;
P.addFishLabels = false;
P.randomSeed = 7;

if ~exist(P.outputDir, 'dir'); mkdir(P.outputDir); end
if ~exist(P.figureDir, 'dir'); mkdir(P.figureDir); end

%% ========================== LOAD INPUT TABLES ============================
allNeuronTables = cell(numel(P.morphs), 1);
allSourceTables = cell(numel(P.morphs), 1);

for m = 1:numel(P.morphs)
    morph = P.morphs(m);
    fprintf('\n==================== %s ====================\n', morph.display);

    if ~isfile(morph.modelFile)
        warning('Missing model file for %s:\n  %s', morph.display, morph.modelFile);
        allNeuronTables{m} = table();
        allSourceTables{m} = table();
    else
        T = loadNeuronTableFromModelFile(morph.modelFile);
        T = standardizeNeuronTable(T, morph.name);
        allNeuronTables{m} = T;
        allSourceTables{m} = buildModelSourceIndex(morph.modelFile, morph.name);
        fprintf('Loaded %d model rows from:\n  %s\n', height(T), morph.modelFile);
    end
end

AllNeurons = vertcatTables(allNeuronTables);
AllNeuronsFull = AllNeurons;  % retain all model rows for exact centering audit
ModelSourceIndex = vertcatTables(allSourceTables);

if isempty(AllNeurons) || height(AllNeurons) == 0
    error('No neuron classification table was loaded. Check P.morphs(i).modelFile.');
end
% Keep only the three populations of interest. No functional coefficient is
% used below.
AllNeurons = AllNeurons(ismember(string(AllNeurons.class3), string(P.classes)), :);
if height(AllNeurons) == 0
    error('No CCW, CW or Symmetric rows were found in the model files.');
end

%% ======================= BUILD ANATOMICAL TABLES =========================
[AnatomyNeurons, AnatomyMatchReport] = ...
    buildAnatomyNeuronTableFromModelCoordinates(AllNeurons, P);
if isempty(AnatomyNeurons) || height(AnatomyNeurons) == 0
    error('No neurons with usable model-table anatomical coordinates were found.');
end

% Count all imaged neurons from the unique original *_RASTER.mat belonging
% to each fish. Population-specific/generated RASTER files are ignored.
FishImagingTotals = buildFishImagingTotals( ...
    AnatomyNeurons, ModelSourceIndex, P);

FishClassMetrics = buildFishClassMetrics(AnatomyNeurons, P);
FishClassMetrics = attachImagingTotals(FishClassMetrics, FishImagingTotals);
MetricDefinitions = getMetricDefinitions();
PairwiseStats = runAdaptiveClassStatistics( ...
    FishClassMetrics, MetricDefinitions, P);

% Overall (not class-separated) comparisons. "Highly tuned" is the pooled
% CCW + CW + Symmetric population; "All cells" comes from each fish's
% original RASTER count and original ALL_CELLS ROI coordinates.
FishOverallMetrics = buildFishOverallMetrics( ...
    AnatomyNeurons, FishImagingTotals, ModelSourceIndex, P);
OverallPairwiseStats = runOverallPopulationStatistics( ...
    FishOverallMetrics, P);

% Fish-level CW-versus-CCW left-right organization. These metrics compare the
% two directional populations within the same fish, then compare the resulting
% fish-level scores across morphs.
FishOpponentMetrics = buildFishOpponentMetrics(AnatomyNeurons, P);
OpponentMetricDefinitions = getOpponentMetricDefinitions();
OpponentPairwiseStats = runAdaptiveOpponentStatistics( ...
    FishOpponentMetrics, OpponentMetricDefinitions, P);

% Central left-right anatomical enrichment and directional KDE peak metrics.
% Every row is one fish/recording.
FishCentralXMetrics = buildFishCentralXMetrics(AnatomyNeurons, P);
CentralXMetricDefinitions = getCentralXMetricDefinitions();
CentralXPairwiseStats = runAdaptiveCentralXStatistics( ...
    FishCentralXMetrics, CentralXMetricDefinitions, P);

% Independently verify that neuronID points to the same original ROI whose
% coordinates were saved in the model table. This audit does not feed the
% main analysis; it is a QC check only.
[MatchingNeuronDiagnostics, MatchingFishDiagnostics] = ...
    buildIndependentAnatomyMatchingAudit(AllNeuronsFull, ModelSourceIndex, P);

%% ================================ PLOTS =================================
plotSpatialOverlays(AnatomyNeurons, P);
plotXAxisDensityCurvesMicrons(AnatomyNeurons, P);
plotOverallNumberOfCells(FishOverallMetrics, OverallPairwiseStats, P);
plotOverallPopulationCompactness( ...
    FishOverallMetrics, OverallPairwiseStats, P);

for i = 1:height(MetricDefinitions)
    plotMetricByClassWithPairwiseBars(FishClassMetrics, PairwiseStats, ...
        MetricDefinitions(i,:), P);
end

fprintf('\nDone. Anatomical-only figures were saved to:\n  %s\n', P.figureDir);

%% =======================================================================
results = struct('outputDir', P.outputDir, ...
    'allNeurons', AllNeurons, 'anatomyNeurons', AnatomyNeurons, ...
    'anatomyMatchReport', AnatomyMatchReport, ...
    'fishImagingTotals', FishImagingTotals, ...
    'fishClassMetrics', FishClassMetrics, ...
    'pairwiseStatistics', PairwiseStats, ...
    'fishOverallMetrics', FishOverallMetrics, ...
    'overallPairwiseStatistics', OverallPairwiseStats, ...
    'fishOpponentMetrics', FishOpponentMetrics, ...
    'opponentPairwiseStatistics', OpponentPairwiseStats, ...
    'fishCentralXMetrics', FishCentralXMetrics, ...
    'centralXPairwiseStatistics', CentralXPairwiseStats, ...
    'matchingNeuronDiagnostics', MatchingNeuronDiagnostics, ...
    'matchingFishDiagnostics', MatchingFishDiagnostics);
end

%% ============================ LOCAL FUNCTIONS ===========================
%% =======================================================================

function T = loadNeuronTableFromModelFile(modelFile)
    S = load(modelFile);

    if isfield(S, 'AllNeurons') && istable(S.AllNeurons)
        T = S.AllNeurons;
        return;
    end

    if isfield(S, 'Results')
        pieces = {};
        R = S.Results;
        for i = 1:numel(R)
            if iscell(R)
                Ri = R{i};
            else
                Ri = R(i);
            end
            if isstruct(Ri) && isfield(Ri, 'fitTable') && ...
                    istable(Ri.fitTable) && height(Ri.fitTable) > 0
                pieces{end+1,1} = Ri.fitTable; %#ok<AGROW>
            end
        end
        if ~isempty(pieces)
            T = vertcatTables(pieces);
            return;
        end
    end

    candidates = findTablesRecursive(S);
    bestScore = -Inf;
    bestTable = table();
    for i = 1:numel(candidates)
        Ti = candidates{i};
        names = lower(string(Ti.Properties.VariableNames));
        score = 0;
        score = score + 10 * any(ismember(names, ...
            {'class','classlabel','neuronclass','modelclass','functionalclass','assignedclass','finalclass'}));
        score = score + 4 * any(ismember(names, ...
            {'session','sessionname','fish','fishid','recording','recordingname'}));
        score = score + height(Ti) / 1e6;
        if score > bestScore
            bestScore = score;
            bestTable = Ti;
        end
    end

    if bestScore > 0
        T = bestTable;
    else
        error('Could not locate a usable neuron table inside %s.', modelFile);
    end
end

function list = findTablesRecursive(x)
    list = {};
    if istable(x)
        list = {x};
        return;
    end
    if isstruct(x)
        fields = fieldnames(x);
        for i = 1:numel(x)
            for f = 1:numel(fields)
                try
                    sub = findTablesRecursive(x(i).(fields{f}));
                    list = [list; sub(:)]; %#ok<AGROW>
                catch
                end
            end
        end
    elseif iscell(x)
        for i = 1:numel(x)
            try
                sub = findTablesRecursive(x{i});
                list = [list; sub(:)]; %#ok<AGROW>
            catch
            end
        end
    end
end

function T = standardizeNeuronTable(T, morphName)
    if isempty(T) || height(T) == 0
        T = table();
        return;
    end

    T.morph = repmat(string(morphName), height(T), 1);

    sessionVar = findFirstVariable(T, ...
        {'session','sessionName','fish','fishID','recording','recordingName','name'});
    if isempty(sessionVar)
        error('No session/fish/recording column was found in the model table.');
    end
    T.session = string(T.(sessionVar));
    T.sessionKey = makeSessionKey(T.session);

    neuronVar = findFirstVariable(T, ...
        {'neuronID','neuronId','roiID','roiId','cellID','cellId','selectedRoiIdx','selectedROI'});
    if isempty(neuronVar)
        % This remains useful for row-order matching.
        T.neuronID = (1:height(T))';
    else
        T.neuronID = convertToNumericColumn(T.(neuronVar), height(T));
    end

    classVar = findFirstVariable(T, ...
        {'class','classLabel','neuronClass','modelClass','functionalClass','assignedClass','finalClass'});
    if isempty(classVar)
        error('No neuron class column was found in the model table.');
    end
    T.classRaw = string(T.(classVar));
    T.class3 = canonicalClass(T.classRaw);
end

function varName = findFirstVariable(T, candidates)
    varName = '';
    vars = string(T.Properties.VariableNames);
    lowerVars = lower(vars);
    for i = 1:numel(candidates)
        idx = find(lowerVars == lower(string(candidates{i})), 1, 'first');
        if ~isempty(idx)
            varName = char(vars(idx));
            return;
        end
    end
end

function x = convertToNumericColumn(v, expectedN)
    if isnumeric(v) || islogical(v)
        x = double(v(:));
    else
        x = str2double(string(v(:)));
    end
    if numel(x) ~= expectedN
        x = nan(expectedN, 1);
    end
end

function out = canonicalClass(x)
    x = lower(strtrim(string(x)));
    out = strings(size(x));
    out(:) = "Other";

    % Check CCW first because the text 'ccw' contains 'cw'.
    out(contains(x, 'ccw') | contains(x, 'counter')) = "CCW";
    out((contains(x, 'cw') | contains(x, 'clock')) & out == "Other") = "CW";
    out(contains(x, 'sym') | strcmp(x, 'hd') | ...
        contains(x, 'hd-like') | contains(x, 'hd_like')) = "Symmetric";
end

function F = discoverSelectedAllCellsFiles(rootDir, morphName, P)
    F = table();
    if ~isfolder(rootDir)
        warning('Root folder does not exist: %s', rootDir);
        return;
    end

    d = dir(fullfile(rootDir, '**', P.selectedAllCellsPattern));
    d = d(~[d.isdir]);

    keep = true(numel(d), 1);
    for i = 1:numel(P.excludeTokens)
        keep = keep & ~contains(lower({d.name})', lower(P.excludeTokens{i}));
    end
    d = d(keep);
    if isempty(d); return; end

    fullPaths = string(fullfile({d.folder}, {d.name}))';
    [~, order] = sort(lower(fullPaths));
    d = d(order);

    rows = cell(numel(d), 1);
    for i = 1:numel(d)
        [~, recFolder] = fileparts(d(i).folder);
        rows{i} = table(string(morphName), string(recFolder), ...
            makeSessionKey(string(recFolder)), ...
            string(fullfile(d(i).folder, d(i).name)), ...
            'VariableNames', {'morph','session','sessionKey','allCellsFile'});
    end
    F = vertcatTables(rows);
end


function [AnatomyNeurons, MatchReport] = ...
        buildAnatomyNeuronTableFromModelCoordinates(AllNeurons, P)
    % Use the exact coordinates already generated by the morph-specific
    % HD/AHV script. This avoids any rematching to SELECTED_ROIS files.
    required = {'morph','session','sessionKey','class3','neuronID', ...
        'xCenteredUm','yCenteredUm'};
    missing = required(~ismember(required, AllNeurons.Properties.VariableNames));
    if ~isempty(missing)
        error(['The model table is missing required variable(s): %s. ' ...
            'Rerun the morph-specific model script that saves ' ...
            'xCenteredUm/yCenteredUm.'], strjoin(missing, ', '));
    end

    n = height(AllNeurons);
    AnatomyNeurons = table();
    AnatomyNeurons.morph = string(AllNeurons.morph);
    AnatomyNeurons.session = string(AllNeurons.session);
    AnatomyNeurons.sessionKey = string(AllNeurons.sessionKey);
    AnatomyNeurons.class3 = string(AllNeurons.class3);
    AnatomyNeurons.neuronID = double(AllNeurons.neuronID);
    AnatomyNeurons.xUm = double(AllNeurons.xCenteredUm);
    AnatomyNeurons.yUm = double(AllNeurons.yCenteredUm);

    % Retain the raw pre-FOV coordinate when available for auditing.
    if ismember('xCenteredRawUm', AllNeurons.Properties.VariableNames)
        AnatomyNeurons.xRawUm = double(AllNeurons.xCenteredRawUm);
    else
        AnatomyNeurons.xRawUm = nan(n,1);
    end
    if ismember('yCenteredRawUm', AllNeurons.Properties.VariableNames)
        AnatomyNeurons.yRawUm = double(AllNeurons.yCenteredRawUm);
    else
        AnatomyNeurons.yRawUm = nan(n,1);
    end

    % Isotropic within-fish normalization. The same positive scale is used
    % for x and y, so left/right signs and 2-D geometry are preserved.
    AnatomyNeurons.xNorm = nan(n,1);
    AnatomyNeurons.yNorm = nan(n,1);
    AnatomyNeurons.normalizationScaleUm = nan(n,1);

    morphOrder = string({P.morphs.name});
    reportRows = {};
    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        keys = unique(AnatomyNeurons.sessionKey( ...
            AnatomyNeurons.morph == morphName), 'stable');

        for k = 1:numel(keys)
            idx = AnatomyNeurons.morph == morphName & ...
                AnatomyNeurons.sessionKey == keys(k);

            x = AnatomyNeurons.xUm(idx);
            y = AnatomyNeurons.yUm(idx);
            finiteXY = isfinite(x) & isfinite(y);

            scale = max([abs(x(finiteXY)); abs(y(finiteXY))], [], 'omitnan');
            if isempty(scale) || ~isfinite(scale) || scale <= 0
                scale = NaN;
            else
                AnatomyNeurons.xNorm(idx) = x ./ scale;
                AnatomyNeurons.yNorm(idx) = y ./ scale;
                AnatomyNeurons.normalizationScaleUm(idx) = scale;
            end

            nRows = sum(idx);
            nUsable = sum(finiteXY);
            sessionName = AnatomyNeurons.session(find(idx,1,'first'));
            reportRows{end+1,1} = table( ...
                morphName, sessionName, keys(k), "", ...
                nRows, nRows, nUsable, ...
                "model_table_exact_coordinates", "not_applicable", ...
                'VariableNames', {'morph','session','sessionKey','allCellsFile', ...
                'nModelRows','nAnatomyRois','nMatched','roiMatchMode', ...
                'fileMatchMode'}); %#ok<AGROW>
        end
    end

    MatchReport = vertcatTables(reportRows);

    % Drop only neurons without a usable FOV-corrected position.
    good = isfinite(AnatomyNeurons.xUm) & isfinite(AnatomyNeurons.yUm) & ...
        isfinite(AnatomyNeurons.xNorm) & isfinite(AnatomyNeurons.yNorm);
    if any(~good)
        warning('Dropping %d/%d neurons with missing anatomical coordinates.', ...
            sum(~good), height(AnatomyNeurons));
        AnatomyNeurons = AnatomyNeurons(good,:);
    end
end

function [AnatomyNeurons, MatchReport] = buildAnatomyNeuronTable(AllNeurons, AllCellsFiles, P)
    neuronRows = {};
    reportRows = {};

    morphList = unique(AllNeurons.morph, 'stable');
    for m = 1:numel(morphList)
        morphName = morphList(m);
        Tmorph = AllNeurons(AllNeurons.morph == morphName, :);
        keys = unique(Tmorph.sessionKey, 'stable');
        transform = getMorphTransform(morphName, P);

        for k = 1:numel(keys)
            key = keys(k);
            Tfish = Tmorph(Tmorph.sessionKey == key, :);
            sessionName = Tfish.session(1);

            [allCellsFile, fileMatchMode] = findAllCellsForSession( ...
                AllCellsFiles, morphName, key, sessionName);

            if strlength(allCellsFile) == 0 || ~isfile(allCellsFile)
                warning('No selected ALL_CELLS file matched %s / %s.', morphName, sessionName);
                reportRows{end+1,1} = makeMatchReportRow(morphName, sessionName, key, ...
                    "", height(Tfish), NaN, 0, "missing_all_cells", fileMatchMode); %#ok<AGROW>
                continue;
            end

            try
                A = extractAnatomyFromAllCellsFile(char(allCellsFile));
                A = applyCoordinateTransform(A, transform);
            catch ME
                warning('Could not extract anatomy from %s: %s', allCellsFile, ME.message);
                reportRows{end+1,1} = makeMatchReportRow(morphName, sessionName, key, ...
                    allCellsFile, height(Tfish), NaN, 0, "read_error", fileMatchMode); %#ok<AGROW>
                continue;
            end

            [roiIdx, roiMatchMode] = matchModelRowsToAnatomyRows(Tfish, A, P);
            good = isfinite(roiIdx) & roiIdx >= 1 & roiIdx <= height(A);
            nMatched = sum(good);

            reportRows{end+1,1} = makeMatchReportRow(morphName, sessionName, key, ...
                allCellsFile, height(Tfish), height(A), nMatched, roiMatchMode, fileMatchMode); %#ok<AGROW>

            if nMatched == 0
                warning('Zero model rows matched anatomy for %s / %s.', morphName, sessionName);
                continue;
            end

            Tused = Tfish(good, :);
            Aused = A(roiIdx(good), :);

            R = table();
            R.morph = repmat(morphName, nMatched, 1);
            R.session = repmat(sessionName, nMatched, 1);
            R.sessionKey = repmat(key, nMatched, 1);
            R.class3 = Tused.class3;
            R.neuronID = Tused.neuronID;
            R.anatomyRoiIdx = Aused.roiIdx;
            R.originalRoiIdx = Aused.originalRoiIdx;
            R.xPixels = Aused.x;
            R.yPixels = Aused.y;
            R.xNormRawImage = Aused.xNormRaw;
            R.yNormRawImage = Aused.yNormRaw;
            R.xNorm = Aused.xNorm;
            R.yNorm = Aused.yNorm;
            R.imageWidth = Aused.imageWidth;
            R.imageHeight = Aused.imageHeight;

            % Reference coordinates already computed by the morph-specific
            % HD/AHV script from the ORIGINAL ALL_CELLS file.
            R.modelXRawUm = numericTableColumnOrNaN(Tused, ...
                'xCenteredRawUm', nMatched);
            R.modelYRawUm = numericTableColumnOrNaN(Tused, ...
                'yCenteredRawUm', nMatched);
            R.modelXFovUm = numericTableColumnOrNaN(Tused, ...
                'xCenteredUm', nMatched);
            R.modelYFovUm = numericTableColumnOrNaN(Tused, ...
                'yCenteredUm', nMatched);
            R.modelFovAngleDeg = numericTableColumnOrNaN(Tused, ...
                'fov_angle_deg', nMatched);
            R.modelPixelLengthX = numericTableColumnOrDefault(Tused, ...
                'pixelLengthX', nMatched, 1);
            R.modelPixelLengthY = numericTableColumnOrDefault(Tused, ...
                'pixelLengthY', nMatched, 1);

            % Rebuild the exact two coordinate stages used by the
            % morph-specific script, but now from the SELECTED ALL_CELLS
            % centroids matched by this cross-morph script.
            xCenterMatchedPx = median(Aused.x, 'omitnan');
            yCenterMatchedPx = median(Aused.y, 'omitnan');
            R.crossCandidateCenteredRawXUm = ...
                (Aused.x - xCenterMatchedPx) .* R.modelPixelLengthX;
            R.crossCandidateCenteredRawYUm = ...
                -(Aused.y - yCenterMatchedPx) .* R.modelPixelLengthY;

            fovAngleThisFish = median(R.modelFovAngleDeg, 'omitnan');
            if isfinite(fovAngleThisFish)
                thetaDiag = deg2rad(fovAngleThisFish);
                Rdiag = [cos(thetaDiag) sin(thetaDiag); ...
                        -sin(thetaDiag) cos(thetaDiag)];
                xyDiag = [ ...
                    (Aused.x - xCenterMatchedPx) .* R.modelPixelLengthX, ...
                    (Aused.y - yCenterMatchedPx) .* R.modelPixelLengthY];
                xyDiagFov = xyDiag * Rdiag';
                R.crossRecomputedFovXUm = xyDiagFov(:,1);
                R.crossRecomputedFovYUm = -xyDiagFov(:,2);
            else
                R.crossRecomputedFovXUm = nan(nMatched,1);
                R.crossRecomputedFovYUm = nan(nMatched,1);
            end
            R.allCellsFile = repmat(allCellsFile, nMatched, 1);
            R.roiMatchMode = repmat(string(roiMatchMode), nMatched, 1);
            R.fileMatchMode = repmat(string(fileMatchMode), nMatched, 1);
            neuronRows{end+1,1} = R; %#ok<AGROW>
        end
    end

    AnatomyNeurons = table();

    AnatomyNeurons.morph     = AllNeurons.morph;
    AnatomyNeurons.session   = AllNeurons.session;
    AnatomyNeurons.sessionKey = AllNeurons.sessionKey;
    AnatomyNeurons.class3    = AllNeurons.class3;
    AnatomyNeurons.neuronID  = AllNeurons.neuronID;
    
    % Exact FOV-corrected coordinates used by the morph-specific plots
    AnatomyNeurons.xUm = double(AllNeurons.xCenteredUm);
    AnatomyNeurons.yUm = double(AllNeurons.yCenteredUm);
    MatchReport = vertcatTables(reportRows);

end

function R = makeMatchReportRow(morph, session, key, file, nModel, nAnatomy, nMatched, roiMode, fileMode)
    R = table(string(morph), string(session), string(key), string(file), ...
        nModel, nAnatomy, nMatched, string(roiMode), string(fileMode), ...
        'VariableNames', {'morph','session','sessionKey','allCellsFile', ...
        'nModelRows','nAnatomyRois','nMatched','roiMatchMode','fileMatchMode'});
end

function transform = getMorphTransform(morphName, P)
    transform = struct('swapXY', false, 'flipX', false, 'flipY', false);
    idx = find(strcmpi(string({P.morphs.name}), string(morphName)), 1, 'first');
    if isempty(idx); return; end
    transform.swapXY = logical(P.morphs(idx).swapXY);
    transform.flipX = logical(P.morphs(idx).flipX);
    transform.flipY = logical(P.morphs(idx).flipY);
end

function [filePath, matchMode] = findAllCellsForSession(AllCellsFiles, morphName, key, sessionName)
    filePath = "";
    matchMode = "unmatched";
    if isempty(AllCellsFiles); return; end

    idx = find(AllCellsFiles.morph == morphName & AllCellsFiles.sessionKey == key);
    if ~isempty(idx)
        filePath = chooseFirstPath(AllCellsFiles.allCellsFile(idx));
        matchMode = "exact_session_key";
        return;
    end

    sameMorph = AllCellsFiles.morph == morphName;
    candidateKeys = string(AllCellsFiles.sessionKey);
    sessionContainsCandidate = false(height(AllCellsFiles),1);
    for i = 1:height(AllCellsFiles)
        sessionContainsCandidate(i) = contains(lower(string(sessionName)), lower(candidateKeys(i)));
    end
    idx = find(sameMorph & ...
        (sessionContainsCandidate | contains(lower(candidateKeys), lower(string(key)))));
    if ~isempty(idx)
        filePath = chooseFirstPath(AllCellsFiles.allCellsFile(idx));
        matchMode = "contains_fallback";
    end
end

function p = chooseFirstPath(paths)
    paths = sort(string(paths));
    if isempty(paths)
        p = "";
    else
        p = paths(1);
        if numel(paths) > 1
            warning('Multiple ALL_CELLS files matched one recording; using: %s', p);
        end
    end
end

function A = extractAnatomyFromAllCellsFile(allCellsFile)
    S = load(allCellsFile);

    if isfield(S, 'cell_number')
        nRois = double(S.cell_number(1));
    elseif isfield(S, 'cell_per')
        nRois = numel(S.cell_per);
    elseif isfield(S, 'cells')
        nRois = numel(S.cells);
    else
        error('Expected cell_number, cell_per or cells.');
    end

    [imageHeight, imageWidth] = inferImageSize(S);
    x = nan(nRois, 1);
    y = nan(nRois, 1);

    if isfield(S, 'cell_per') && numel(S.cell_per) >= nRois
        for i = 1:nRois
            boundary = S.cell_per{i};
            if isempty(boundary); continue; end
            boundary = double(boundary);
            if size(boundary, 2) >= 2
                x(i) = mean(boundary(:,1), 'omitnan');
                y(i) = mean(boundary(:,2), 'omitnan');
            elseif size(boundary, 1) >= 2
                x(i) = mean(boundary(1,:), 'omitnan');
                y(i) = mean(boundary(2,:), 'omitnan');
            end
        end
    elseif isfield(S, 'cells') && isfinite(imageHeight) && isfinite(imageWidth)
        for i = 1:nRois
            pix = S.cells{i};
            if isempty(pix); continue; end
            [yy, xx] = ind2sub([imageHeight imageWidth], double(pix(:)));
            x(i) = mean(xx, 'omitnan');
            y(i) = mean(yy, 'omitnan');
        end
    else
        error('Could not calculate ROI centroids.');
    end

    originalRoiIdx = nan(nRois, 1);
    originalCandidates = {'subsetOriginalRoiIdx','selectedOriginalRoiIdx', ...
        'originalRoiIdx','originalROIIdx','keptOriginalRoiIdx'};
    for i = 1:numel(originalCandidates)
        if isfield(S, originalCandidates{i})
            candidate = double(S.(originalCandidates{i})(:));
            if numel(candidate) == nRois
                originalRoiIdx = candidate;
                break;
            end
        end
    end

    if isfinite(imageWidth) && imageWidth > 1
        xNorm = (x - (imageWidth + 1)/2) ./ ((imageWidth - 1)/2);
    else
        xNorm = minMaxNormalizeToMinusOnePlusOne(x);
    end
    if isfinite(imageHeight) && imageHeight > 1
        yNorm = (y - (imageHeight + 1)/2) ./ ((imageHeight - 1)/2);
    else
        yNorm = minMaxNormalizeToMinusOnePlusOne(y);
    end

    xNormRaw = xNorm;
    yNormRaw = yNorm;
    A = table((1:nRois)', originalRoiIdx, x, y, ...
        xNormRaw, yNormRaw, xNorm, yNorm, ...
        repmat(imageWidth, nRois, 1), repmat(imageHeight, nRois, 1), ...
        'VariableNames', {'roiIdx','originalRoiIdx','x','y', ...
        'xNormRaw','yNormRaw','xNorm','yNorm', ...
        'imageWidth','imageHeight'});
end

function [h, w] = inferImageSize(S)
    h = NaN;
    w = NaN;
    candidates = {'avg','bkg','imageAvg','im','background'};
    for i = 1:numel(candidates)
        if isfield(S, candidates{i})
            value = S.(candidates{i});
            if isnumeric(value) && ndims(value) >= 2
                h = size(value, 1);
                w = size(value, 2);
                return;
            end
        end
    end
end

function z = minMaxNormalizeToMinusOnePlusOne(x)
    x = double(x(:));
    mn = min(x, [], 'omitnan');
    mx = max(x, [], 'omitnan');
    if ~isfinite(mn) || ~isfinite(mx) || mx <= mn
        z = nan(size(x));
    else
        z = 2 * (x - mn) ./ (mx - mn) - 1;
    end
end

function A = applyCoordinateTransform(A, transform)
    x = A.xNorm;
    y = A.yNorm;
    if transform.swapXY
        tmp = x;
        x = y;
        y = tmp;
    end
    if transform.flipX; x = -x; end
    if transform.flipY; y = -y; end
    A.xNorm = x;
    A.yNorm = y;
end

function [roiIdx, matchMode] = matchModelRowsToAnatomyRows(Tfish, A, P)
    nModel = height(Tfish);
    nAnatomy = height(A);
    roiIdx = nan(nModel, 1);
    matchMode = "unmatched";

    neuronID = double(Tfish.neuronID(:));

    % 1) neuronID already indexes the selected anatomy file.
    validSelectedIndex = all(isfinite(neuronID)) && ...
        all(neuronID >= 1) && all(neuronID <= nAnatomy) && ...
        numel(unique(neuronID)) == numel(neuronID);
    if validSelectedIndex
        roiIdx = neuronID;
        matchMode = "neuronID_selected_roi_index";
        return;
    end

    % 2) neuronID refers to original ROI indices saved in the selected file.
    if any(isfinite(A.originalRoiIdx)) && any(isfinite(neuronID))
        [tf, loc] = ismember(neuronID, A.originalRoiIdx);
        if sum(tf) >= max(1, round(0.75 * nModel))
            roiIdx(tf) = loc(tf);
            matchMode = "neuronID_original_roi_index";
            return;
        end
    end

    % 3) Explicitly permitted row-order fallback.
    if strcmpi(P.anatomy.matchingFallback, 'row_order')
        n = min(nModel, nAnatomy);
        roiIdx(1:n) = (1:n)';
        matchMode = "row_order_fallback";
    end
end


function x = numericTableColumnOrNaN(T, variableName, n)
    if ismember(variableName, T.Properties.VariableNames)
        x = double(T.(variableName)(:));
        if numel(x) ~= n
            x = nan(n,1);
        end
    else
        x = nan(n,1);
    end
end

function x = numericTableColumnOrDefault(T, variableName, n, defaultValue)
    x = numericTableColumnOrNaN(T, variableName, n);
    x(~isfinite(x)) = defaultValue;
end

function [NeuronDiag, FishDiag] = buildCoordinateSystemDiagnostics(AnatomyNeurons, P)
    NeuronDiag = table();
    FishDiag = table();

    if ~P.coordinateDiagnostic.enabled || isempty(AnatomyNeurons)
        return;
    end

    required = {'xNorm','crossCandidateCenteredRawXUm', ...
        'crossRecomputedFovXUm','modelXRawUm','modelXFovUm'};
    if ~all(ismember(required, AnatomyNeurons.Properties.VariableNames))
        warning(['Coordinate diagnostic skipped because one or more required ' ...
            'coordinate columns are missing.']);
        return;
    end

    keepVars = {'morph','session','sessionKey','class3','neuronID', ...
        'anatomyRoiIdx','originalRoiIdx','roiMatchMode','fileMatchMode', ...
        'xPixels','yPixels','xNorm','xNormRawImage', ...
        'crossCandidateCenteredRawXUm','crossRecomputedFovXUm', ...
        'modelXRawUm','modelXFovUm','modelFovAngleDeg', ...
        'modelPixelLengthX','modelPixelLengthY','allCellsFile'};
    keepVars = keepVars(ismember(keepVars, AnatomyNeurons.Properties.VariableNames));
    NeuronDiag = AnatomyNeurons(:, keepVars);

    NeuronDiag.rawCoordinateErrorUm = ...
        NeuronDiag.crossCandidateCenteredRawXUm - NeuronDiag.modelXRawUm;
    NeuronDiag.fovCoordinateErrorUm = ...
        NeuronDiag.crossRecomputedFovXUm - NeuronDiag.modelXFovUm;

    morphOrder = string({P.morphs.name});
    rows = {};
    minCells = P.anatomy.minCellsPerDirectionForLROpposition;

    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        Tmorph = AnatomyNeurons(string(AnatomyNeurons.morph) == morphName, :);
        keys = unique(Tmorph.sessionKey, 'stable');

        for k = 1:numel(keys)
            T = Tmorph(Tmorph.sessionKey == keys(k), :);
            if isempty(T); continue; end

            [scoreCurrent, deltaCurrent, nCW, nCCW] = ...
                opponentAndDeltaFromX(double(T.xNorm), T.class3, minCells);
            [scoreCrossRaw, deltaCrossRaw] = ...
                opponentAndDeltaFromX(double(T.crossCandidateCenteredRawXUm), ...
                T.class3, minCells);
            [scoreCrossFov, deltaCrossFov] = ...
                opponentAndDeltaFromX(double(T.crossRecomputedFovXUm), ...
                T.class3, minCells);
            [scoreModelRaw, deltaModelRaw] = ...
                opponentAndDeltaFromX(double(T.modelXRawUm), T.class3, minCells);
            [scoreModelFov, deltaModelFov] = ...
                opponentAndDeltaFromX(double(T.modelXFovUm), T.class3, minCells);

            rSelectedVsModelRaw = safeCorrelation( ...
                T.crossCandidateCenteredRawXUm, T.modelXRawUm, ...
                P.coordinateDiagnostic.minNeuronsForCorrelation);
            rRecomputedVsModelFov = safeCorrelation( ...
                T.crossRecomputedFovXUm, T.modelXFovUm, ...
                P.coordinateDiagnostic.minNeuronsForCorrelation);
            rCurrentVsModelFov = safeCorrelation( ...
                T.xNorm, T.modelXFovUm, ...
                P.coordinateDiagnostic.minNeuronsForCorrelation);

            maeRawUm = median(abs( ...
                T.crossCandidateCenteredRawXUm - T.modelXRawUm), 'omitnan');
            maeFovUm = median(abs( ...
                T.crossRecomputedFovXUm - T.modelXFovUm), 'omitnan');

            scoreDifferenceCurrentMinusModelFov = scoreCurrent - scoreModelFov;
            scoreDifferenceRecomputedMinusModelFov = scoreCrossFov - scoreModelFov;

            roiMode = string(T.roiMatchMode(1));
            fileMode = string(T.fileMatchMode(1));
            fovAngle = median(double(T.modelFovAngleDeg), 'omitnan');

            diagnosis = diagnoseCoordinateDiscrepancy( ...
                rSelectedVsModelRaw, rRecomputedVsModelFov, ...
                scoreDifferenceCurrentMinusModelFov, ...
                scoreDifferenceRecomputedMinusModelFov, ...
                P.coordinateDiagnostic);

            rows{end+1,1} = table( ...
                morphName, string(T.session(1)), string(T.sessionKey(1)), ...
                height(T), nCW, nCCW, roiMode, fileMode, fovAngle, ...
                rSelectedVsModelRaw, rRecomputedVsModelFov, ...
                rCurrentVsModelFov, maeRawUm, maeFovUm, ...
                scoreCurrent, scoreCrossRaw, scoreCrossFov, ...
                scoreModelRaw, scoreModelFov, ...
                deltaCurrent, deltaCrossRaw, deltaCrossFov, ...
                deltaModelRaw, deltaModelFov, ...
                scoreDifferenceCurrentMinusModelFov, ...
                scoreDifferenceRecomputedMinusModelFov, diagnosis, ...
                'VariableNames', { ...
                'morph','session','sessionKey','nMatchedNeurons','nCW','nCCW', ...
                'roiMatchMode','fileMatchMode','fovAngleDeg', ...
                'rSelectedVsModelRaw','rRecomputedVsModelFov', ...
                'rCurrentNormVsModelFov','medianAbsErrorRawUm', ...
                'medianAbsErrorFovUm','opponentCurrentCrossMorph', ...
                'opponentCrossCandidateCenteredRaw', ...
                'opponentCrossRecomputedFov','opponentModelRaw', ...
                'opponentModelFov','deltaCurrentCrossMorph', ...
                'deltaCrossCandidateCenteredRawUm', ...
                'deltaCrossRecomputedFovUm','deltaModelRawUm', ...
                'deltaModelFovUm','scoreDiffCurrentMinusModelFov', ...
                'scoreDiffRecomputedMinusModelFov','diagnosis'}); %#ok<AGROW>
        end
    end

    FishDiag = vertcatTables(rows);
    if ~isempty(FishDiag)
        FishDiag.morph = categorical(FishDiag.morph, morphOrder, 'Ordinal', true);
    end
end

function [score, deltaMedian, nCW, nCCW] = opponentAndDeltaFromX(x, classes, minCells)
    x = double(x(:));
    classes = string(classes(:));

    xCW = x(classes == "CW" & isfinite(x));
    xCCW = x(classes == "CCW" & isfinite(x));
    nCW = numel(xCW);
    nCCW = numel(xCCW);

    score = NaN;
    deltaMedian = NaN;
    if nCW < minCells || nCCW < minCells
        return;
    end

    score = sideFractionWithHalfMidline(xCW, "right") + ...
        sideFractionWithHalfMidline(xCCW, "left") - 1;
    deltaMedian = median(xCW, 'omitnan') - median(xCCW, 'omitnan');
end

function r = safeCorrelation(x, y, minN)
    x = double(x(:));
    y = double(y(:));
    good = isfinite(x) & isfinite(y);
    if sum(good) < minN
        r = NaN;
        return;
    end
    if std(x(good)) == 0 || std(y(good)) == 0
        r = NaN;
        return;
    end
    C = corrcoef(x(good), y(good));
    r = C(1,2);
end

function diagnosis = diagnoseCoordinateDiscrepancy( ...
        rRaw, rFov, scoreDiffCurrent, scoreDiffRecomputed, D)

    diagnosis = "inconclusive";

    if ~isfinite(rRaw)
        diagnosis = "reference model coordinates missing or too few matched neurons";
    elseif rRaw < D.goodCorrelation
        diagnosis = "ROI identity/source mismatch likely before FOV correction";
    elseif ~isfinite(rFov) || rFov < D.goodCorrelation
        diagnosis = "FOV transform could not be reproduced from matched selected ROIs";
    elseif abs(scoreDiffRecomputed) <= D.scoreTolerance && ...
            abs(scoreDiffCurrent) > D.scoreTolerance
        diagnosis = "different centering/FOV coordinate system explains score mismatch";
    elseif abs(scoreDiffCurrent) <= D.scoreTolerance
        diagnosis = "coordinate systems agree at opponent-score level";
    elseif abs(scoreDiffRecomputed) > D.scoreTolerance
        diagnosis = "coordinate match is good but inclusion/midline rules still differ";
    else
        diagnosis = "inspect per-neuron diagnostic table";
    end
end

function plotCoordinateSystemDiagnostics(FishDiag, P)
    if ~P.coordinateDiagnostic.enabled || isempty(FishDiag)
        return;
    end

    morphOrder = string({P.morphs.name});
    morphLabels = string({P.morphs.display});

    fig = figure('Color','w','Visible',P.figureVisible, ...
        'Position',[70 80 1450 470]);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    % Panel 1: score actually used by cross-morph script versus the exact
    % FOV-corrected score already stored by the morph-specific script.
    ax1 = nexttile(tl); hold(ax1,'on');
    addIdentityAndReverseLines(ax1, [-1 1]);
    for m = 1:numel(morphOrder)
        idx = string(FishDiag.morph) == morphOrder(m);
        scatter(ax1, FishDiag.opponentModelFov(idx), ...
            FishDiag.opponentCurrentCrossMorph(idx), 55, ...
            P.morphColors(m,:), 'filled', 'MarkerEdgeColor','k', ...
            'DisplayName', morphLabels(m));
        if P.coordinateDiagnostic.addSessionLabels
            addPointLabels(ax1, FishDiag.opponentModelFov(idx), ...
                FishDiag.opponentCurrentCrossMorph(idx), FishDiag.session(idx));
        end
    end
    xlabel(ax1, 'Morph-specific score: model xCenteredUm');
    ylabel(ax1, 'Cross-morph score: current xNorm');
    title(ax1, 'Current coordinate comparison');
    axis(ax1,'square'); xlim(ax1,[-1 1]); ylim(ax1,[-1 1]);
    grid(ax1,'on'); set(ax1,'TickDir','out');
    legend(ax1,'Location','best');

    % Panel 2: same selected-file coordinates after reproducing the exact
    % candidate-median centering and FOV rotation.
    ax2 = nexttile(tl); hold(ax2,'on');
    addIdentityAndReverseLines(ax2, [-1 1]);
    for m = 1:numel(morphOrder)
        idx = string(FishDiag.morph) == morphOrder(m);
        scatter(ax2, FishDiag.opponentModelFov(idx), ...
            FishDiag.opponentCrossRecomputedFov(idx), 55, ...
            P.morphColors(m,:), 'filled', 'MarkerEdgeColor','k', ...
            'HandleVisibility','off');
        if P.coordinateDiagnostic.addSessionLabels
            addPointLabels(ax2, FishDiag.opponentModelFov(idx), ...
                FishDiag.opponentCrossRecomputedFov(idx), FishDiag.session(idx));
        end
    end
    xlabel(ax2, 'Morph-specific score: model xCenteredUm');
    ylabel(ax2, 'Recomputed score from selected ROIs + FOV');
    title(ax2, 'After reproducing FOV correction');
    axis(ax2,'square'); xlim(ax2,[-1 1]); ylim(ax2,[-1 1]);
    grid(ax2,'on'); set(ax2,'TickDir','out');

    % Panel 3: neuron-level correspondence within each fish.
    ax3 = nexttile(tl); hold(ax3,'on');
    offsets = linspace(-0.12,0.12,numel(morphOrder));
    for m = 1:numel(morphOrder)
        idx = string(FishDiag.morph) == morphOrder(m);
        n = sum(idx);
        scatter(ax3, 1 + offsets(m) + zeros(n,1), ...
            FishDiag.rSelectedVsModelRaw(idx), 48, P.morphColors(m,:), ...
            'filled', 'MarkerEdgeColor','k', 'HandleVisibility','off');
        scatter(ax3, 2 + offsets(m) + zeros(n,1), ...
            FishDiag.rRecomputedVsModelFov(idx), 48, P.morphColors(m,:), ...
            'filled', 'MarkerEdgeColor','k', 'HandleVisibility','off');
    end
    yline(ax3, P.coordinateDiagnostic.goodCorrelation, 'k:', ...
        'LineWidth',1.2,'HandleVisibility','off');
    xlim(ax3,[0.5 2.5]); ylim(ax3,[-1.05 1.05]);
    xticks(ax3,[1 2]);
    xticklabels(ax3,{'Selected vs model raw','Recomputed vs model FOV'});
    ylabel(ax3,'Neuron-level Pearson r within fish');
    title(ax3,'ROI matching and transform checks');
    grid(ax3,'on'); set(ax3,'TickDir','out');

    title(tl,'Diagnostic: are the two analyses using the same x coordinate?');
    subtitle(tl,['Panel 1 should lie on y=x only if current xNorm matches ' ...
        'the morph-specific FOV-corrected coordinate. Panel 2 tests the exact reconstruction.']);

    saveFigure(fig, P.figureDir, ...
        '13_coordinate_system_diagnostic_cross_vs_morph_specific', P);
end

function addIdentityAndReverseLines(ax, lim)
    plot(ax, lim, lim, 'k--', 'LineWidth',1.2, ...
        'DisplayName','y = x');
    plot(ax, lim, -lim, ':', 'Color',[0.55 0.55 0.55], ...
        'LineWidth',1.0, 'DisplayName','y = -x');
    xline(ax,0,'k:','HandleVisibility','off');
    yline(ax,0,'k:','HandleVisibility','off');
end

function addPointLabels(ax, x, y, labels)
    for i = 1:numel(x)
        if isfinite(x(i)) && isfinite(y(i))
            text(ax, x(i)+0.025, y(i)+0.025, string(labels(i)), ...
                'FontSize',7,'Interpreter','none');
        end
    end
end

function FishImagingTotals = buildFishImagingTotals( ...
        AnatomyNeurons, ModelSourceIndex, P)
    % Locate one ORIGINAL *_RASTER.mat per fish and read only the variable
    % metadata. The raster itself is not loaded into memory.

    rows = {};
    morphOrder = string({P.morphs.name});

    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        rootDir = string(P.morphs(m).rootDir);
        Tmorph = AnatomyNeurons(string(AnatomyNeurons.morph) == morphName,:);
        keys = unique(string(Tmorph.sessionKey),'stable');

        rasterFiles = discoverOriginalRasterFiles(rootDir,P);

        for k = 1:numel(keys)
            key = keys(k);
            Tfish = Tmorph(string(Tmorph.sessionKey) == key,:);
            sessionName = string(Tfish.session(1));
            nClassifiedNeuronsInAnalysis = height(Tfish);

            sourceFolder = "";
            if ~isempty(ModelSourceIndex)
                src = string(ModelSourceIndex.morph) == morphName & ...
                    string(ModelSourceIndex.sessionKey) == key;
                srcIdx = find(src,1,'first');
                if ~isempty(srcIdx)
                    originalAllCellsPath = string( ...
                        ModelSourceIndex.originalAllCellsPath(srcIdx));
                    if strlength(originalAllCellsPath) > 0
                        sourceFolder = string(fileparts(char(originalAllCellsPath)));
                    end
                end
            end

            [selectedFile,nCandidates,matchMode] = selectRasterForFish( ...
                rasterFiles,key,sourceFolder);
            [nTotal,status] = readTotalNeuronCount(selectedFile,P);

            if isfinite(nTotal) && nTotal < nClassifiedNeuronsInAnalysis
                status = "INVALID_TOTAL_SMALLER_THAN_CLASSIFIED_COUNT";
                nTotal = NaN;
            end

            if ~isfinite(nTotal)
                warning(['No valid original RASTER neuron count for %s | %s ' ...
                    '(status: %s). Count/rate tests will omit this fish.'], ...
                    morphName,sessionName,status);
            end

            rows{end+1,1} = table(morphName,sessionName,key, ...
                nClassifiedNeuronsInAnalysis,nTotal, ...
                selectedFile,nCandidates,matchMode,status, ...
                'VariableNames',{'morph','session','sessionKey', ...
                'nClassifiedNeuronsInAnalysis','nTotalImagedNeurons', ...
                'originalRasterFile', ...
                'nOriginalRasterCandidates','rasterMatchMode', ...
                'rasterCountStatus'}); %#ok<AGROW>
        end
    end

    FishImagingTotals = vertcatTables(rows);
    if P.totalRaster.failOnMissing && ...
            any(~isfinite(FishImagingTotals.nTotalImagedNeurons))
        error(['At least one fish has no unique valid original *_RASTER.mat. ' ...
            'Inspect results.fishImagingTotals.']);
    end
end

function files = discoverOriginalRasterFiles(rootDir,P)
    files = strings(0,1);
    if ~isfolder(rootDir)
        warning('Morph root folder does not exist: %s',rootDir);
        return;
    end

    d = dir(fullfile(char(rootDir),'**',P.totalRaster.pattern));
    if isempty(d); return; end

    fullPaths = string(fullfile({d.folder},{d.name}))';
    upperPaths = upper(fullPaths);
    keep = true(size(fullPaths));
    for q = 1:numel(P.totalRaster.excludeTokens)
        keep = keep & ~contains(upperPaths, ...
            upper(string(P.totalRaster.excludeTokens{q})));
    end
    files = unique(fullPaths(keep),'stable');
end

function [selectedFile,nCandidates,matchMode] = selectRasterForFish( ...
        rasterFiles,sessionKey,sourceFolder)
    selectedFile = "";
    nCandidates = 0;
    matchMode = "none";
    if isempty(rasterFiles); return; end

    candidateMask = false(size(rasterFiles));

    % Strongest match: the RASTER is in the same recording folder as the
    % original ALL_CELLS file recorded in Results.paths.allCellsPath.
    if strlength(sourceFolder) > 0
        rasterFolders = strings(size(rasterFiles));
        for i = 1:numel(rasterFiles)
            rasterFolders(i) = string(fileparts(char(rasterFiles(i))));
        end
        candidateMask = strcmpi(rasterFolders,sourceFolder);
        if any(candidateMask); matchMode = "same_source_folder"; end
    end

    % Fallback: match the recXX key to the immediate parent folder.
    if ~any(candidateMask)
        for i = 1:numel(rasterFiles)
            folder = string(fileparts(char(rasterFiles(i))));
            [~,folderBase] = fileparts(char(folder));
            candidateMask(i) = makeSessionKey(string(folderBase)) == sessionKey;
        end
        if any(candidateMask); matchMode = "session_folder_key"; end
    end

    % Last conservative fallback: the session key occurs in the full path.
    if ~any(candidateMask) && sessionKey ~= "unknown"
        candidateMask = contains(lower(rasterFiles),lower(sessionKey));
        if any(candidateMask); matchMode = "session_key_in_path"; end
    end

    candidates = rasterFiles(candidateMask);
    nCandidates = numel(candidates);
    if nCandidates == 1
        selectedFile = candidates(1);
    elseif nCandidates > 1
        matchMode = matchMode + "_AMBIGUOUS";
    end
end

function [nTotal,status] = readTotalNeuronCount(rasterFile,P)
    nTotal = NaN;
    status = "NO_UNIQUE_FILE";
    if strlength(rasterFile) == 0; return; end
    if ~isfile(rasterFile)
        status = "FILE_NOT_FOUND";
        return;
    end

    try
        info = whos('-file',char(rasterFile),P.totalRaster.neuronVariable);
    catch ME
        status = "WHOS_FAILED: " + string(ME.message);
        return;
    end
    if isempty(info)
        status = "MISSING_VARIABLE_" + string(P.totalRaster.neuronVariable);
        return;
    end

    dims = double(info(1).size);
    neuronDim = P.totalRaster.neuronDimension;
    if numel(dims) < neuronDim || dims(neuronDim) < 1
        status = "INVALID_RASTER_DIMENSIONS";
        return;
    end

    nTotal = dims(neuronDim);
    status = "OK";
end

function FishClassMetrics = attachImagingTotals( ...
        FishClassMetrics,FishImagingTotals)
    nRows = height(FishClassMetrics);
    FishClassMetrics.nTotalImagedNeurons = nan(nRows,1);
    FishClassMetrics.pctOfAllImagedNeurons = nan(nRows,1);
    FishClassMetrics.originalRasterFile = strings(nRows,1);
    FishClassMetrics.rasterCountStatus = strings(nRows,1);

    for i = 1:nRows
        idx = string(FishImagingTotals.morph) == ...
            string(FishClassMetrics.morph(i)) & ...
            string(FishImagingTotals.sessionKey) == ...
            string(FishClassMetrics.sessionKey(i));
        j = find(idx,1,'first');
        if isempty(j); continue; end

        total = double(FishImagingTotals.nTotalImagedNeurons(j));
        FishClassMetrics.nTotalImagedNeurons(i) = total;
        FishClassMetrics.originalRasterFile(i) = ...
            string(FishImagingTotals.originalRasterFile(j));
        FishClassMetrics.rasterCountStatus(i) = ...
            string(FishImagingTotals.rasterCountStatus(j));
        if isfinite(total) && total > 0
            FishClassMetrics.pctOfAllImagedNeurons(i) = ...
                100 * FishClassMetrics.nCells(i) / total;
        end
    end
end

function FishOverallMetrics = buildFishOverallMetrics( ...
        AnatomyNeurons,FishImagingTotals,ModelSourceIndex,P)
    rows = {};
    morphOrder = string({P.morphs.name});
    populationOrder = ["All cells" "Highly tuned cells"];

    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        Tmorph = AnatomyNeurons(string(AnatomyNeurons.morph) == morphName,:);
        keys = unique(string(Tmorph.sessionKey),'stable');

        for k = 1:numel(keys)
            key = keys(k);
            Tfish = Tmorph(string(Tmorph.sessionKey) == key,:);
            if isempty(Tfish); continue; end
            sessionName = string(Tfish.session(1));

            totalIdx = string(FishImagingTotals.morph) == morphName & ...
                string(FishImagingTotals.sessionKey) == key;
            totalRow = find(totalIdx,1,'first');
            nTotal = NaN;
            if ~isempty(totalRow)
                nTotal = double(FishImagingTotals.nTotalImagedNeurons(totalRow));
            end

            sourceRow = [];
            requiredSourceVars = {'morph','sessionKey','originalAllCellsPath'};
            if ~isempty(ModelSourceIndex) && all(ismember(requiredSourceVars, ...
                    ModelSourceIndex.Properties.VariableNames))
                sourceIdx = string(ModelSourceIndex.morph) == morphName & ...
                    string(ModelSourceIndex.sessionKey) == key;
                sourceRow = find(sourceIdx,1,'first');
            end
            allCellsPath = "";
            if ~isempty(sourceRow)
                allCellsPath = string( ...
                    ModelSourceIndex.originalAllCellsPath(sourceRow));
            end

            [allXUm,allYUm,allCellScaleUm,compactnessStatus] = ...
                loadOverallAllCellCoordinates(allCellsPath);
            nAllSpatial = sum(isfinite(allXUm) & isfinite(allYUm));
            allRadius = normalizedRadiusOfGyration( ...
                allXUm,allYUm,allCellScaleUm);

            tunedXUm = double(Tfish.xUm);
            tunedYUm = double(Tfish.yUm);
            tunedGood = isfinite(tunedXUm) & isfinite(tunedYUm);
            nTuned = sum(tunedGood);
            tunedRadius = normalizedRadiusOfGyration( ...
                tunedXUm,tunedYUm,allCellScaleUm);

            rows{end+1,1} = table(morphName,sessionName,key, ...
                populationOrder(1),nTotal,nTotal,nAllSpatial,allRadius, ...
                allCellsPath,compactnessStatus, ...
                'VariableNames',{'morph','session','sessionKey', ...
                'populationName','nCells','nTotalImagedNeurons', ...
                'nSpatialCells','radiusOfGyrationNorm', ...
                'originalAllCellsFile','compactnessStatus'}); %#ok<AGROW>
            rows{end+1,1} = table(morphName,sessionName,key, ...
                populationOrder(2),nTuned,nTotal,nTuned,tunedRadius, ...
                allCellsPath,compactnessStatus, ...
                'VariableNames',{'morph','session','sessionKey', ...
                'populationName','nCells','nTotalImagedNeurons', ...
                'nSpatialCells','radiusOfGyrationNorm', ...
                'originalAllCellsFile','compactnessStatus'}); %#ok<AGROW>
        end
    end

    FishOverallMetrics = vertcatTables(rows);
    FishOverallMetrics.morph = categorical( ...
        FishOverallMetrics.morph,morphOrder,'Ordinal',true);
    FishOverallMetrics.populationName = categorical( ...
        FishOverallMetrics.populationName,populationOrder,'Ordinal',true);
end

function [xUm,yUm,scaleUm,status] = ...
        loadOverallAllCellCoordinates(allCellsPath)
    xUm = [];
    yUm = [];
    scaleUm = NaN;
    status = "missing original ALL_CELLS file";
    if strlength(allCellsPath) == 0 || ~isfile(allCellsPath)
        return;
    end

    try
        info = whos('-file',char(allCellsPath),'cell_per','cells');
        names = string({info.name});
        idx = find(names == "cell_per",1,'first');
        if isempty(idx)
            idx = find(names == "cells",1,'first');
        end
        if isempty(idx)
            status = "original ALL_CELLS lacks cell_per/cells";
            return;
        end
        nAll = prod(double(info(idx).size));
        [xPx,yPx,validID,pixelLengthX,pixelLengthY] = ...
            extractOriginalRoiCentroidsForAudit( ...
            char(allCellsPath),(1:nAll)');
        good = validID & isfinite(xPx) & isfinite(yPx);
        if ~any(good)
            status = "no valid all-cell ROI centroids";
            return;
        end
        if ~isfinite(pixelLengthX) || pixelLengthX <= 0 || ...
                ~isfinite(pixelLengthY) || pixelLengthY <= 0
            status = "missing pixel dimensions in original ALL_CELLS";
            return;
        end

        xUm = (xPx(good) - median(xPx(good),'omitnan')) .* pixelLengthX;
        yUm = (yPx(good) - median(yPx(good),'omitnan')) .* pixelLengthY;
        scaleUm = max([abs(xUm); abs(yUm)],[],'omitnan');
        if ~isfinite(scaleUm) || scaleUm <= 0
            scaleUm = NaN;
            status = "invalid all-cell anatomical scale";
        else
            status = "OK";
        end
    catch ME
        xUm = [];
        yUm = [];
        scaleUm = NaN;
        status = "ALL_CELLS read failed: " + string(ME.message);
        warning('Overall compactness unavailable for %s: %s', ...
            allCellsPath,ME.message);
    end
end

function radius = normalizedRadiusOfGyration(x,y,scaleUm)
    x = double(x(:));
    y = double(y(:));
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);
    radius = NaN;
    if isempty(x) || ~isfinite(scaleUm) || scaleUm <= 0
        return;
    end
    centroidX = mean(x,'omitnan');
    centroidY = mean(y,'omitnan');
    distanceSquared = (x-centroidX).^2 + (y-centroidY).^2;
    radius = sqrt(mean(distanceSquared,'omitnan')) ./ scaleUm;
end

function StatsTable = runOverallPopulationStatistics(FishOverallMetrics,P)
    rows = {};
    morphOrder = string({P.morphs.name});
    populationOrder = ["All cells" "Highly tuned cells"];
    metricOrder = ["nCells" "radiusOfGyrationNorm"];
    pairs = nchoosek(1:numel(morphOrder),2);

    for mi = 1:numel(metricOrder)
        metricName = metricOrder(mi);
        for pi = 1:numel(populationOrder)
            populationName = populationOrder(pi);
            T = FishOverallMetrics(string(FishOverallMetrics.populationName) == ...
                populationName,:);
            y = double(T.(metricName));
            valid = isfinite(y);
            if metricName == "radiusOfGyrationNorm"
                valid = valid & double(T.nSpatialCells) > 0;
            end

            pRaw = nan(size(pairs,1),1);
            tempRows = cell(size(pairs,1),1);
            for p = 1:size(pairs,1)
                g1 = morphOrder(pairs(p,1));
                g2 = morphOrder(pairs(p,2));
                pairMask = valid & ismember(string(T.morph),[g1 g2]);
                Tp = T(pairMask,:);
                yp = double(Tp.(metricName));
                group = string(Tp.morph) == g2;

                if metricName == "nCells"
                    exposure = ones(size(yp));
                    [pValue,~,~,~,~,~] = quasiPoissonRateTest( ...
                        yp,exposure,group,P.stats.minFishPerMorph);
                    testMethod = "quasi-Poisson fish-level count GLM";
                else
                    weights = normalizePrecisionWeights( ...
                        double(Tp.nSpatialCells),P);
                    [pValue,~,~,~,~] = weightedPermutationTest( ...
                        yp,weights,group,P, ...
                        P.stats.permutationRandomSeed + 800000 + ...
                        10000*mi + 100*pi + p);
                    testMethod = "weighted fish-level permutation";
                end
                pRaw(p) = pValue;
                tempRows{p} = table(metricName,populationName,g1,g2, ...
                    pairs(p,1),pairs(p,2),sum(~group),sum(group), ...
                    pValue,NaN,testMethod, ...
                    'VariableNames',{'metric','populationName', ...
                    'group1','group2','group1Index', ...
                    'group2Index','n1','n2','pRaw', ...
                    'pAdjusted','testMethod'});
            end

            for p = 1:numel(tempRows)
                % Retain the legacy column name for output compatibility.
                tempRows{p}.pAdjusted = pRaw(p);
                rows{end+1,1} = tempRows{p}; %#ok<AGROW>
            end
        end
    end
    StatsTable = vertcatTables(rows);
end

function plotOverallNumberOfCells(FishOverallMetrics,StatsTable,P)
    plotOverallPopulationMetric(FishOverallMetrics,StatsTable, ...
        "nCells","Number of cells", ...
        "17_number_of_cells_overall_across_morphs",P);
end

function plotOverallPopulationCompactness(FishOverallMetrics,StatsTable,P)
    plotOverallPopulationMetric(FishOverallMetrics,StatsTable, ...
        "radiusOfGyrationNorm", ...
        "Radius of gyration (all-cell-normalized units; lower = more compact)", ...
        "18_overall_population_compactness_across_morphs",P);
end

function plotOverallPopulationMetric( ...
        FishOverallMetrics,StatsTable,metricName,yLabel,fileBase,P)
    populationOrder = ["All cells" "Highly tuned cells"];
    morphOrder = string({P.morphs.name});
    morphLabels = string({P.morphs.display});

    fig = figure('Color','w','Visible',P.figureVisible, ...
        'Position',[80 100 980 470]);
    tl = tiledlayout(1,numel(populationOrder), ...
        'TileSpacing','compact','Padding','compact');

    for pi = 1:numel(populationOrder)
        populationName = populationOrder(pi);
        ax = nexttile(tl,pi);
        hold(ax,'on');
        T = FishOverallMetrics(string(FishOverallMetrics.populationName) == ...
            populationName,:);
        valuesAll = double(T.(metricName));
        validAll = isfinite(valuesAll);

        for m = 1:numel(morphOrder)
            idx = validAll & string(T.morph) == morphOrder(m);
            drawBoxAndFishPoints(ax,m,valuesAll(idx), ...
                P.morphColors(m,:),T.session(idx),P);
        end

        xlim(ax,[0.45 numel(morphOrder)+0.55]);
        xticks(ax,1:numel(morphOrder));
        xticklabels(ax,morphLabels);
        ylabel(ax,yLabel,'Interpreter','none');
        title(ax,populationName,'Interpreter','none');
        set(ax,'TickDir','out','LineWidth',1,'Box','off', ...
            'FontSize',10);
        grid(ax,'on');
        ax.XGrid = 'off';

        S = StatsTable(StatsTable.metric == metricName & ...
            string(StatsTable.populationName) == populationName,:);
        addPairwiseStatBars(ax,S,valuesAll(validAll),P);
    end

    title(tl,strrep(fileBase,'_',' '), ...
        'Interpreter','none','FontWeight','bold');
    subtitle(tl,['Fish-level comparisons across morphs; ' ...
        'unadjusted pairwise p-values']);
    saveFigure(fig,P.figureDir,char(fileBase),P);
end

function FishClassMetrics = buildFishClassMetrics(AnatomyNeurons, P)
    rows = {};
    morphOrder = string({P.morphs.name});
    classOrder = string(P.classes);

    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        Tmorph = AnatomyNeurons(string(AnatomyNeurons.morph) == morphName, :);
        keys = unique(Tmorph.sessionKey, 'stable');

        for k = 1:numel(keys)
            key = keys(k);
            Tfish = Tmorph(Tmorph.sessionKey == key, :);
            if isempty(Tfish); continue; end
            sessionName = Tfish.session(1);

            allX = double(Tfish.xNorm);
            allY = double(Tfish.yNorm);
            goodAll = isfinite(allX) & isfinite(allY);
            nAll = sum(goodAll);
            referenceHullArea = convexHullAreaSafe(allX(goodAll), allY(goodAll), ...
                P.anatomy.minCellsForHull);

            for c = 1:numel(classOrder)
                className = classOrder(c);
                T = Tfish(string(Tfish.class3) == className, :);
                x = double(T.xNorm);
                y = double(T.yNorm);
                good = isfinite(x) & isfinite(y);
                x = x(good);
                y = y(good);
                n = numel(x);

                pctClassified = safePercent(n, nAll);
                densityPerNormalizedFovArea = n / 4; % normalized FOV is 2 x 2
                densityPerReferenceHullArea = safeDivide(n, referenceHullArea);

                centroidX = mean(x, 'omitnan');
                centroidY = mean(y, 'omitnan');
                medianX = median(x, 'omitnan');
                medianY = median(y, 'omitnan');
                xIQR = iqrFinite(x);
                yIQR = iqrFinite(y);
                xSpan = spanFinite(x, P.anatomy.minCellsForSpread);
                ySpan = spanFinite(y, P.anatomy.minCellsForSpread);
                xOccupancy = axisOccupancyFraction(x, P.anatomy.axisEdges);
                yOccupancy = axisOccupancyFraction(y, P.anatomy.axisEdges);

                radiusOfGyration = NaN;
                medianDistanceToCentroid = NaN;
                if n >= 1 && isfinite(centroidX) && isfinite(centroidY)
                    dCentroid = sqrt((x - centroidX).^2 + (y - centroidY).^2);
                    radiusOfGyration = sqrt(mean(dCentroid.^2, 'omitnan'));
                    medianDistanceToCentroid = median(dCentroid, 'omitnan');
                end

                hullArea = convexHullAreaSafe(x, y, P.anatomy.minCellsForHull);
                hullAreaPerCell = safeDivide(hullArea, n);

                medianNN = medianNearestNeighbourDistance(x, y, ...
                    P.anatomy.minCellsForNearestNeighbour);
                localDensityIndex = safeDivide(1, pi * medianNN.^2);

                fractionXNegative = safeDivide(sum(x < 0), n);
                fractionXPositive = safeDivide(sum(x > 0), n);
                fractionYNegative = safeDivide(sum(y < 0), n);
                fractionYPositive = safeDivide(sum(y > 0), n);

                rows{end+1,1} = table(morphName, sessionName, key, className, ...
                    n, nAll, pctClassified, referenceHullArea, ...
                    densityPerNormalizedFovArea, densityPerReferenceHullArea, ...
                    centroidX, centroidY, medianX, medianY, xIQR, yIQR, ...
                    xSpan, ySpan, xOccupancy, yOccupancy, ...
                    radiusOfGyration, medianDistanceToCentroid, ...
                    hullArea, hullAreaPerCell, medianNN, localDensityIndex, ...
                    fractionXNegative, fractionXPositive, ...
                    fractionYNegative, fractionYPositive, ...
                    'VariableNames', {'morph','session','sessionKey','className', ...
                    'nCells','nAllClassifiedCells','pctOfClassifiedCells', ...
                    'referenceHullAreaNorm','densityPerNormalizedFovArea', ...
                    'densityPerReferenceHullArea','centroidXNorm','centroidYNorm', ...
                    'medianXNorm','medianYNorm','iqrXNorm','iqrYNorm', ...
                    'xSpanNorm','ySpanNorm','xOccupancyFraction','yOccupancyFraction', ...
                    'radiusOfGyrationNorm','medianDistanceToCentroidNorm', ...
                    'convexHullAreaNorm','hullAreaPerCellNorm', ...
                    'medianNearestNeighbourDistanceNorm','localDensityIndex', ...
                    'fractionXNegative','fractionXPositive', ...
                    'fractionYNegative','fractionYPositive'}); %#ok<AGROW>
            end
        end
    end

    FishClassMetrics = vertcatTables(rows);
    FishClassMetrics.morph = categorical(FishClassMetrics.morph, morphOrder, 'Ordinal', true);
    FishClassMetrics.className = categorical(FishClassMetrics.className, classOrder, 'Ordinal', true);
end


function FishOpponentMetrics = buildFishOpponentMetrics(AnatomyNeurons, P)
    rows = {};
    morphOrder = string({P.morphs.name});
    edges = double(P.anatomy.axisEdges(:)');
    minCells = P.anatomy.minCellsPerDirectionForLROpposition;

    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        Tmorph = AnatomyNeurons(string(AnatomyNeurons.morph) == morphName, :);
        keys = unique(Tmorph.sessionKey, 'stable');

        for k = 1:numel(keys)
            key = keys(k);
            Tfish = Tmorph(Tmorph.sessionKey == key, :);
            if isempty(Tfish); continue; end
            sessionName = Tfish.session(1);

            xCW = double(Tfish.xNorm(string(Tfish.class3) == "CW"));
            xCCW = double(Tfish.xNorm(string(Tfish.class3) == "CCW"));
            xCW = xCW(isfinite(xCW));
            xCCW = xCCW(isfinite(xCCW));

            nCW = numel(xCW);
            nCCW = numel(xCCW);
            enoughCells = nCW >= minCells && nCCW >= minCells;

            pCWright = NaN;
            pCWleft = NaN;
            pCCWright = NaN;
            pCCWleft = NaN;
            opponentLR = NaN;
            medianXCW = NaN;
            medianXCCW = NaN;
            deltaMedianX_CWminusCCW = NaN;
            densityOverlap = NaN;
            densitySegregation = NaN;
            signedDensitySegregation = NaN;

            if enoughCells
                % Cells exactly on x=0 contribute half to each side rather
                % than being discarded.
                pCWright = sideFractionWithHalfMidline(xCW, "right");
                pCWleft = sideFractionWithHalfMidline(xCW, "left");
                pCCWright = sideFractionWithHalfMidline(xCCW, "right");
                pCCWleft = sideFractionWithHalfMidline(xCCW, "left");

                % +1 = all CW right and all CCW left.
                %  0 = no net opponent left-right organization.
                % -1 = completely reversed organization.
                opponentLR = pCWright + pCCWleft - 1;

                medianXCW = median(xCW, 'omitnan');
                medianXCCW = median(xCCW, 'omitnan');

                % Positive means that the CW population lies to the right of
                % the CCW population.
                deltaMedianX_CWminusCCW = medianXCW - medianXCCW;

                % Compare the full fish-level CW and CCW distributions.
                % Histograms are normalized to probability mass, smoothed,
                % and renormalized before calculating their overlap.
                dCW = normalizedSmoothedHistogram(xCW, edges, ...
                    P.anatomy.densitySmoothingSigmaBins);
                dCCW = normalizedSmoothedHistogram(xCCW, edges, ...
                    P.anatomy.densitySmoothingSigmaBins);

                if all(isfinite(dCW)) && all(isfinite(dCCW))
                    densityOverlap = sum(min(dCW, dCCW));
                    densityOverlap = min(max(densityOverlap, 0), 1);
                    densitySegregation = 1 - densityOverlap;

                    % Direction is supplied by the continuous median shift.
                    % Positive values correspond to CW-right / CCW-left.
                    signedDensitySegregation = ...
                        sign(deltaMedianX_CWminusCCW) * densitySegregation;
                end
            end

            rows{end+1,1} = table(morphName, sessionName, key, ...
                nCW, nCCW, enoughCells, ...
                pCWleft, pCWright, pCCWleft, pCCWright, ...
                opponentLR, medianXCW, medianXCCW, ...
                deltaMedianX_CWminusCCW, densityOverlap, ...
                densitySegregation, signedDensitySegregation, ...
                'VariableNames', {'morph','session','sessionKey', ...
                'nCW','nCCW','enoughCells', ...
                'fractionCWLeft','fractionCWRight', ...
                'fractionCCWLeft','fractionCCWRight', ...
                'opponentLR','medianXCWNorm','medianXCCWNorm', ...
                'deltaMedianX_CWminusCCW', ...
                'densityOverlap','densitySegregation', ...
                'signedDensitySegregation'}); %#ok<AGROW>
        end
    end

    FishOpponentMetrics = vertcatTables(rows);
    if ~isempty(FishOpponentMetrics)
        FishOpponentMetrics.morph = categorical(FishOpponentMetrics.morph, ...
            morphOrder, 'Ordinal', true);
    end
end

function p = sideFractionWithHalfMidline(x, requestedSide)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        p = NaN;
        return;
    end

    nMidline = sum(x == 0);
    switch lower(string(requestedSide))
        case "right"
            p = (sum(x > 0) + 0.5*nMidline) / numel(x);
        case "left"
            p = (sum(x < 0) + 0.5*nMidline) / numel(x);
        otherwise
            error('requestedSide must be "left" or "right".');
    end
end

function h = normalizedSmoothedHistogram(values, edges, sigmaBins)
    values = double(values(:));
    values = values(isfinite(values));
    counts = histcounts(values, edges);

    if isempty(values) || sum(counts) <= 0
        h = nan(1, numel(edges)-1);
        return;
    end

    h = counts ./ sum(counts);
    h = smoothVectorGaussian(h, sigmaBins);
    total = sum(h, 'omitnan');

    if ~isfinite(total) || total <= 0
        h(:) = NaN;
    else
        h = h ./ total;
    end
end

function FishCentralXMetrics = buildFishCentralXMetrics(AnatomyNeurons, P)
    % Compute central-region probabilities and fixed-bandwidth KDE peaks.
    % Every output row is one fish/recording.

    rows = {};
    morphOrder = string({P.morphs.name});
    C = P.anatomy.centralX;

    validateCentralXParameters(C);

    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        Tmorph = AnatomyNeurons(string(AnatomyNeurons.morph) == morphName, :);
        keys = unique(string(Tmorph.sessionKey), 'stable');

        for k = 1:numel(keys)
            key = keys(k);
            Tfish = Tmorph(string(Tmorph.sessionKey) == key, :);
            if isempty(Tfish); continue; end

            sessionName = string(Tfish.session(1));
            xCCW = finiteClassXUm(Tfish, "CCW");
            xCW = finiteClassXUm(Tfish, "CW");
            xSym = finiteClassXUm(Tfish, "Symmetric");

            nCCW = numel(xCCW);
            nCW = numel(xCW);
            nSymmetric = numel(xSym);

            enoughCCWForFraction = nCCW >= C.minCellsForFraction;
            enoughCWForFraction = nCW >= C.minCellsForFraction;
            enoughSymmetricForFraction = nSymmetric >= C.minCellsForFraction;
            enoughDirectionalForPeaks = ...
                nCCW >= C.minCellsForPeak && nCW >= C.minCellsForPeak;

            ccwCentralLeftFraction = NaN;
            cwCentralRightFraction = NaN;
            opponentCentralEnrichment = NaN;
            symmetricCentralFraction = NaN;
            nCCWCentralLeft = 0;
            nCWCentralRight = 0;
            nSymmetricCentral = 0;

            % Calculate the fractions for every fish that contains at least
            % one neuron in the relevant class. The enough* flags below are
            % used only to decide which fish enter morph-level statistics.
            if nCCW > 0
                nCCWCentralLeft = sum( ...
                    xCCW >= C.ccwRangeUm(1) & xCCW < C.ccwRangeUm(2));
                ccwCentralLeftFraction = nCCWCentralLeft / nCCW;
            end
            if nCW > 0
                nCWCentralRight = sum( ...
                    xCW > C.cwRangeUm(1) & xCW <= C.cwRangeUm(2));
                cwCentralRightFraction = nCWCentralRight / nCW;
            end
            if isfinite(ccwCentralLeftFraction) && ...
                    isfinite(cwCentralRightFraction)
                opponentCentralEnrichment = ...
                    ccwCentralLeftFraction + cwCentralRightFraction;
            end
            if nSymmetric > 0
                nSymmetricCentral = sum( ...
                    xSym >= C.symmetricRangeUm(1) & ...
                    xSym <= C.symmetricRangeUm(2));
                symmetricCentralFraction = nSymmetricCentral / nSymmetric;
            end

            xPeakCCWUm = NaN;
            xPeakCCW_CILowUm = NaN;
            xPeakCCW_CIHighUm = NaN;
            xPeakCWUm = NaN;
            xPeakCW_CILowUm = NaN;
            xPeakCW_CIHighUm = NaN;
            deltaXPeakUm = NaN;
            deltaXPeak_CILowUm = NaN;
            deltaXPeak_CIHighUm = NaN;

            if enoughDirectionalForPeaks
                seed = C.bootstrapRandomSeed + 10000*m + k;
                [xPeakCCWUm, xPeakCCW_CILowUm, xPeakCCW_CIHighUm, ...
                    bootCCW] = estimateKdePeakWithBootstrap(xCCW, C, seed);
                [xPeakCWUm, xPeakCW_CILowUm, xPeakCW_CIHighUm, ...
                    bootCW] = estimateKdePeakWithBootstrap(xCW, C, seed + 1000);

                if isfinite(xPeakCCWUm) && isfinite(xPeakCWUm)
                    deltaXPeakUm = xPeakCWUm - xPeakCCWUm;
                    deltaBoot = bootCW - bootCCW;
                    deltaBoot = deltaBoot(isfinite(deltaBoot));
                    if ~isempty(deltaBoot)
                        ciTail = (100 - C.bootstrapCI) / 2;
                        deltaCI = prctile(deltaBoot, [ciTail 100-ciTail]);
                        deltaXPeak_CILowUm = deltaCI(1);
                        deltaXPeak_CIHighUm = deltaCI(2);
                    end
                end
            end

            rows{end+1,1} = table( ...
                morphName, sessionName, key, nCCW, nCW, nSymmetric, ...
                nCCWCentralLeft, nCWCentralRight, nSymmetricCentral, ...
                enoughCCWForFraction, enoughCWForFraction, ...
                enoughSymmetricForFraction, enoughDirectionalForPeaks, ...
                ccwCentralLeftFraction, cwCentralRightFraction, ...
                opponentCentralEnrichment, symmetricCentralFraction, ...
                xPeakCCWUm, xPeakCCW_CILowUm, xPeakCCW_CIHighUm, ...
                xPeakCWUm, xPeakCW_CILowUm, xPeakCW_CIHighUm, ...
                deltaXPeakUm, deltaXPeak_CILowUm, deltaXPeak_CIHighUm, ...
                C.kdeBandwidthUm, C.kdeSearchRangeUm(1), ...
                C.kdeSearchRangeUm(2), C.bootstrapIterations, ...
                C.bootstrapCI, ...
                'VariableNames', {'morph','session','sessionKey', ...
                'nCCW','nCW','nSymmetric', ...
                'nCCWCentralLeft','nCWCentralRight','nSymmetricCentral', ...
                'enoughCCWForFraction','enoughCWForFraction', ...
                'enoughSymmetricForFraction','enoughDirectionalForPeaks', ...
                'ccwCentralLeftFraction','cwCentralRightFraction', ...
                'opponentCentralEnrichment','symmetricCentralFraction', ...
                'xPeakCCWUm','xPeakCCW_CILowUm','xPeakCCW_CIHighUm', ...
                'xPeakCWUm','xPeakCW_CILowUm','xPeakCW_CIHighUm', ...
                'deltaXPeakUm','deltaXPeak_CILowUm','deltaXPeak_CIHighUm', ...
                'kdeBandwidthUm','kdeSearchLowUm','kdeSearchHighUm', ...
                'bootstrapIterations','bootstrapCI'}); %#ok<AGROW>
        end
    end

    FishCentralXMetrics = vertcatTables(rows);
    if ~isempty(FishCentralXMetrics)
        FishCentralXMetrics.morph = categorical( ...
            FishCentralXMetrics.morph, morphOrder, 'Ordinal', true);
    end
end

function validateCentralXParameters(C)
    assert(numel(C.ccwRangeUm) == 2 && ...
        C.ccwRangeUm(1) < C.ccwRangeUm(2), ...
        'P.anatomy.centralX.ccwRangeUm must be [low high].');
    assert(numel(C.cwRangeUm) == 2 && ...
        C.cwRangeUm(1) < C.cwRangeUm(2), ...
        'P.anatomy.centralX.cwRangeUm must be [low high].');
    assert(numel(C.symmetricRangeUm) == 2 && ...
        C.symmetricRangeUm(1) < C.symmetricRangeUm(2), ...
        'P.anatomy.centralX.symmetricRangeUm must be [low high].');
    assert(C.minCellsForFraction >= 1 && C.minCellsForPeak >= 1, ...
        'Central-x minimum cell counts must be positive.');
    assert(isfinite(C.kdeBandwidthUm) && C.kdeBandwidthUm > 0, ...
        'P.anatomy.centralX.kdeBandwidthUm must be positive.');
    assert(numel(C.kdeSearchRangeUm) == 2 && ...
        all(isfinite(C.kdeSearchRangeUm)) && ...
        C.kdeSearchRangeUm(1) < C.kdeSearchRangeUm(2), ...
        'P.anatomy.centralX.kdeSearchRangeUm must be [low high].');
    assert(isfinite(C.kdeGridStepUm) && C.kdeGridStepUm > 0, ...
        'P.anatomy.centralX.kdeGridStepUm must be positive.');
    assert(C.bootstrapIterations >= 1 && ...
        C.bootstrapIterations == round(C.bootstrapIterations), ...
        'P.anatomy.centralX.bootstrapIterations must be a positive integer.');
    assert(C.bootstrapCI > 0 && C.bootstrapCI < 100, ...
        'P.anatomy.centralX.bootstrapCI must be between 0 and 100.');
end

function x = finiteClassXUm(Tfish, className)
    idx = string(Tfish.class3) == string(className);
    x = double(Tfish.xUm(idx));
    x = x(isfinite(x));
end

function [peakUm, ciLowUm, ciHighUm, bootPeaks] = ...
        estimateKdePeakWithBootstrap(x, C, randomSeed)
    x = double(x(:));
    x = x(isfinite(x));

    peakUm = NaN;
    ciLowUm = NaN;
    ciHighUm = NaN;
    bootPeaks = nan(C.bootstrapIterations,1);

    if numel(x) < C.minCellsForPeak
        return;
    end

    gridUm = (C.kdeSearchRangeUm(1):C.kdeGridStepUm: ...
        C.kdeSearchRangeUm(2))';
    n = numel(x);

    % Precompute every neuron's Gaussian contribution to every grid point.
    % Bootstrap resamples are then evaluated together with one matrix
    % multiplication, which is much faster than rerunning a KDE in a loop.
    z = (gridUm - x') ./ C.kdeBandwidthUm;
    kernelMatrix = exp(-0.5 .* z.^2) ./ ...
        (sqrt(2*pi) .* C.kdeBandwidthUm);

    density = mean(kernelMatrix,2);
    peakUm = peakFromKdeDensity(gridUm,density);

    previousRng = rng;
    cleaner = onCleanup(@() rng(previousRng)); %#ok<NASGU>
    rng(randomSeed,'twister');

    B = C.bootstrapIterations;
    sampledIndices = randi(n,n,B);
    bootstrapCounts = zeros(n,B);
    for b = 1:B
        bootstrapCounts(:,b) = accumarray( ...
            sampledIndices(:,b),1,[n 1]);
    end

    bootstrapDensity = (kernelMatrix * bootstrapCounts) ./ n;
    [~,peakIndex] = max(bootstrapDensity,[],1);
    bootPeaks = gridUm(peakIndex(:));

    finiteBoot = bootPeaks(isfinite(bootPeaks));
    if ~isempty(finiteBoot)
        ciTail = (100 - C.bootstrapCI) / 2;
        ci = prctile(finiteBoot,[ciTail 100-ciTail]);
        ciLowUm = ci(1);
        ciHighUm = ci(2);
    end
end

function peakUm = peakFromKdeDensity(gridUm,density)
    maxDensity = max(density,[],'omitnan');
    maxIdx = find(density == maxDensity);
    if isempty(maxIdx)
        peakUm = NaN;
    else
        % Average tied grid maxima to avoid an arbitrary left/right choice.
        peakUm = mean(gridUm(maxIdx),'omitnan');
    end
end

function M = getCentralXMetricDefinitions()
    metricName = [ ...
        "ccwCentralLeftFraction"; ...
        "cwCentralRightFraction"; ...
        "opponentCentralEnrichment"; ...
        "deltaXPeakUm"; ...
        "symmetricCentralFraction"];

    yLabel = [ ...
        "P(-20 <= x < 0 | CCW)"; ...
        "P(0 < x <= 20 | CW)"; ...
        "CCW-left + CW-right fraction"; ...
        "KDE peak x(CW) - peak x(CCW) (micrometres)"; ...
        "P(-10 <= x <= 10 | Symmetric)"];

    panelTitle = [ ...
        "CCW central-left fraction"; ...
        "CW central-right fraction"; ...
        "Opponent central enrichment"; ...
        "Directional peak separation"; ...
        "Symmetric central fraction"];

    explanation = [ ...
        "Fraction of each fish's CCW cells in [-20,0) micrometres"; ...
        "Fraction of each fish's CW cells in (0,20] micrometres"; ...
        "Sum of the CCW central-left and CW central-right fractions"; ...
        "Fixed-bandwidth KDE peak(CW) minus peak(CCW); positive is opponent ordering"; ...
        "Fraction of each fish's Symmetric cells in [-10,10] micrometres"];

    referenceLine = [NaN; NaN; NaN; 0; NaN];
    M = table(metricName, yLabel, panelTitle, explanation, referenceLine);
end

function valid = centralMetricValidMask(T, metricName)
    metricName = string(metricName);
    switch metricName
        case "ccwCentralLeftFraction"
            valid = isfinite(T.ccwCentralLeftFraction) & ...
                logical(T.enoughCCWForFraction);
        case "cwCentralRightFraction"
            valid = isfinite(T.cwCentralRightFraction) & ...
                logical(T.enoughCWForFraction);
        case "opponentCentralEnrichment"
            valid = isfinite(T.opponentCentralEnrichment) & ...
                logical(T.enoughCCWForFraction) & ...
                logical(T.enoughCWForFraction);
        case "deltaXPeakUm"
            valid = isfinite(T.deltaXPeakUm) & ...
                isfinite(T.deltaXPeak_CILowUm) & ...
                isfinite(T.deltaXPeak_CIHighUm) & ...
                logical(T.enoughDirectionalForPeaks);
        case "symmetricCentralFraction"
            valid = isfinite(T.symmetricCentralFraction) & ...
                logical(T.enoughSymmetricForFraction);
        otherwise
            error('Unknown central-x metric: %s',metricName);
    end
end

function StatsTable = runAdaptiveCentralXStatistics( ...
        FishCentralXMetrics, CentralXMetricDefinitions, P)

    rows = {};
    morphOrder = string({P.morphs.name});
    pairs = nchoosek(1:numel(morphOrder), 2);

    for mi = 1:height(CentralXMetricDefinitions)
        metricName = CentralXMetricDefinitions.metricName(mi);
        y = double(FishCentralXMetrics.(metricName));
        valid = centralMetricValidMask(FishCentralXMetrics, metricName);

        pRaw = nan(size(pairs,1),1);
        tempRows = cell(size(pairs,1),1);

        for p = 1:size(pairs,1)
            i1 = pairs(p,1);
            i2 = pairs(p,2);
            g1 = morphOrder(i1);
            g2 = morphOrder(i2);

            pairMask = valid & ismember(string(FishCentralXMetrics.morph), ...
                [g1 g2]);
            Tp = FishCentralXMetrics(pairMask,:);
            yp = double(Tp.(metricName));
            group = string(Tp.morph) == g2;

            switch metricName
                case "ccwCentralLeftFraction"
                    success = double(Tp.nCCWCentralLeft);
                    total = double(Tp.nCCW);
                    [pValue,effect,testStatistic,dispersion,info1,info2] = ...
                        quasiBinomialFishTest(success,total,group, ...
                        P.stats.minFishPerMorph);
                    testMethod = "quasi-binomial GLM";
                    precisionBasis = "CCW neurons per fish";
                    effectType = "odds ratio (group2/group1)";

                case "cwCentralRightFraction"
                    success = double(Tp.nCWCentralRight);
                    total = double(Tp.nCW);
                    [pValue,effect,testStatistic,dispersion,info1,info2] = ...
                        quasiBinomialFishTest(success,total,group, ...
                        P.stats.minFishPerMorph);
                    testMethod = "quasi-binomial GLM";
                    precisionBasis = "CW neurons per fish";
                    effectType = "odds ratio (group2/group1)";

                case "symmetricCentralFraction"
                    success = double(Tp.nSymmetricCentral);
                    total = double(Tp.nSymmetric);
                    [pValue,effect,testStatistic,dispersion,info1,info2] = ...
                        quasiBinomialFishTest(success,total,group, ...
                        P.stats.minFishPerMorph);
                    testMethod = "quasi-binomial GLM";
                    precisionBasis = "Symmetric neurons per fish";
                    effectType = "odds ratio (group2/group1)";

                case "opponentCentralEnrichment"
                    weights = centralOpponentPrecisionWeights(Tp,P);
                    [pValue,effect,testStatistic,info1,info2] = ...
                        weightedPermutationTest(yp,weights,group,P, ...
                        P.stats.permutationRandomSeed + 1000*mi + 10*p);
                    dispersion = NaN;
                    testMethod = "weighted fish-level permutation";
                    precisionBasis = "inverse binomial variance of both fractions";
                    effectType = "weighted mean difference (group2-group1)";

                case "deltaXPeakUm"
                    weights = peakPrecisionWeights(Tp,P);
                    [pValue,effect,testStatistic,info1,info2] = ...
                        weightedPermutationTest(yp,weights,group,P, ...
                        P.stats.permutationRandomSeed + 1000*mi + 10*p);
                    dispersion = NaN;
                    testMethod = "weighted fish-level permutation";
                    precisionBasis = "inverse variance from bootstrap CI";
                    effectType = "weighted mean difference (group2-group1)";

                otherwise
                    error('No adaptive test assigned to %s.',metricName);
            end
            pRaw(p) = pValue;

            v1 = yp(~group & isfinite(yp));
            v2 = yp(group & isfinite(yp));

            tempRows{p} = table(metricName, g1, g2, i1, i2, ...
                numel(v1), numel(v2), medianOrNan(v1), medianOrNan(v2), ...
                meanOrNan(v1), meanOrNan(v2), testMethod,precisionBasis, ...
                effect,effectType,testStatistic,dispersion,info1,info2, ...
                pValue,NaN, ...
                'VariableNames', {'metric','group1','group2', ...
                'group1Index','group2Index','n1','n2','median1','median2', ...
                'mean1','mean2','testMethod','precisionBasis', ...
                'effectEstimate','effectType','testStatistic','dispersion', ...
                'information1','information2','pRaw','pAdjusted'});
        end

        for p = 1:numel(tempRows)
            % Retain the legacy column name for output compatibility.
            tempRows{p}.pAdjusted = pRaw(p);
            rows{end+1,1} = tempRows{p}; %#ok<AGROW>
        end
    end

    StatsTable = vertcatTables(rows);
end

function addBootstrapCIToFishPoints(ax, xPosition, values, ciLow, ciHigh, P)
    values = double(values(:));
    ciLow = double(ciLow(:));
    ciHigh = double(ciHigh(:));
    good = isfinite(values) & isfinite(ciLow) & isfinite(ciHigh) & ...
        ciLow <= ciHigh;
    if ~any(good); return; end

    % Reproduce the exact jitter generated by drawBoxAndFishPoints.
    previousRng = rng;
    cleaner = onCleanup(@() rng(previousRng)); %#ok<NASGU>
    rng(P.randomSeed + xPosition);
    jitter = 0.13 * (rand(numel(values),1) - 0.5);
    x = xPosition + jitter;

    capHalfWidth = 0.025;
    for i = find(good(:))'
        plot(ax,[x(i) x(i)],[ciLow(i) ciHigh(i)],'k-', ...
            'LineWidth',0.9,'HandleVisibility','off');
        plot(ax,x(i)+[-capHalfWidth capHalfWidth],[ciLow(i) ciLow(i)], ...
            'k-','LineWidth',0.9,'HandleVisibility','off');
        plot(ax,x(i)+[-capHalfWidth capHalfWidth],[ciHigh(i) ciHigh(i)], ...
            'k-','LineWidth',0.9,'HandleVisibility','off');
    end
end

function plotHorizontalCI(ax, estimate, ciLow, ciHigh, y, color)
    if ~all(isfinite([estimate ciLow ciHigh]))
        return;
    end
    plot(ax,[ciLow ciHigh],[y y],'-','Color',color,'LineWidth',2, ...
        'HandleVisibility','off');
    plot(ax,[ciLow ciLow],[y-0.08 y+0.08],'-','Color',color, ...
        'LineWidth',1,'HandleVisibility','off');
    plot(ax,[ciHigh ciHigh],[y-0.08 y+0.08],'-','Color',color, ...
        'LineWidth',1,'HandleVisibility','off');
    scatter(ax,estimate,y,38,color,'filled','MarkerEdgeColor','k', ...
        'LineWidth',0.5,'HandleVisibility','off');
end

function M = getOpponentMetricDefinitions()
    metricName = [ ...
        "opponentLR"; ...
        "deltaMedianX_CWminusCCW"; ...
        "signedDensitySegregation"];

    yLabel = [ ...
        "Opponent left-right score"; ...
        "Median x(CW) - median x(CCW), normalized units"; ...
        "Signed CW/CCW density segregation"];

    panelTitle = [ ...
        "Side-based opponent score"; ...
        "Continuous CW-CCW displacement"; ...
        "Full-density signed segregation"];

    explanation = [ ...
        "+1: CW entirely right and CCW entirely left; 0: no net opposition; -1: reversed"; ...
        "Positive: CW population is to the right of CCW"; ...
        "Magnitude = 1 - density overlap; positive sign = CW right of CCW"];

    referenceLine = [0; 0; 0];

    M = table(metricName, yLabel, panelTitle, explanation, referenceLine);
end

function StatsTable = runAdaptiveOpponentStatistics( ...
        FishOpponentMetrics, OpponentMetricDefinitions, P)

    rows = {};
    morphOrder = string({P.morphs.name});
    pairs = nchoosek(1:numel(morphOrder), 2);

    for mi = 1:height(OpponentMetricDefinitions)
        metricName = OpponentMetricDefinitions.metricName(mi);
        y = double(FishOpponentMetrics.(metricName));
        valid = isfinite(y) & logical(FishOpponentMetrics.enoughCells);

        pRaw = nan(size(pairs,1), 1);
        tempRows = cell(size(pairs,1), 1);

        for p = 1:size(pairs,1)
            i1 = pairs(p,1);
            i2 = pairs(p,2);
            g1 = morphOrder(i1);
            g2 = morphOrder(i2);

            pairMask = valid & ismember(string(FishOpponentMetrics.morph), ...
                [g1 g2]);
            Tp = FishOpponentMetrics(pairMask,:);
            yp = double(Tp.(metricName));
            group = string(Tp.morph) == g2;
            effectiveCount = 2 .* double(Tp.nCW) .* double(Tp.nCCW) ./ ...
                (double(Tp.nCW) + double(Tp.nCCW));
            weights = normalizePrecisionWeights(effectiveCount,P);
            [pValue,effect,testStatistic,info1,info2] = ...
                weightedPermutationTest(yp,weights,group,P, ...
                P.stats.permutationRandomSeed + 10000 + 1000*mi + 10*p);
            pRaw(p) = pValue;

            v1 = yp(~group & isfinite(yp));
            v2 = yp(group & isfinite(yp));
            testMethod = "weighted fish-level permutation";
            precisionBasis = "harmonic mean of nCW and nCCW";
            effectType = "weighted mean difference (group2-group1)";

            tempRows{p} = table(metricName, g1, g2, i1, i2, ...
                numel(v1), numel(v2), medianOrNan(v1), medianOrNan(v2), ...
                meanOrNan(v1), meanOrNan(v2), testMethod,precisionBasis, ...
                effect,effectType,testStatistic,NaN,info1,info2, ...
                pValue,NaN, ...
                'VariableNames', {'metric','group1','group2', ...
                'group1Index','group2Index','n1','n2','median1','median2', ...
                'mean1','mean2','testMethod','precisionBasis', ...
                'effectEstimate','effectType','testStatistic','dispersion', ...
                'information1','information2','pRaw','pAdjusted'});
        end

        for p = 1:numel(tempRows)
            % Retain the legacy column name for output compatibility.
            tempRows{p}.pAdjusted = pRaw(p);
            rows{end+1,1} = tempRows{p}; %#ok<AGROW>
        end
    end

    StatsTable = vertcatTables(rows);
end

function M = getMetricDefinitions()
    metricName = [ ...
        "nCells"; ...
        "pctOfClassifiedCells"; ...
        "densityPerReferenceHullArea"; ...
        "localDensityIndex"; ...
        "hullAreaPerCellNorm"; ...
        "medianXNorm"];

    yLabel = [ ...
        "Number of cells (test offset: all imaged neurons)"; ...
        "Population size (% of classified cells)"; ...
        "Cells per sampled-region hull area"; ...
        "Local density index, 1 / (pi x median NN distance^2)"; ...
        "Convex-hull area per cell (lower = more compact)"; ...
        "Median normalized x position"];

    fileBase = [ ...
        "01_number_of_cells"; ...
        "02_population_fraction"; ...
        "03_spatial_density"; ...
        "04_local_nearest_neighbour_density"; ...
        "06_compactness_hull_area_per_cell"; ...
        "07_position_along_x"];

    minCells = [0; 0; 0; 2; 3; 1];
    referenceLine = [NaN; NaN; NaN; NaN; NaN; 0];

    M = table(metricName, yLabel, fileBase, minCells, referenceLine);
end

function StatsTable = runAdaptiveClassStatistics( ...
        FishClassMetrics,MetricDefinitions,P)
    rows = {};
    morphOrder = string({P.morphs.name});
    classOrder = string(P.classes);
    pairs = nchoosek(1:numel(morphOrder), 2);

    for mi = 1:height(MetricDefinitions)
        metricName = MetricDefinitions.metricName(mi);
        minCells = MetricDefinitions.minCells(mi);

        for c = 1:numel(classOrder)
            className = classOrder(c);
            T = FishClassMetrics(string(FishClassMetrics.className) == className, :);
            y = double(T.(metricName));
            valid = isfinite(y) & T.nCells >= minCells;
            if metricName == "nCells"
                valid = valid & isfinite(T.nTotalImagedNeurons) & ...
                    T.nTotalImagedNeurons > 0;
            elseif metricName == "pctOfClassifiedCells"
                valid = valid & isfinite(T.nAllClassifiedCells) & ...
                    T.nAllClassifiedCells > 0;
            elseif metricName == "densityPerReferenceHullArea"
                valid = valid & isfinite(T.referenceHullAreaNorm) & ...
                    T.referenceHullAreaNorm > 0;
            end

            pRaw = nan(size(pairs,1), 1);
            tempRows = cell(size(pairs,1), 1);
            for p = 1:size(pairs,1)
                i1 = pairs(p,1);
                i2 = pairs(p,2);
                g1 = morphOrder(i1);
                g2 = morphOrder(i2);
                pairMask = valid & ismember(string(T.morph),[g1 g2]);
                Tp = T(pairMask,:);
                yp = double(Tp.(metricName));
                group = string(Tp.morph) == g2;

                switch metricName
                    case "nCells"
                        count = double(Tp.nCells);
                        exposure = double(Tp.nTotalImagedNeurons);
                        glmValid = isfinite(count) & count >= 0 & ...
                            isfinite(exposure) & exposure > 0;
                        [pValue,effect,testStatistic,dispersion,info1,info2] = ...
                            quasiPoissonRateTest(count(glmValid), ...
                            exposure(glmValid),group(glmValid), ...
                            P.stats.minFishPerMorph);
                        testMethod = "quasi-Poisson GLM with offset";
                        precisionBasis = "total imaged neurons from original RASTER";
                        effectType = "rate ratio per imaged neuron (group2/group1)";

                    case "pctOfClassifiedCells"
                        success = double(Tp.nCells);
                        total = double(Tp.nAllClassifiedCells);
                        [pValue,effect,testStatistic,dispersion,info1,info2] = ...
                            quasiBinomialFishTest(success,total,group, ...
                            P.stats.minFishPerMorph);
                        testMethod = "quasi-binomial GLM";
                        precisionBasis = "number of classified neurons per fish";
                        effectType = "odds ratio (group2/group1)";

                    case "densityPerReferenceHullArea"
                        count = double(Tp.nCells);
                        exposure = double(Tp.referenceHullAreaNorm);
                        glmValid = isfinite(count) & count >= 0 & ...
                            isfinite(exposure) & exposure > 0;
                        [pValue,effect,testStatistic,dispersion,info1,info2] = ...
                            quasiPoissonRateTest(count(glmValid), ...
                            exposure(glmValid),group(glmValid), ...
                            P.stats.minFishPerMorph);
                        testMethod = "quasi-Poisson GLM with area offset";
                        precisionBasis = "sampled reference-hull area and cell count";
                        effectType = "density rate ratio (group2/group1)";

                    otherwise
                        weights = normalizePrecisionWeights( ...
                            double(Tp.nCells),P);
                        [pValue,effect,testStatistic,info1,info2] = ...
                            weightedPermutationTest(yp,weights,group,P, ...
                            P.stats.permutationRandomSeed + ...
                            100000 + 10000*mi + 100*c + p);
                        dispersion = NaN;
                        testMethod = "weighted fish-level permutation";
                        precisionBasis = "number of contributing neurons per fish";
                        effectType = "weighted mean difference (group2-group1)";
                end
                pRaw(p) = pValue;

                v1 = yp(~group & isfinite(yp));
                v2 = yp(group & isfinite(yp));

                tempRows{p} = table(metricName, className, g1, g2, i1, i2, ...
                    numel(v1), numel(v2), medianOrNan(v1), medianOrNan(v2), ...
                    meanOrNan(v1), meanOrNan(v2),testMethod,precisionBasis, ...
                    effect,effectType,testStatistic,dispersion,info1,info2, ...
                    pValue,NaN, ...
                    'VariableNames', {'metric','className','group1','group2', ...
                    'group1Index','group2Index','n1','n2','median1','median2', ...
                    'mean1','mean2','testMethod','precisionBasis', ...
                    'effectEstimate','effectType','testStatistic','dispersion', ...
                    'information1','information2','pRaw','pAdjusted'});
            end

            for p = 1:numel(tempRows)
                % Retain the legacy column name for output compatibility.
                tempRows{p}.pAdjusted = pRaw(p);
                rows{end+1,1} = tempRows{p}; %#ok<AGROW>
            end
        end
    end

    StatsTable = vertcatTables(rows);
end

function [pValue,effect,testStatistic,dispersion,info1,info2] = ...
        quasiPoissonRateTest(count,exposure,group,minN)
    pValue = NaN;
    effect = NaN;
    testStatistic = NaN;
    dispersion = NaN;
    info1 = NaN;
    info2 = NaN;

    count = double(count(:));
    exposure = double(exposure(:));
    group = logical(group(:));
    valid = isfinite(count) & count >= 0 & isfinite(exposure) & exposure > 0;
    count = count(valid);
    exposure = exposure(valid);
    group = group(valid);
    if sum(~group) < minN || sum(group) < minN; return; end

    info1 = sum(exposure(~group));
    info2 = sum(exposure(group));
    try
        [b,~,stats] = glmfit(double(group),count,'poisson', ...
            'link','log','offset',log(exposure),'estdisp','on');
        effect = exp(b(2));
        testStatistic = b(2) / stats.se(2);
        pValue = stats.p(2);
        if isfield(stats,'s') && isfinite(stats.s)
            dispersion = stats.s.^2;
        end
    catch ME
        warning('Quasi-Poisson GLM failed: %s',ME.message);
    end
end

function [pValue,effect,testStatistic,dispersion,info1,info2] = ...
        quasiBinomialFishTest(success,total,group,minN)
    pValue = NaN;
    effect = NaN;
    testStatistic = NaN;
    dispersion = NaN;
    info1 = NaN;
    info2 = NaN;

    success = double(success(:));
    total = double(total(:));
    group = logical(group(:));
    valid = isfinite(success) & success >= 0 & isfinite(total) & ...
        total > 0 & success <= total;
    success = success(valid);
    total = total(valid);
    group = group(valid);
    if sum(~group) < minN || sum(group) < minN; return; end

    info1 = sum(total(~group));
    info2 = sum(total(group));
    try
        [b,~,stats] = glmfit(double(group),[success total], ...
            'binomial','link','logit','estdisp','on');
        effect = exp(b(2));
        testStatistic = b(2) / stats.se(2);
        pValue = stats.p(2);
        if isfield(stats,'s') && isfinite(stats.s)
            dispersion = stats.s.^2;
        end
    catch ME
        warning('Quasi-binomial GLM failed: %s',ME.message);
    end
end

function weights = normalizePrecisionWeights(rawPrecision,P)
    weights = double(rawPrecision(:));
    valid = isfinite(weights) & weights > 0;
    weights(~valid) = NaN;
    if ~any(valid); return; end

    scale = median(weights(valid));
    if ~isfinite(scale) || scale <= 0; scale = mean(weights(valid)); end
    weights(valid) = weights(valid) ./ scale;
    weights(valid) = min(weights(valid),P.stats.maxRelativeWeight);
end

function weights = centralOpponentPrecisionWeights(T,P)
    nCCW = double(T.nCCW);
    nCW = double(T.nCW);
    % Jeffreys-smoothed proportions prevent zero estimated variance when a
    % fish has all or none of its cells inside the central interval.
    pCCW = (double(T.nCCWCentralLeft) + 0.5) ./ (nCCW + 1);
    pCW = (double(T.nCWCentralRight) + 0.5) ./ (nCW + 1);
    variance = pCCW .* (1-pCCW) ./ max(nCCW,1) + ...
        pCW .* (1-pCW) ./ max(nCW,1);
    weights = normalizePrecisionWeights(1 ./ variance,P);
end

function weights = peakPrecisionWeights(T,P)
    ciLow = double(T.deltaXPeak_CILowUm);
    ciHigh = double(T.deltaXPeak_CIHighUm);
    ciLevel = double(T.bootstrapCI);
    pUpper = (1 + ciLevel./100) ./ 2;
    z = -sqrt(2) .* erfcinv(2 .* pUpper);
    se = (ciHigh-ciLow) ./ (2 .* z);
    se(~isfinite(se)) = NaN;
    finiteNonpositive = isfinite(se) & se <= 0;
    se(finiteNonpositive) = P.stats.peakSEFloorUm;
    se(isfinite(se)) = max(se(isfinite(se)),P.stats.peakSEFloorUm);
    weights = normalizePrecisionWeights(1 ./ (se.^2),P);
end

function [pValue,effect,testStatistic,info1,info2] = ...
        weightedPermutationTest(y,weights,group,P,randomSeed)
    pValue = NaN;
    effect = NaN;
    testStatistic = NaN;
    info1 = NaN;
    info2 = NaN;

    y = double(y(:));
    weights = double(weights(:));
    group = logical(group(:));
    valid = isfinite(y) & isfinite(weights) & weights > 0;
    y = y(valid);
    weights = weights(valid);
    group = group(valid);
    n1 = sum(~group);
    n2 = sum(group);
    if n1 < P.stats.minFishPerMorph || n2 < P.stats.minFishPerMorph
        return;
    end

    effect = weightedGroupDifference(y,weights,group);
    testStatistic = effect;
    info1 = kishEffectiveSampleSize(weights(~group));
    info2 = kishEffectiveSampleSize(weights(group));
    observed = abs(effect);
    n = numel(y);

    logCombinations = gammaln(n+1)-gammaln(n1+1)-gammaln(n2+1);
    if logCombinations <= log(P.stats.maxExactPermutations)
        combinations = nchoosek(1:n,n2);
        permuted = nan(size(combinations,1),1);
        for b = 1:size(combinations,1)
            permGroup = false(n,1);
            permGroup(combinations(b,:)) = true;
            permuted(b) = abs(weightedGroupDifference( ...
                y,weights,permGroup));
        end
        pValue = mean(permuted >= observed-10*eps(max(1,observed)));
    else
        previousRng = rng;
        cleaner = onCleanup(@() rng(previousRng)); %#ok<NASGU>
        rng(randomSeed,'twister');
        B = P.stats.permutationIterations;
        exceed = 0;
        for b = 1:B
            order = randperm(n);
            permGroup = false(n,1);
            permGroup(order(1:n2)) = true;
            permuted = abs(weightedGroupDifference(y,weights,permGroup));
            exceed = exceed + (permuted >= observed-10*eps(max(1,observed)));
        end
        pValue = (exceed+1) / (B+1);
    end
end

function d = weightedGroupDifference(y,weights,group)
    mean1 = sum(weights(~group).*y(~group)) / sum(weights(~group));
    mean2 = sum(weights(group).*y(group)) / sum(weights(group));
    d = mean2-mean1;
end

function nEff = kishEffectiveSampleSize(weights)
    weights = double(weights(:));
    nEff = (sum(weights).^2) ./ sum(weights.^2);
end

function plotMetricByClassWithPairwiseBars(FishClassMetrics, PairwiseStats, metricDef, P)
    metricName = metricDef.metricName;
    yLabel = metricDef.yLabel;
    fileBase = metricDef.fileBase;
    minCells = metricDef.minCells;
    referenceLine = metricDef.referenceLine;

    classOrder = string(P.classes);
    morphOrder = string({P.morphs.name});
    morphLabels = string({P.morphs.display});

    fig = figure('Color','w','Visible',P.figureVisible, ...
        'Position',[80 100 1320 470]);
    tl = tiledlayout(1, numel(classOrder), 'TileSpacing','compact', 'Padding','compact');

    for c = 1:numel(classOrder)
        ax = nexttile(tl, c);
        hold(ax, 'on');
        className = classOrder(c);
        T = FishClassMetrics(string(FishClassMetrics.className) == className, :);
        yAll = double(T.(metricName));
        validAll = isfinite(yAll) & T.nCells >= minCells;
        if metricName == "nCells"
            validAll = validAll & isfinite(T.nTotalImagedNeurons) & ...
                T.nTotalImagedNeurons > 0;
        elseif metricName == "pctOfClassifiedCells"
            validAll = validAll & T.nAllClassifiedCells > 0;
        elseif metricName == "densityPerReferenceHullArea"
            validAll = validAll & isfinite(T.referenceHullAreaNorm) & ...
                T.referenceHullAreaNorm > 0;
        end

        for m = 1:numel(morphOrder)
            idx = validAll & string(T.morph) == morphOrder(m);
            values = yAll(idx);
            labels = T.session(idx);
            drawBoxAndFishPoints(ax, m, values, P.morphColors(m,:), labels, P);
        end

        xlim(ax, [0.45 numel(morphOrder)+0.55]);
        xticks(ax, 1:numel(morphOrder));
        xticklabels(ax, morphLabels);
        xtickangle(ax, 0);
        ylabel(ax, yLabel, 'Interpreter','none');
        title(ax, className, 'Interpreter','none');
        set(ax, 'TickDir','out', 'LineWidth',1, 'Box','off', 'FontSize',10);
        grid(ax, 'on');
        ax.XGrid = 'off';

        if isfinite(referenceLine)
            yline(ax, referenceLine, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
        end

        S = PairwiseStats(PairwiseStats.metric == metricName & ...
            PairwiseStats.className == className, :);
        addPairwiseStatBars(ax, S, yAll(validAll), P);
    end

    title(tl, strrep(fileBase, '_', ' '), 'Interpreter','none', 'FontWeight','bold');
    subtitle(tl, ['Distribution-aware fish-level tests; neuron-count ' ...
        'precision modeled; unadjusted p-values']);
    saveFigure(fig, P.figureDir, char(fileBase), P);
end

function drawBoxAndFishPoints(ax, xPosition, values, color, labels, P)
    values = double(values(:));
    values = values(isfinite(values));
    if isempty(values); return; end

    q1 = prctile(values, 25);
    med = median(values);
    q3 = prctile(values, 75);
    iqrValue = q3 - q1;
    lowerFence = q1 - 1.5 * iqrValue;
    upperFence = q3 + 1.5 * iqrValue;
    lowerWhisker = min(values(values >= lowerFence));
    upperWhisker = max(values(values <= upperFence));
    if isempty(lowerWhisker); lowerWhisker = min(values); end
    if isempty(upperWhisker); upperWhisker = max(values); end

    boxWidth = 0.54;
    patch(ax, xPosition + [-boxWidth/2 boxWidth/2 boxWidth/2 -boxWidth/2], ...
        [q1 q1 q3 q3], color, 'FaceAlpha',0.70, ...
        'EdgeColor','k', 'LineWidth',1.1, 'HandleVisibility','off');
    plot(ax, xPosition + [-boxWidth/2 boxWidth/2], [med med], ...
        'k-', 'LineWidth',2, 'HandleVisibility','off');
    plot(ax, [xPosition xPosition], [q3 upperWhisker], ...
        'k-', 'LineWidth',1.1, 'HandleVisibility','off');
    plot(ax, [xPosition xPosition], [q1 lowerWhisker], ...
        'k-', 'LineWidth',1.1, 'HandleVisibility','off');
    capWidth = 0.22;
    plot(ax, xPosition + [-capWidth/2 capWidth/2], [upperWhisker upperWhisker], ...
        'k-', 'LineWidth',1.1, 'HandleVisibility','off');
    plot(ax, xPosition + [-capWidth/2 capWidth/2], [lowerWhisker lowerWhisker], ...
        'k-', 'LineWidth',1.1, 'HandleVisibility','off');

    rng(P.randomSeed + xPosition);
    jitter = 0.13 * (rand(numel(values),1) - 0.5);
    scatter(ax, xPosition + jitter, values, 42, color, 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',0.55, 'MarkerFaceAlpha',0.85, ...
        'HandleVisibility','off');

    if P.addFishLabels && numel(labels) == numel(values)
        for i = 1:numel(values)
            text(ax, xPosition + jitter(i) + 0.035, values(i), string(labels(i)), ...
                'FontSize',7, 'Interpreter','none');
        end
    end
end

function addPairwiseStatBars(ax, S, plottedValues, P)
    if isempty(S) || height(S) == 0 || isempty(plottedValues)
        return;
    end

    if ~P.stats.showNonsignificant
        S = S(isfinite(S.pAdjusted) & S.pAdjusted < P.stats.alpha, :);
    else
        S = S(isfinite(S.pAdjusted), :);
    end
    if isempty(S); return; end

    % Put adjacent comparisons low and the widest comparison highest.
    spans = S.group2Index - S.group1Index;
    [~, order] = sortrows([spans S.group1Index], [1 2]);
    S = S(order,:);

    dataMin = min(plottedValues, [], 'omitnan');
    dataMax = max(plottedValues, [], 'omitnan');
    dataRange = dataMax - dataMin;
    if ~isfinite(dataRange) || dataRange <= 0
        dataRange = max(1, abs(dataMax));
    end

    currentLimits = ylim(ax);
    dataMax = max(dataMax, currentLimits(2));
    base = dataMax + 0.08 * dataRange;
    step = 0.13 * dataRange;
    tick = 0.035 * dataRange;

    for i = 1:height(S)
        x1 = S.group1Index(i);
        x2 = S.group2Index(i);
        y = base + (i-1) * step;
        plot(ax, [x1 x1 x2 x2], [y-tick y y y-tick], ...
            'k-', 'LineWidth',1.35, 'HandleVisibility','off');
        label = significanceLabel(S.pAdjusted(i), P.stats.labelMode);
        text(ax, mean([x1 x2]), y + 0.025*dataRange, label, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontWeight','bold', 'FontSize',11, 'Interpreter','none');
    end

    newTop = base + max(height(S)-1,0)*step + 0.18*dataRange;
    lower = currentLimits(1);
    if ~isfinite(lower); lower = dataMin - 0.08*dataRange; end
    ylim(ax, [lower newTop]);
    addPairwisePValueText(ax,S);
end

function addPairwisePValueText(ax,S)
    pieces = strings(height(S),1);
    for i = 1:height(S)
        pieces(i) = string(S.group1(i)) + "-" + string(S.group2(i)) + ...
            ": p=" + string(sprintf('%.3g',S.pAdjusted(i)));
    end
    text(ax,0.02,0.98,strjoin(cellstr(pieces),newline), ...
        'Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontSize',8,'Interpreter','none', ...
        'BackgroundColor','w','Margin',2,'Clipping','off');
end

function label = significanceLabel(p, mode)
    if ~isfinite(p)
        label = 'n/a';
        return;
    end
    if strcmpi(mode, 'exact')
        if p < 1e-4
            label = 'p < 0.0001';
        else
            label = sprintf('p = %.3g', p);
        end
        return;
    end

    if p < 1e-4
        label = '****';
    elseif p < 1e-3
        label = '***';
    elseif p < 1e-2
        label = '**';
    elseif p < 0.05
        label = '*';
    else
        label = 'n.s.';
    end
end


function plotXAxisDensityCurvesMicrons(AnatomyNeurons, P)
    % Keep only the signed-log display of the full x-density range.
    % The nonlinear transform affects the horizontal display coordinate only;
    % density values remain probability densities per micrometre.

    classOrder = string(P.classes);
    morphOrder = string({P.morphs.name});
    morphLabels = string({P.morphs.display});

    xAll = double(AnatomyNeurons.xUm);
    xAll = xAll(isfinite(xAll));
    if isempty(xAll)
        warning('No finite xUm coordinates; skipping x-density plot.');
        return;
    end

    binWidthUm = P.anatomy.xDensity.fullBinWidthUm;
    sigmaUm = P.anatomy.xDensity.smoothingSigmaUm;
    signedLogScale = P.anatomy.xDensity.signedLogScaleUm;
    assert(isfinite(binWidthUm) && binWidthUm > 0, ...
        'P.anatomy.xDensity.fullBinWidthUm must be positive.');
    assert(isfinite(sigmaUm) && sigmaUm >= 0, ...
        'P.anatomy.xDensity.smoothingSigmaUm must be nonnegative.');
    assert(isfinite(signedLogScale) && signedLogScale > 0, ...
        'P.anatomy.xDensity.signedLogScaleUm must be positive.');

    xMax = ceil(max(abs(xAll)) / binWidthUm) * binWidthUm;
    if ~isfinite(xMax) || xMax <= 0
        xMax = binWidthUm;
    end
    edges = -xMax:binWidthUm:xMax;
    if edges(end) < xMax
        edges(end+1) = xMax;
    end
    centers = edges(1:end-1) + diff(edges)/2;
    displayX = signedLogTransform(centers, signedLogScale);

    fig = figure('Color','w','Visible',P.figureVisible, ...
        'Position',[45 45 1420 390]);
    tl = tiledlayout(1, numel(classOrder), ...
        'TileSpacing','compact', 'Padding','compact');

    for c = 1:numel(classOrder)
        ax = nexttile(tl, c);
        hold(ax, 'on');
        leg = gobjects(numel(morphOrder),1);
        legText = strings(numel(morphOrder),1);

        for m = 1:numel(morphOrder)
            curves = fishDensityCurvesPerMicron(AnatomyNeurons, ...
                morphOrder(m), classOrder(c), 'xUm', edges, sigmaUm);
            if isempty(curves); continue; end

            mu = mean(curves, 1, 'omitnan');
            nAtBin = sum(isfinite(curves), 1);
            se = std(curves, 0, 1, 'omitnan') ./ sqrt(max(nAtBin,1));
            leg(m) = plotMeanSemCurve(ax, displayX, mu, se, ...
                P.morphColors(m,:));
            legText(m) = sprintf('%s (n=%d fish)', ...
                morphLabels(m), size(curves,1));
        end

        xline(ax, 0, 'k:', 'HandleVisibility','off');
        applySignedLogTicks(ax, xMax, signedLogScale);
        xlabel(ax, sprintf(['Signed-log x display (tick labels in um; ' ...
            'scale = %.0f um)'], signedLogScale), 'Interpreter','none');
        ylabel(ax, 'Probability density (1/um)', 'Interpreter','none');
        title(ax, classOrder(c), 'Interpreter','none');
        styleDensityAxis(ax);

        validLeg = isgraphics(leg);
        if c == numel(classOrder) && any(validLeg)
            legend(ax, leg(validLeg), legText(validLeg), ...
                'Location','eastoutside');
        end
    end

    title(tl, 'Anatomical density along the left-right axis', ...
        'FontWeight','bold');
    subtitle(tl, ['Each fish is normalized before morph averaging. ' ...
        'Signed-log x display expands the region around zero.']);

    saveFigure(fig, P.figureDir, ...
        '10_density_along_x_axis_um_full_zoom_signedlog', P);
end

function plotMorphDensitySet(ax, AnatomyNeurons, className, ...
        morphOrder, morphLabels, edges, sigmaUm, P)

    leg = gobjects(numel(morphOrder),1);
    legText = strings(numel(morphOrder),1);

    for m = 1:numel(morphOrder)
        curves = fishDensityCurvesPerMicron(AnatomyNeurons, ...
            morphOrder(m), className, 'xUm', edges, sigmaUm);
        if isempty(curves); continue; end

        centers = edges(1:end-1) + diff(edges)/2;
        mu = mean(curves, 1, 'omitnan');
        nAtBin = sum(isfinite(curves), 1);
        se = std(curves, 0, 1, 'omitnan') ./ sqrt(max(nAtBin,1));
        leg(m) = plotMeanSemCurve(ax, centers, mu, se, P.morphColors(m,:));
        legText(m) = sprintf('%s (n=%d fish)', ...
            morphLabels(m), size(curves,1));
    end

    validLeg = isgraphics(leg);
    if any(validLeg)
        legend(ax, leg(validLeg), legText(validLeg), 'Location','best');
    end
end

function curves = fishDensityCurvesPerMicron(AnatomyNeurons, ...
        morphName, className, axisVar, edges, sigmaUm)

    T = AnatomyNeurons(string(AnatomyNeurons.morph) == morphName & ...
        string(AnatomyNeurons.class3) == className, :);
    curves = [];
    if isempty(T); return; end

    widths = diff(edges);
    if any(~isfinite(widths)) || any(widths <= 0)
        error('Density edges must be strictly increasing.');
    end

    % The bins are uniform in the supplied full and zoom configurations.
    binWidth = median(widths);
    sigmaBins = sigmaUm / binWidth;

    keys = unique(T.sessionKey, 'stable');
    for i = 1:numel(keys)
        values = double(T.(axisVar)(T.sessionKey == keys(i)));
        values = values(isfinite(values));
        if isempty(values); continue; end

        counts = histcounts(values, edges);

        % Normalize by every class neuron in that fish, including neurons
        % outside a zoomed range. Therefore the zoomed curve preserves the
        % absolute probability density rather than renormalizing the center.
        h = counts ./ (numel(values) * binWidth);
        h = smoothVectorGaussian(h, sigmaBins);
        curves(end+1,:) = h; %#ok<AGROW>
    end
end

function u = signedLogTransform(x, scaleUm)
    x = double(x);
    u = sign(x) .* log10(1 + abs(x) ./ scaleUm);
end

function applySignedLogTicks(ax, xMax, scaleUm)
    candidateTicksUm = [-200 -150 -100 -50 -20 -10 0 10 20 50 100 150 200];
    tickValuesUm = candidateTicksUm(abs(candidateTicksUm) <= xMax + eps);
    if ~any(tickValuesUm == 0)
        tickValuesUm = sort([tickValuesUm 0]);
    end

    xticks(ax, signedLogTransform(tickValuesUm, scaleUm));
    xticklabels(ax, string(tickValuesUm));
    xlim(ax, signedLogTransform([-xMax xMax], scaleUm));
end

function styleDensityAxis(ax)
    set(ax, 'TickDir','out', 'LineWidth',1, 'Box','off');
    grid(ax, 'on');
end

function h = plotMeanSemCurve(ax, x, mu, se, color)
    x = double(x(:));
    mu = double(mu(:));
    se = double(se(:));
    good = isfinite(x) & isfinite(mu);
    x = x(good);
    mu = mu(good);
    se = se(good);

    if isempty(x)
        h = gobjects(1);
        return;
    end

    fill(ax, [x; flipud(x)], [mu-se; flipud(mu+se)], color, ...
        'FaceAlpha',0.18, 'EdgeColor','none', 'HandleVisibility','off');
    h = plot(ax, x, mu, 'Color',color, 'LineWidth',2.2);
end

function plotSpatialOverlays(AnatomyNeurons, P)
    morphOrder = string({P.morphs.name});
    morphLabels = string({P.morphs.display});
    classOrder = string(P.classes);

    fig = figure('Color','w','Visible',P.figureVisible, ...
        'Position',[80 80 1250 430]);
    tl = tiledlayout(1, numel(morphOrder), 'TileSpacing','compact', 'Padding','compact');
    legendHandles = gobjects(numel(classOrder),1);

    for m = 1:numel(morphOrder)
        ax = nexttile(tl, m);
        hold(ax, 'on');
        for c = 1:numel(classOrder)
            idx = string(AnatomyNeurons.morph) == morphOrder(m) & ...
                string(AnatomyNeurons.class3) == classOrder(c);
            x = AnatomyNeurons.xNorm(idx);
            y = AnatomyNeurons.yNorm(idx);
            h = scatter(ax, x, y, 16, P.classColors(c,:), 'filled', ...
                'MarkerFaceAlpha',0.35, 'MarkerEdgeAlpha',0.05);
            if m == numel(morphOrder)
                legendHandles(c) = h;
            end
        end
        xline(ax, 0, 'k:', 'HandleVisibility','off');
        yline(ax, 0, 'k:', 'HandleVisibility','off');
        axis(ax, 'equal');
        xlim(ax, [-1.08 1.08]);
        ylim(ax, [-1.08 1.08]);
        set(ax, 'YDir','reverse'); % preserve image-coordinate orientation
        xlabel(ax, 'Normalized x');
        ylabel(ax, 'Normalized y');
        title(ax, morphLabels(m), 'Interpreter','none');
        set(ax, 'TickDir','out', 'LineWidth',1, 'Box','off');
        grid(ax, 'on');
    end

    validLegend = isgraphics(legendHandles);
    if any(validLegend)
        legend(legendHandles(validLegend), classOrder(validLegend), 'Location','eastoutside');
    end
    title(tl, 'Pooled anatomical ROI centroids (descriptive)', 'FontWeight','bold');
    subtitle(tl, 'FOV-corrected model coordinates, isotropically normalized within each recording');
    saveFigure(fig, P.figureDir, '00_spatial_overlay_by_morph_and_class', P);
end

function area = convexHullAreaSafe(x, y, minCells)
    area = NaN;
    x = double(x(:));
    y = double(y(:));
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);
    if numel(x) < minCells; return; end
    try
        k = convhull(x, y);
        area = polyarea(x(k), y(k));
        if ~isfinite(area) || area <= 0
            area = NaN;
        end
    catch
        area = NaN;
    end
end

function d = medianNearestNeighbourDistance(x, y, minCells)
    d = NaN;
    x = double(x(:));
    y = double(y(:));
    good = isfinite(x) & isfinite(y);
    xy = [x(good) y(good)];
    if size(xy,1) < minCells; return; end

    D = squareform(pdist(xy, 'euclidean'));
    D(1:size(D,1)+1:end) = Inf;
    nearest = min(D, [], 2);
    d = median(nearest, 'omitnan');
    if ~isfinite(d) || d <= 0
        d = NaN;
    end
end

function occupancy = axisOccupancyFraction(x, edges)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        occupancy = NaN;
        return;
    end
    counts = histcounts(x, edges);
    occupancy = sum(counts > 0) / numel(counts);
end

function q = iqrFinite(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        q = NaN;
    else
        q = prctile(x,75) - prctile(x,25);
    end
end

function s = spanFinite(x, minCells)
    x = double(x(:));
    x = x(isfinite(x));
    if numel(x) < minCells
        s = NaN;
    else
        s = max(x) - min(x);
    end
end

function y = smoothVectorGaussian(x, sigmaBins)
    x = double(x(:))';
    if sigmaBins <= 0 || numel(x) < 3
        y = x;
        return;
    end
    radius = max(1, ceil(3*sigmaBins));
    t = -radius:radius;
    kernel = exp(-(t.^2) ./ (2*sigmaBins^2));
    kernel = kernel ./ sum(kernel);
    y = conv(x, kernel, 'same');
end

function p = safePercent(num, den)
    p = 100 * safeDivide(num, den);
end

function x = safeDivide(num, den)
    if ~isfinite(den) || den == 0
        x = NaN;
    else
        x = num ./ den;
    end
end

function m = medianOrNan(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x); m = NaN; else; m = median(x); end
end

function m = meanOrNan(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x); m = NaN; else; m = mean(x); end
end

function key = makeSessionKey(x)
    x = string(x);
    key = strings(size(x));
    for i = 1:numel(x)
        value = lower(strtrim(x(i)));
        if strlength(value) == 0 || ismissing(value)
            key(i) = "unknown";
            continue;
        end
        recToken = regexp(char(value), 'rec\d+', 'match', 'once');
        if ~isempty(recToken)
            key(i) = string(recToken);
        else
            [~, base] = fileparts(char(value));
            if isempty(base); base = char(value); end
            key(i) = string(regexprep(lower(base), '[^a-z0-9]+', ''));
        end
    end
end

function saveFigure(fig, outDir, fileBase, P)
    outPath = fullfile(outDir, [fileBase '.' P.figureFormat]);
    try
        exportgraphics(fig, outPath, 'Resolution', P.figureResolution);
    catch
        saveas(fig, outPath);
    end
    fprintf('Saved figure: %s\n', outPath);
    svgPath = fullfile(outDir, [fileBase '.svg']);
    if ~strcmpi(outPath, svgPath)
        try
            exportgraphics(fig, svgPath, 'ContentType', 'vector');
        catch
            print(fig, svgPath, '-dsvg');
        end
        fprintf('Saved figure: %s\n', svgPath);
    end
end

function T = vertcatTables(tables)
    if isempty(tables)
        T = table();
        return;
    end
    tables = tables(:);
    keep = false(numel(tables),1);
    for i = 1:numel(tables)
        keep(i) = istable(tables{i}) && height(tables{i}) > 0;
    end
    tables = tables(keep);
    if isempty(tables)
        T = table();
        return;
    end

    T = tables{1};
    for i = 2:numel(tables)
        T = vertcatCompatible(T, tables{i});
    end
end

function T = vertcatCompatible(A, B)
    vars = union(A.Properties.VariableNames, B.Properties.VariableNames, 'stable');
    for i = 1:numel(vars)
        if ~ismember(vars{i}, A.Properties.VariableNames)
            A.(vars{i}) = missingColumnLike(B.(vars{i}), height(A));
        end
        if ~ismember(vars{i}, B.Properties.VariableNames)
            B.(vars{i}) = missingColumnLike(A.(vars{i}), height(B));
        end
    end
    A = A(:,vars);
    B = B(:,vars);
    T = [A; B];
end

function col = missingColumnLike(example, n)
    if isnumeric(example)
        col = nan(n, size(example,2));
    elseif islogical(example)
        col = false(n, size(example,2));
    elseif isstring(example)
        col = strings(n, size(example,2));
        col(:) = "";
    elseif iscategorical(example)
        col = categorical(repmat({''}, n, 1));
    elseif iscell(example)
        col = cell(n, size(example,2));
    else
        col = cell(n, 1);
    end
end


%% =======================================================================
%% =============== INDEPENDENT ANATOMY-MATCHING AUDIT ====================
%% =======================================================================

function SourceIndex = buildModelSourceIndex(modelFile, morphName)
    % Recover the exact source files used when the morph-specific model was
    % fitted. Result.paths.allCellsPath points to the ORIGINAL ALL_CELLS file.
    SourceIndex = table();

    try
        S = load(modelFile, 'Results');
    catch ME
        warning('Could not read Results from %s: %s', modelFile, ME.message);
        return;
    end
    if ~isfield(S, 'Results') || isempty(S.Results)
        warning('No Results structure found in %s; matching audit unavailable.', modelFile);
        return;
    end

    rows = {};
    R = S.Results;
    for i = 1:numel(R)
        if iscell(R)
            Ri = R{i};
        else
            Ri = R(i);
        end
        if ~isstruct(Ri) || ~isfield(Ri, 'name') || isempty(Ri.name)
            continue;
        end

        session = string(Ri.name);
        allCellsPath = "";
        fovPath = "";
        if isfield(Ri, 'paths') && isstruct(Ri.paths)
            if isfield(Ri.paths, 'allCellsPath') && ~isempty(Ri.paths.allCellsPath)
                allCellsPath = string(Ri.paths.allCellsPath);
            end
            if isfield(Ri.paths, 'fovPath') && ~isempty(Ri.paths.fovPath)
                fovPath = string(Ri.paths.fovPath);
            end
        end

        rows{end+1,1} = table( ...
            string(morphName), session, makeSessionKey(session), ...
            allCellsPath, fovPath, string(modelFile), ...
            'VariableNames', {'morph','session','sessionKey', ...
            'originalAllCellsPath','fovPath','modelFile'}); %#ok<AGROW>
    end

    SourceIndex = vertcatTables(rows);
end

function [NeuronAudit, FishAudit] = ...
        buildIndependentAnatomyMatchingAudit(AllNeuronsFull, SourceIndex, P)

    NeuronAudit = table();
    FishAudit = table();

    if ~P.matchAudit.enabled
        return;
    end
    if isempty(AllNeuronsFull) || height(AllNeuronsFull) == 0
        warning('Matching audit skipped: AllNeuronsFull is empty.');
        return;
    end
    if isempty(SourceIndex) || height(SourceIndex) == 0
        warning(['Matching audit skipped: no source-file index could be recovered ' ...
            'from the model Results structures.']);
        return;
    end

    required = {'morph','session','sessionKey','neuronID','class3', ...
        'x','y','xCenteredUm','yCenteredUm'};
    missing = required(~ismember(required, AllNeuronsFull.Properties.VariableNames));
    if ~isempty(missing)
        warning('Matching audit skipped; model table lacks: %s', strjoin(missing, ', '));
        return;
    end

    neuronRows = {};
    fishRows = {};
    morphOrder = string({P.morphs.name});

    for m = 1:numel(morphOrder)
        morphName = morphOrder(m);
        Tmorph = AllNeuronsFull(string(AllNeuronsFull.morph) == morphName, :);
        keys = unique(string(Tmorph.sessionKey), 'stable');

        for k = 1:numel(keys)
            key = keys(k);
            T = Tmorph(string(Tmorph.sessionKey) == key, :);
            if isempty(T); continue; end

            sessionName = string(T.session(1));
            srcIdx = find(string(SourceIndex.morph) == morphName & ...
                string(SourceIndex.sessionKey) == key, 1, 'first');

            if isempty(srcIdx)
                fishRows{end+1,1} = makeMissingMatchingAuditFishRow( ...
                    morphName, sessionName, key, height(T), ...
                    "missing source path in model Results"); %#ok<AGROW>
                continue;
            end

            allCellsPath = string(SourceIndex.originalAllCellsPath(srcIdx));
            fovPath = string(SourceIndex.fovPath(srcIdx));

            if strlength(allCellsPath) == 0 || ~isfile(allCellsPath)
                fishRows{end+1,1} = makeMissingMatchingAuditFishRow( ...
                    morphName, sessionName, key, height(T), ...
                    "original ALL_CELLS file missing"); %#ok<AGROW>
                continue;
            end

            ids = double(T.neuronID(:));
            try
                [xRePx, yRePx, validID, pxXFile, pxYFile] = ...
                    extractOriginalRoiCentroidsForAudit(char(allCellsPath), ids);
            catch ME
                warning('Matching audit failed for %s: %s', sessionName, ME.message);
                fishRows{end+1,1} = makeMissingMatchingAuditFishRow( ...
                    morphName, sessionName, key, height(T), ...
                    "could not read original ALL_CELLS"); %#ok<AGROW>
                continue;
            end

            % Use the full, unfiltered model table for the center, matching
            % the original morph-specific model script.
            xCenterPx = median(xRePx, 'omitnan');
            yCenterPx = median(yRePx, 'omitnan');

            pxXStored = medianNumericAudit(T, 'pixelLengthX', 1);
            pxYStored = medianNumericAudit(T, 'pixelLengthY', 1);
            if ~isfinite(pxXFile); pxXFile = pxXStored; end
            if ~isfinite(pxYFile); pxYFile = pxYStored; end

            angleStored = medianNumericAudit(T, 'fov_angle_deg', NaN);
            angleFile = readFovAngleForAudit(char(fovPath));
            angleUsed = angleFile;
            if ~isfinite(angleUsed)
                angleUsed = angleStored;
            end

            xRawRe = (xRePx - xCenterPx) .* pxXFile;
            yRawRe = -(yRePx - yCenterPx) .* pxYFile;

            if isfinite(angleUsed)
                theta = deg2rad(angleUsed);
                Rfov = [cos(theta) sin(theta); -sin(theta) cos(theta)];
                xyImage = [(xRePx - xCenterPx) .* pxXFile, ...
                           (yRePx - yCenterPx) .* pxYFile];
                xyFov = xyImage * Rfov';
                xFovRe = xyFov(:,1);
                yFovRe = -xyFov(:,2);
            else
                xFovRe = nan(height(T),1);
                yFovRe = nan(height(T),1);
            end

            xSavedPx = double(T.x(:));
            ySavedPx = double(T.y(:));
            xSavedFov = double(T.xCenteredUm(:));
            ySavedFov = double(T.yCenteredUm(:));
            xSavedRaw = numericColumnAudit(T, 'xCenteredRawUm');
            ySavedRaw = numericColumnAudit(T, 'yCenteredRawUm');

            usedInAnalysis = ismember(string(T.class3), string(P.classes));
            pixelError = hypot(xRePx - xSavedPx, yRePx - ySavedPx);
            rawError = hypot(xRawRe - xSavedRaw, yRawRe - ySavedRaw);
            fovError = hypot(xFovRe - xSavedFov, yFovRe - ySavedFov);

            N = table();
            N.morph = repmat(morphName, height(T), 1);
            N.session = repmat(sessionName, height(T), 1);
            N.sessionKey = repmat(key, height(T), 1);
            N.neuronID = ids;
            N.class3 = string(T.class3);
            N.usedInCrossMorphAnalysis = usedInAnalysis;
            N.validOriginalRoiID = validID;
            N.savedPixelX = xSavedPx;
            N.savedPixelY = ySavedPx;
            N.recomputedPixelX = xRePx;
            N.recomputedPixelY = yRePx;
            N.pixelCoordinateErrorPx = pixelError;
            N.savedRawXUm = xSavedRaw;
            N.savedRawYUm = ySavedRaw;
            N.recomputedRawXUm = xRawRe;
            N.recomputedRawYUm = yRawRe;
            N.rawCoordinateErrorUm = rawError;
            N.savedFovXUm = xSavedFov;
            N.savedFovYUm = ySavedFov;
            N.recomputedFovXUm = xFovRe;
            N.recomputedFovYUm = yFovRe;
            N.fovCoordinateErrorUm = fovError;
            N.originalAllCellsPath = repmat(allCellsPath, height(T), 1);
            N.fovPath = repmat(fovPath, height(T), 1);
            neuronRows{end+1,1} = N; %#ok<AGROW>

            validFraction = mean(validID);
            nDuplicateIDs = numel(ids(isfinite(ids))) - numel(unique(ids(isfinite(ids))));
            rPixelX = safeCorrelation(xSavedPx, xRePx, ...
                P.matchAudit.minNeuronsForCorrelation);
            rPixelY = safeCorrelation(ySavedPx, yRePx, ...
                P.matchAudit.minNeuronsForCorrelation);
            rFovX = safeCorrelation(xSavedFov, xFovRe, ...
                P.matchAudit.minNeuronsForCorrelation);
            rFovY = safeCorrelation(ySavedFov, yFovRe, ...
                P.matchAudit.minNeuronsForCorrelation);

            medianPixelError = median(pixelError, 'omitnan');
            maxPixelError = max(pixelError, [], 'omitnan');
            medianFovError = median(fovError, 'omitnan');
            maxFovError = max(fovError, [], 'omitnan');

            angleDifference = circularAngleDifferenceDegAudit(angleFile, angleStored);
            pxXDifference = pxXFile - pxXStored;
            pxYDifference = pxYFile - pxYStored;

            status = diagnoseIndependentMatchingAudit( ...
                validFraction, nDuplicateIDs, rPixelX, rPixelY, ...
                medianPixelError, rFovX, rFovY, medianFovError, P.matchAudit);

            fishRows{end+1,1} = table( ...
                morphName, sessionName, key, height(T), ...
                sum(usedInAnalysis), sum(validID), validFraction, nDuplicateIDs, ...
                rPixelX, rPixelY, medianPixelError, maxPixelError, ...
                rFovX, rFovY, medianFovError, maxFovError, ...
                pxXStored, pxYStored, pxXFile, pxYFile, ...
                pxXDifference, pxYDifference, ...
                angleStored, angleFile, angleUsed, angleDifference, ...
                allCellsPath, fovPath, status, ...
                'VariableNames', {'morph','session','sessionKey', ...
                'nModelRows','nRowsUsedInAnalysis','nValidOriginalRoiIDs', ...
                'validOriginalRoiFraction','nDuplicateNeuronIDs', ...
                'rSavedVsRecomputedPixelX','rSavedVsRecomputedPixelY', ...
                'medianPixelCoordinateErrorPx','maxPixelCoordinateErrorPx', ...
                'rSavedVsRecomputedFovX','rSavedVsRecomputedFovY', ...
                'medianFovCoordinateErrorUm','maxFovCoordinateErrorUm', ...
                'pixelLengthXStored','pixelLengthYStored', ...
                'pixelLengthXFromFile','pixelLengthYFromFile', ...
                'pixelLengthXDifference','pixelLengthYDifference', ...
                'fovAngleStoredDeg','fovAngleFileDeg','fovAngleUsedDeg', ...
                'fovAngleDifferenceDeg','originalAllCellsPath','fovPath', ...
                'status'}); %#ok<AGROW>
        end
    end

    NeuronAudit = vertcatTables(neuronRows);
    FishAudit = vertcatTables(fishRows);
    if ~isempty(FishAudit)
        FishAudit.morph = categorical(FishAudit.morph, morphOrder, 'Ordinal', true);
    end
end

function [x, y, validID, pixelLengthX, pixelLengthY] = ...
        extractOriginalRoiCentroidsForAudit(allCellsPath, neuronIDs)

    S = load(allCellsPath);
    neuronIDs = double(neuronIDs(:));

    if isfield(S, 'cell_per')
        nAll = numel(S.cell_per);
    elseif isfield(S, 'cells')
        nAll = numel(S.cells);
    else
        error('Original ALL_CELLS file contains neither cell_per nor cells.');
    end

    validID = isfinite(neuronIDs) & neuronIDs == round(neuronIDs) & ...
        neuronIDs >= 1 & neuronIDs <= nAll;
    x = nan(size(neuronIDs));
    y = nan(size(neuronIDs));

    [imageHeight, imageWidth] = inferImageSize(S);

    for i = find(validID(:))'
        id = neuronIDs(i);
        if isfield(S, 'cell_per') && ~isempty(S.cell_per{id})
            per = double(S.cell_per{id});
            if size(per,2) >= 2
                x(i) = mean(per(:,1), 'omitnan');
                y(i) = mean(per(:,2), 'omitnan');
            elseif size(per,1) >= 2
                x(i) = mean(per(1,:), 'omitnan');
                y(i) = mean(per(2,:), 'omitnan');
            end
        elseif isfield(S, 'cells') && ~isempty(S.cells{id}) && ...
                isfinite(imageHeight) && isfinite(imageWidth)
            pix = double(S.cells{id}(:));
            pix = pix(isfinite(pix) & pix >= 1 & pix <= imageHeight*imageWidth);
            if ~isempty(pix)
                [yy, xx] = ind2sub([imageHeight imageWidth], pix);
                x(i) = mean(xx, 'omitnan');
                y(i) = mean(yy, 'omitnan');
            end
        end
    end

    pixelLengthX = scalarFieldAudit(S, 'pixelLengthX', NaN);
    pixelLengthY = scalarFieldAudit(S, 'pixelLengthY', NaN);
end

function value = scalarFieldAudit(S, fieldName, defaultValue)
    value = defaultValue;
    if isfield(S, fieldName)
        tmp = double(S.(fieldName));
        tmp = tmp(isfinite(tmp));
        if ~isempty(tmp)
            value = tmp(1);
        end
    end
end

function angle = readFovAngleForAudit(fovPath)
    angle = NaN;
    if isempty(fovPath) || ~isfile(fovPath)
        return;
    end
    try
        F = load(fovPath);
        if isfield(F, 'fov_correction_saved')
            F = F.fov_correction_saved;
        end
        if isfield(F, 'fov_angle_deg')
            tmp = double(F.fov_angle_deg);
            tmp = tmp(isfinite(tmp));
            if ~isempty(tmp)
                angle = tmp(1);
            end
        end
    catch
        angle = NaN;
    end
end

function x = numericColumnAudit(T, variableName)
    if ismember(variableName, T.Properties.VariableNames)
        x = double(T.(variableName)(:));
    else
        x = nan(height(T),1);
    end
end

function value = medianNumericAudit(T, variableName, defaultValue)
    value = defaultValue;
    if ismember(variableName, T.Properties.VariableNames)
        x = double(T.(variableName)(:));
        x = x(isfinite(x));
        if ~isempty(x)
            value = median(x);
        end
    end
end

function d = circularAngleDifferenceDegAudit(a, b)
    if ~isfinite(a) || ~isfinite(b)
        d = NaN;
    else
        d = mod((a-b)+180, 360) - 180;
    end
end

function status = diagnoseIndependentMatchingAudit( ...
        validFraction, nDuplicateIDs, rPixelX, rPixelY, medianPixelError, ...
        rFovX, rFovY, medianFovError, D)

    if ~isfinite(validFraction)
        status = "FAIL: no valid original ROI IDs";
    elseif validFraction < D.minValidOriginalRoiFraction
        status = "FAIL: some neuronID values are invalid for original ALL_CELLS";
    elseif nDuplicateIDs > 0
        status = "FAIL: duplicate neuronID values";
    elseif ~isfinite(rPixelX) || ~isfinite(rPixelY)
        status = "CHECK: too few neurons or constant pixel coordinates";
    elseif rPixelX < D.goodCorrelation || rPixelY < D.goodCorrelation
        status = "FAIL: neuronID does not reproduce saved ROI centroids";
    elseif ~isfinite(medianPixelError) || ...
            medianPixelError > D.maxMedianPixelErrorPx
        status = "FAIL: saved and recomputed pixel centroids differ";
    elseif ~isfinite(rFovX) || ~isfinite(rFovY)
        status = "CHECK: FOV-coordinate correlation unavailable";
    elseif rFovX < D.goodCorrelation || rFovY < D.goodCorrelation
        status = "FAIL: reconstructed FOV coordinates do not match saved coordinates";
    elseif ~isfinite(medianFovError) || ...
            medianFovError > D.maxMedianFovErrorUm
        status = "FAIL: reconstructed FOV positions differ from saved positions";
    else
        status = "PASS: neuron IDs and anatomical coordinates match";
    end
end

function R = makeMissingMatchingAuditFishRow(morph, session, key, nRows, status)
    R = table( ...
        string(morph), string(session), string(key), nRows, 0, 0, ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        "", "", string(status), ...
        'VariableNames', {'morph','session','sessionKey', ...
        'nModelRows','nRowsUsedInAnalysis','nValidOriginalRoiIDs', ...
        'validOriginalRoiFraction','nDuplicateNeuronIDs', ...
        'rSavedVsRecomputedPixelX','rSavedVsRecomputedPixelY', ...
        'medianPixelCoordinateErrorPx','maxPixelCoordinateErrorPx', ...
        'rSavedVsRecomputedFovX','rSavedVsRecomputedFovY', ...
        'medianFovCoordinateErrorUm','maxFovCoordinateErrorUm', ...
        'pixelLengthXStored','pixelLengthYStored', ...
        'pixelLengthXFromFile','pixelLengthYFromFile', ...
        'pixelLengthXDifference','pixelLengthYDifference', ...
        'fovAngleStoredDeg','fovAngleFileDeg','fovAngleUsedDeg', ...
        'fovAngleDifferenceDeg','originalAllCellsPath','fovPath','status'});
end

function scatterAuditByMorph(ax, x, y, morph, P)
    morphOrder = string({P.morphs.name});
    morphLabels = string({P.morphs.display});
    for m = 1:numel(morphOrder)
        idx = morph == morphOrder(m) & isfinite(x) & isfinite(y);
        if any(idx)
            scatter(ax, x(idx), y(idx), 14, P.morphColors(m,:), ...
                'filled', 'MarkerFaceAlpha',0.28, ...
                'DisplayName',morphLabels(m));
        end
    end
    legend(ax,'Location','best');
end

function addIdentityLineFromDataAudit(ax, x, y)
    values = [double(x(:)); double(y(:))];
    values = values(isfinite(values));
    if isempty(values); return; end
    lo = min(values); hi = max(values);
    if lo == hi
        lo = lo - 1; hi = hi + 1;
    end
    plot(ax,[lo hi],[lo hi],'k--','LineWidth',1.3,'HandleVisibility','off');
    xlim(ax,[lo hi]); ylim(ax,[lo hi]);
end

function addAuditLabels(ax, x, y, labels)
    for i = 1:numel(y)
        if isfinite(x(i)) && isfinite(y(i))
            text(ax,x(i)+0.015,y(i),string(labels(i)), ...
                'FontSize',7,'Interpreter','none');
        end
    end
end

