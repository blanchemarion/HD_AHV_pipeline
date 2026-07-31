function pipelineResults = fit_neuron_HD_AHV_four_models_behavior_tuned_cells_only_ROC(cfg)
%FIT_NEURON_HD_AHV_FOUR_MODELS_BEHAVIOR_TUNED_CELLS_ONLY_ROC
% Fit Surface, Molino and Pachon in one reproducible pipeline run.
% Pass the struct returned by pipeline_config; outputs are centralized.

% First select network-phase-tuned candidate neurons using the exact
% ORI_V15 STEP 17 information/shuffle logic, then fit the augmented HD/AHV
% model only to the selected neurons.
%
% Tuning selection is performed separately within each fish:
%   1) z-score every candidate-neuron trace over time;
%   2) at every frame, clip across neurons to the 2nd--98th percentiles and
%      subtract the across-neuron mean (ORI_V15 clip_center_framewise);
%   3) keep the positive part of each trace;
%   4) compute 36-bin, circularly smoothed Skaggs information about phi(t);
%   5) compare observed information with 500 circular shifts (10--90% of
%      the valid sequence), then apply BH-FDR across candidate cells;
%   6) fit the model below only when q_FDR < 0.05.
%
% IMPORTANT: This reproduces ORI_V15's original population-phase test.
% Because each candidate cell may contribute to the population phase, the
% selection remains partly circular unless leave-one-out phase is used.
%
% The fitted model is augmented with forward-bout and continuous-vigor controls.
%
% For every phase-tuned candidate neuron i:
%   y_i(t) = a_i*cos(phi(t)) + d_i*sin(phi(t)) ...
%          + b1_i*abs(AHV(t)) + b2_i*AHV(t) ...
%          + bF_i*F(t) + bV_i*V(t) + c0_i + epsilon_i(t)
%
% where:
%   F(t) = forward-bout impulses, preferably weighted by amp70, convolved
%          with the same causal calcium kernel used for AHV.
%   V(t) = continuous vigor, median-centered, convolved with a unit-area
%          causal calcium kernel, and sampled at calcium times.
%
% Derived after fitting:
%   b0_i        = hypot(a_i,d_i)
%   prefPhase_i = atan2(d_i,a_i)
%
% The script fits four models on the exact same samples. All four models use
% identical valid-frame masks and identical contiguous CV folds:
%   phase-only: cos(phi)+sin(phi)
%   original: cos(phi)+sin(phi)+|AHV|+AHV
%   behavior: |AHV|+AHV+F+V
%   augmented: cos(phi)+sin(phi)+|AHV|+AHV+F+V
%
% The phase-only model is:
%   y_i(t) = a_i*cos(phi(t)) + d_i*sin(phi(t)) + c0_i + epsilon_i(t)
%
% Its separately saved parameters are phaseOnlyA, phaseOnlyD, phaseOnlyB0,
% phaseOnlyPrefRad/Deg and phaseOnlyC0. phaseOnlyB0 is derived as
% hypot(phaseOnlyA,phaseOnlyD), and phaseOnlyPrefRad as
% atan2(phaseOnlyD,phaseOnlyA).
%
% The behavior-only model is:
%   y_i(t) = b1_i*abs(AHV(t)) + b2_i*AHV(t) ...
%          + bF_i*F(t) + bV_i*V(t) + c0_i + epsilon_i(t)
%
% Unlike earlier versions, this script retains the behavior-only
% coefficients, standard errors, p-values, standardized coefficients,
% in-sample R2/adjusted R2, AIC/BIC, residual scale, blocked-CV R2/SSE and
% aligned out-of-fold predictions. The original model also retains AIC/BIC
% and standardized coefficients. Thus, all four models are available as
% complete fitted models rather than the behavior-only model existing only
% inside the cross-validation loop.
%
% Primary control metric:
%   cvUniquePhase = 1 - SSE_CV(augmented)/SSE_CV(behavior)
% Negative out-of-sample values are retained; they are not clipped to zero.
%
% Additional phase-only comparisons:
%   cvUniqueAHVBeyondPhase =
%       1 - SSE_CV(original)/SSE_CV(phase-only)
%   cvUniqueAllBehaviorBeyondPhase =
%       1 - SSE_CV(augmented)/SSE_CV(phase-only)
%
% For downstream ROC/AUC analysis, the script also saves the aligned
% out-of-fold observed activity and predictions from all four models in
% fitTable.cvObserved, fitTable.cvPredPhaseOnly, fitTable.cvPredOriginal,
% fitTable.cvPredBehavior and fitTable.cvPredFull (the augmented-model
% compatibility name). Predictions at a frame
% always come from a model trained without that frame's contiguous test fold.
%
% Fish/session folders are discovered automatically for all three morphs.

if nargin < 1 || isempty(cfg)
    cfg = pipeline_config();
end
assert(isstruct(cfg) && isfield(cfg, 'Morphs') && isfield(cfg, 'ModelDataDir'), ...
    'Pass the struct returned by pipeline_config.');

%% ======================== USER PARAMETERS ========================
P = struct();

P.rootDataDir = '';

P.saveFigures = true;
P.figureVisible = cfg.FigureVisible;
P.figureFormat = 'png';                 % Primary format; SVG is also always saved

% Session discovery.
P.sessionDiscovery.includeRootAsSession = false;
P.sessionDiscovery.candidatePattern = '*candidate_neurons_clean.mat';
P.sessionDiscovery.excludeCandidateTokens = {'_classified','_direct_model'};
P.sessionDiscovery.behaviorPatterns = { ...
    '*swimResults_pass2*.mat', '*swimResults*.mat', ...
    '*BEHAVIOR*.mat', '*behavior*.mat'};
P.sessionDiscovery.preferBehaviorToken = 'PRE_FIRST_EVENT';
P.sessionDiscovery.skipSessionNames = {};
P.sessionDiscovery.skipNameContains = {};

% Common calcium/behavior support.
P.useTimeOverlapOnly = true;
P.behaviorTimeMode = 'auto';  % 'auto', 'local', 'startFrame', 'behaviorTimeSec'
P.model.zscoreActivityWithinFitWindow = true;
P.minSamplesForNeuronFit = 100;

% ORI_V15 STEP 17 phase-tuning selection. Do not change these values if
% exact comparability with that script is required.
P.tuning.nPhaseBins = 36;
P.tuning.smoothBins = 3;
P.tuning.nShuffles = 500;
P.tuning.minShiftFraction = 0.10;
P.tuning.maxShiftFraction = 0.90;
P.tuning.fdrAlpha = 0.05;
P.tuning.percentileLow = 2;
P.tuning.percentileHigh = 98;
P.tuning.minValidPhaseFrames = 20;
P.tuning.randomSeed = 1;      % makes the shuffle selection reproducible
P.tuning.usePhaseOkIfAvailable = true;
P.tuning.requirePhaseOk = false; % ORI_V15 falls back to all finite phases
P.tuning.saveShuffleMatrix = true;

% Causal calcium kernel. AHV and F use event impulses; V uses a unit-area
% continuous convolution. Keep these values identical to the direct model.
P.kernel.tauSeconds = 3.0;
P.kernel.cutoffTau = 6;
P.ahvBoutTime = 'start';       % 'start', 'center', or 'end'

% Forward-bout definition and weighting.
% First choice: |LI| <= LI_threshold. If unavailable, use dtheta threshold.
P.forward.useLIThresholdWhenAvailable = true;
P.forward.fallbackAbsDthetaDeg = 5;
P.forward.useAmp70Weights = true;
P.forward.allowUnitWeightFallback = true;
P.forward.replaceMissingAmp70WithMedian = true;

% Continuous vigor control.
P.vigor.requireField = true;            % skip fish if vigor is unavailable
P.vigor.centerBeforeFiltering = true;   % removes baseline/intercept redundancy
P.vigor.maximumMissingFraction = 0.20;

% Same low-turn inclusion rule used in the behavioral summary.
P.lowTurnQC.enabled = true;
P.lowTurnQC.minRotationsEachDirection = 0; % matches the attached direct-model script

% Blocked time-series cross-validation. All nested models use identical folds.
P.cv.enabled = true;
P.cv.nBlockedFolds = 5;
P.cv.minTestSamplesPerFold = 20;
P.cv.saveOOFTraces = true; % required by the downstream ROC/AUC analysis

% CW/CCW/Symmetric labels from the augmented model. These remain secondary:
% the main purpose of this script is continuous-parameter/variance analysis.
P.class.alpha = 0.05;
P.class.bonferroniMode = 'fixed';       % 'fixed' or 'per_fish'
P.class.fixedNTests = 309;
P.ahvPositiveIsCCW = true;

% Optional non-destructive classified candidate copy.
P.classificationExport.enabled = true;
P.classificationExport.overwrite = true;
P.classificationExport.suffix = '_classified';
P.classificationExport.outputDir = '';

% QC plots.
P.qc.maxTimePoints = 12000;
P.qc.vifReference = 5;

%% ======================== RUN ALL MORPHS ========================
baseP = P;
pipelineResults = struct();
pipelineResults.config = cfg;
pipelineResults.morphs = repmat(struct( ...
    'morph', '', 'resultFile', '', 'sessionDataDir', '', ...
    'figureDir', '', 'nSessions', 0, 'nFittedNeurons', 0), ...
    numel(cfg.Morphs), 1);

for m = 1:numel(cfg.Morphs)
    morph = cfg.Morphs(m);
    P = baseP;
    P.rootDataDir = morph.dataDir;
    P.morphName = morph.name;
    P.morphDisplayName = morph.display;
    P.figureNamePrefix = [morph.name '_'];
    P.classificationExport.outputDir = fullfile( ...
        cfg.ClassifiedCandidateDir, morph.name);
    sessionDataDir = fullfile(cfg.SessionModelDataDir, morph.name);
    figureDir = fullfile(cfg.ModelFigureDir, morph.name);
    ensureDirectory(sessionDataDir);
    ensureDirectory(figureDir);
    ensureDirectory(P.classificationExport.outputDir);

    disp(' ');
    disp('############################################################');
    fprintf('Fitting morph %s (%d/%d)', morph.display, m, numel(cfg.Morphs));
    fprintf('%s', newline);
    fprintf('Raw input: %s', P.rootDataDir); fprintf('%s', newline);
    fprintf('Processed data: %s', cfg.ModelDataDir); fprintf('%s', newline);
    fprintf('Figures: %s', figureDir); fprintf('%s', newline);
    pipelineResults.morphs(m) = runOneMorph( ...
        P, cfg.ModelDataDir, sessionDataDir, figureDir);
end

save(fullfile(cfg.ModelDataDir, 'four_model_fit_pipeline_summary.mat'), ...
    'pipelineResults', '-v7.3');
disp(' ');
fprintf('All morphs complete. Processed data: %s', cfg.ModelDataDir);
fprintf('%s', newline);
end

function morphOutput = runOneMorph(P, modelDataDir, sessionDataDir, figureDir)
%% ======================== DISCOVER / RUN ========================
[sessions, P.rootDataDir] = discoverSessions(P.rootDataDir, P);

Results = cell(numel(sessions),1);
allTables = cell(numel(sessions),1);
fishRows = cell(numel(sessions),1);
allTuningTables = cell(numel(sessions),1);
tuningFishRows = cell(numel(sessions),1);

for s = 1:numel(sessions)
    fprintf('\n================ %s (%d/%d) ================\n', ...
        sessions(s).name, s, numel(sessions));
    try
        Results{s} = analyzeOneSession(sessions(s), P, sessionDataDir, figureDir);
        allTables{s} = Results{s}.fitTable;
        fishRows{s} = Results{s}.fishSummary;
        allTuningTables{s} = Results{s}.tuningTable;
        tuningFishRows{s} = Results{s}.tuningSummary;
    catch ME
        warning('Session %s failed: %s', sessions(s).name, ME.message);
        Results{s} = struct('name',sessions(s).name,'error',ME);
        allTables{s} = table();
        fishRows{s} = table();
        allTuningTables{s} = table();
        tuningFishRows{s} = table();
    end
end

useT = ~cellfun(@isempty,allTables) & cellfun(@(x) istable(x) && height(x)>0,allTables);
if any(useT); AllNeurons = vertcat(allTables{useT}); else; AllNeurons = table(); end
useF = ~cellfun(@isempty,fishRows) & cellfun(@(x) istable(x) && height(x)>0,fishRows);
if any(useF); FishSummary = vertcat(fishRows{useF}); else; FishSummary = table(); end
useQ = ~cellfun(@isempty,allTuningTables) & ...
    cellfun(@(x) istable(x) && height(x)>0,allTuningTables);
if any(useQ)
    AllCandidateTuning = vertcat(allTuningTables{useQ});
else
    AllCandidateTuning = table();
end
useQS = ~cellfun(@isempty,tuningFishRows) & ...
    cellfun(@(x) istable(x) && height(x)>0,tuningFishRows);
if any(useQS)
    TuningFishSummary = vertcat(tuningFishRows{useQS});
else
    TuningFishSummary = table();
end

nameStem = ['HD_AHV_behavior_augmented_phase_tuned_only_' ...
    'with_phase_only'];
resultFile = fullfile(modelDataDir, sprintf('%s_results_%s.mat', ...
    nameStem, P.morphName));
save(resultFile, ...
    'Results','AllNeurons','FishSummary','AllCandidateTuning', ...
    'TuningFishSummary','P','sessions','-v7.3');
if ~isempty(AllNeurons)
    % Variable-length OOF traces are retained in the MAT file but omitted
    % from the flat CSV, which is intended for scalar neuron summaries.
    AllNeuronsCSV = AllNeurons;
    traceVariables = intersect({'cvObserved','cvPredPhaseOnly', ...
        'cvPredOriginal','cvPredBehavior','cvPredFull','cvPredAugmented'}, ...
        AllNeuronsCSV.Properties.VariableNames,'stable');
    if ~isempty(traceVariables)
        AllNeuronsCSV = removevars(AllNeuronsCSV,traceVariables);
    end
    writetable(AllNeuronsCSV,fullfile(modelDataDir, sprintf('%s_neurons_%s.csv', nameStem, P.morphName)));
end
if ~isempty(FishSummary)
    writetable(FishSummary,fullfile(modelDataDir, sprintf('%s_fish_summary_%s.csv', nameStem, P.morphName)));
end
if ~isempty(AllCandidateTuning)
    writetable(AllCandidateTuning,fullfile(modelDataDir, sprintf('HD_AHV_phase_tuning_all_candidates_%s.csv', P.morphName)));
end
if ~isempty(TuningFishSummary)
    writetable(TuningFishSummary,fullfile(modelDataDir, sprintf('HD_AHV_phase_tuning_fish_summary_%s.csv', P.morphName)));
end

if ~isempty(AllNeurons)
    plotAllFishModelComparison(AllNeurons,P,figureDir);
    plotAllFishCoefficientSummary(AllNeurons,P,figureDir);
end
if ~isempty(FishSummary)
    plotFishLevelSummary(FishSummary,P,figureDir);
end

morphOutput = struct('morph', P.morphDisplayName, ...
    'resultFile', resultFile, 'sessionDataDir', sessionDataDir, ...
    'figureDir', figureDir, 'nSessions', numel(sessions), ...
    'nFittedNeurons', height(AllNeurons));
fprintf('Morph %s complete: %s', P.morphDisplayName, resultFile);
fprintf('%s', newline);
end

%% ======================== LOCAL FUNCTIONS ========================

function ensureDirectory(pathName)
if exist(pathName, 'dir') ~= 7
    mkdir(pathName);
end
end

function [sessions,rootDir] = discoverSessions(rootDir,P)
if isempty(rootDir) || exist(rootDir,'dir') ~= 7
    rootDir = uigetdir(pwd,'Select morph root folder');
    assert(~isequal(rootDir,0),'No root folder selected.');
end

d = dir(rootDir);
d = d([d.isdir] & ~ismember({d.name},{'.','..'}));
if P.sessionDiscovery.includeRootAsSession
    rootEntry = struct('name',getFolderName(rootDir),'folder',fileparts(rootDir), ...
        'date','','bytes',0,'isdir',true,'datenum',0);
    d = [rootEntry; d(:)]; %#ok<AGROW>
end

sessions = struct('name',{},'dataDir',{},'candidateFile',{},'behaviorFile',{});
skipped = strings(0,1);
for k = 1:numel(d)
    name = d(k).name;
    if any(strcmpi(name,P.sessionDiscovery.skipSessionNames)) || ...
            any(cellfun(@(x) contains(lower(name),lower(x)), ...
            P.sessionDiscovery.skipNameContains))
        skipped(end+1,1) = string(name) + " (manual skip)"; %#ok<AGROW>
        continue;
    end
    if P.sessionDiscovery.includeRootAsSession && strcmp(name,getFolderName(rootDir))
        folder = rootDir;
    else
        folder = fullfile(rootDir,name);
    end
    try
        candidate = pickCandidateFile(folder,P);
        behavior = pickBehaviorFile(folder,P);
        sessions(end+1) = struct('name',name,'dataDir',folder, ...
            'candidateFile',candidate,'behaviorFile',behavior); %#ok<AGROW>
    catch ME
        skipped(end+1,1) = string(name) + ": " + string(ME.message); %#ok<AGROW>
    end
end

assert(~isempty(sessions),'No valid session folders found under %s.',rootDir);
fprintf('Found %d valid sessions under %s\n',numel(sessions),rootDir);
for k = 1:numel(sessions)
    fprintf('  %s | candidate: %s | behavior: %s\n',sessions(k).name, ...
        sessions(k).candidateFile,sessions(k).behaviorFile);
end
if ~isempty(skipped)
    fprintf('Skipped %d folder(s):\n',numel(skipped));
    for k=1:numel(skipped); fprintf('  %s\n',char(skipped(k))); end
end
end

function name = pickCandidateFile(folder,P)
f = dir(fullfile(folder,P.sessionDiscovery.candidatePattern));
f = f(~[f.isdir]);
for j = numel(f):-1:1
    low = lower(f(j).name);
    if any(cellfun(@(x) contains(low,lower(x)), ...
            P.sessionDiscovery.excludeCandidateTokens))
        f(j) = [];
    end
end
assert(~isempty(f),'No clean candidate-neuron file.');
if numel(f)>1
    [~,ord] = sort([f.datenum],'descend'); f = f(ord);
    warning('Multiple candidate files in %s; using %s.',folder,f(1).name);
end
name = f(1).name;
end

function name = pickBehaviorFile(folder,P)
f = struct([]);
for j = 1:numel(P.sessionDiscovery.behaviorPatterns)
    q = dir(fullfile(folder,P.sessionDiscovery.behaviorPatterns{j}));
    q = q(~[q.isdir]);
    if ~isempty(q); f = [f; q(:)]; end %#ok<AGROW>
end
assert(~isempty(f),'No behavior file matching configured patterns.');
fullPaths=fullfile({f.folder},{f.name});
lowPaths=cellfun(@lower,fullPaths,'UniformOutput',false);
[~,ia] = unique(lowPaths,'stable'); f = f(ia);

lowNames=cellfun(@lower,{f.name},'UniformOutput',false);
isPass2 = contains(lowNames,'swimresults_pass2');
isPreferred = contains(lowNames,lower(P.sessionDiscovery.preferBehaviorToken));
score = 10*double(isPass2) + 3*double(isPreferred);
[~,ord] = sortrows([score(:),[f.datenum]'],[-1 -2]); f = f(ord);

valid = false(numel(f),1);
for j = 1:numel(f)
    try
        w = whos('-file',fullfile(f(j).folder,f(j).name));
        valid(j) = any(strcmp({w.name},'pass2FileResult'));
    catch
        valid(j) = false;
    end
end
f = f(valid);
assert(~isempty(f),'No file containing pass2FileResult.');
if numel(f)>1
    warning('Multiple behavior files in %s; using %s.',folder,f(1).name);
end
name = f(1).name;
end

function Result = analyzeOneSession(sess,P,dataOutDir,figureOutDir)
candidatePath = fullfile(sess.dataDir,sess.candidateFile);
behaviorPath = fullfile(sess.dataDir,sess.behaviorFile);

C = load(candidatePath);
[Y,tCa,phase,phaseOk,sourcePrefRad,candidateIDs,loadDiag] = loadCandidate(C,P);

% Selection is run on every candidate cell and on the full candidate-file
% time support, before behavior alignment or model fitting.
Tuning = selectPhaseTunedCellsORI(Y,phase,phaseOk,candidateIDs,sess.name,P);
Tuning.phaseOkSource = loadDiag.phaseOkSource;
tuningTable = makeTuningTable(sess.name,candidateIDs,Tuning);
tuningSummary = makeTuningSummary(sess.name,Tuning);
plotPhaseTuningSelection(sess.name,Tuning,P,figureOutDir);

allCandidateIDs = candidateIDs;
tunedCandidateOrder = find(Tuning.isPhaseTuned);
fprintf('ORI_V15 phase-tuned selection: %d/%d candidates (BH-FDR q < %.3g).\n', ...
    numel(tunedCandidateOrder),numel(candidateIDs),P.tuning.fdrAlpha);

if isempty(tunedCandidateOrder)
    warning('No phase-tuned cells in %s; no model was fitted.',sess.name);
    Result = struct('name',sess.name,'excludedByLowTurnQC',false, ...
        'excludedBecauseNoTunedCells',true, ...
        'paths',struct('candidatePath',candidatePath,'behaviorPath',behaviorPath), ...
        'candidateDiagnostics',loadDiag,'tuning',Tuning, ...
        'tuningTable',tuningTable,'tuningSummary',tuningSummary, ...
        'fitTable',table(),'fishSummary',table());
    save(fullfile(dataOutDir,[sess.name '_phase_tuned_only_with_phase_only_model.mat']), ...
        'Result','P','-v7.3');
    return;
end

Y = Y(:,tunedCandidateOrder);
sourcePrefRad = sourcePrefRad(tunedCandidateOrder);
candidateIDs = candidateIDs(tunedCandidateOrder);

S = load(behaviorPath,'pass2FileResult');
assert(isfield(S,'pass2FileResult'),'Behavior file lacks pass2FileResult.');
beh = S.pass2FileResult;

[B,behaviorDiag] = buildBehaviorRegressors(beh,tCa,P);
if P.lowTurnQC.enabled && ...
        behaviorDiag.minDirectionalRotations < P.lowTurnQC.minRotationsEachDirection
    fprintf('EXCLUDED: positive rotations %.2f, negative rotations %.2f.\n', ...
        behaviorDiag.totalPositiveRotations,behaviorDiag.totalNegativeRotations);
    Result = struct('name',sess.name,'excludedByLowTurnQC',true, ...
        'excludedBecauseNoTunedCells',false, ...
        'paths',struct('candidatePath',candidatePath,'behaviorPath',behaviorPath), ...
        'behaviorDiagnostics',behaviorDiag,'candidateDiagnostics',loadDiag, ...
        'tuning',Tuning,'tuningTable',tuningTable, ...
        'tuningSummary',tuningSummary,'fitTable',table(),'fishSummary',table());
    save(fullfile(dataOutDir,[sess.name '_phase_tuned_only_with_phase_only_model.mat']), ...
        'Result','P','-v7.3');
    return;
end

common = isfinite(tCa) & isfinite(phase) & isfinite(B.ahv) & ...
         isfinite(B.forward) & isfinite(B.vigor);
if P.useTimeOverlapOnly
    Y = Y(common,:); tUsed = tCa(common); phase = phase(common);
    ahv = B.ahv(common); forward = B.forward(common); vigor = B.vigor(common);
else
    tUsed = tCa; ahv = B.ahv; forward = B.forward; vigor = B.vigor;
end
assert(sum(all(isfinite([phase,ahv,forward,vigor]),2)) >= P.minSamplesForNeuronFit, ...
    'Too few common finite samples after behavior alignment.');

if P.model.zscoreActivityWithinFitWindow; Y = zscoreColumnsFinite(Y); end

predictors = [cos(phase),sin(phase),abs(ahv),ahv,forward,vigor];
predictorNames = {'cosPhi','sinPhi','absAHV','AHV','forward','vigor'};
[predictorCorr,predictorVIF,conditionNumberStandardized] = ...
    predictorDiagnostics(predictors);

fprintf('Fitting phase-only, original, augmented and behavior-only models: %d neurons, %d frames.\n', ...
    size(Y,2),size(Y,1));
Fit = fitAllNeurons(Y,phase,ahv,forward,vigor,sourcePrefRad,P);
[classLabel,ahsLabel,classInfo] = classifyAugmented(Fit,P);

N = size(Y,2);
session = repmat({sess.name},N,1);
fitTable = makeFitTable(session,tunedCandidateOrder,candidateIDs,Fit, ...
    classLabel,ahsLabel,classInfo,predictorVIF,conditionNumberStandardized);
fitTable.phaseInfoBits = Tuning.infoBits(tunedCandidateOrder);
fitTable.phaseInfoShuffleP = Tuning.pShuffle(tunedCandidateOrder);
fitTable.phaseInfoFDRQ = Tuning.qFDR(tunedCandidateOrder);
fitTable.phaseTuningPrefRad = Tuning.prefRad(tunedCandidateOrder);
fitTable.phaseTuningPrefDeg = rad2deg(Tuning.prefRad(tunedCandidateOrder));
fitTable.phaseTuningVectorStrength = Tuning.vectorStrength(tunedCandidateOrder);
fitTable.selectedByPhaseTuning = true(N,1);

fishSummary = makeFishSummary(sess.name,fitTable,behaviorDiag, ...
    predictorCorr,predictorVIF,conditionNumberStandardized);
fishSummary.nCandidatesPhaseTested = numel(allCandidateIDs);
fishSummary.nPhaseTunedSelected = numel(tunedCandidateOrder);
fishSummary.phaseTunedFraction = numel(tunedCandidateOrder)/numel(allCandidateIDs);
fishSummary.phaseTuningFDRAlpha = P.tuning.fdrAlpha;

Result = struct();
Result.name = sess.name;
Result.excludedByLowTurnQC = false;
Result.excludedBecauseNoTunedCells = false;
Result.paths = struct('candidatePath',candidatePath,'behaviorPath',behaviorPath);
Result.tCa = tUsed;
Result.phaseCa = phase;
Result.ahvCa = ahv;
Result.forwardCa = forward;
Result.vigorCa = vigor;
Result.predictorNames = predictorNames;
Result.predictorCorrelation = predictorCorr;
Result.predictorVIF = predictorVIF;
Result.conditionNumberStandardized = conditionNumberStandardized;
Result.behaviorDiagnostics = behaviorDiag;
Result.candidateDiagnostics = loadDiag;
Result.tuning = Tuning;
Result.tuningTable = tuningTable;
Result.tuningSummary = tuningSummary;
Result.tunedCandidateOrder = tunedCandidateOrder;
Result.tunedCandidateIDs = candidateIDs;
Result.fit = Fit;
Result.classLabel = classLabel;
Result.classInfo = classInfo;
Result.fitTable = fitTable;
Result.fishSummary = fishSummary;

perFishFile = fullfile(dataOutDir,[sess.name '_phase_tuned_only_with_phase_only_model.mat']);
save(perFishFile,'Result','P','-v7.3');
plotRegressorQC(Result,P,figureOutDir);
plotPredictorCorrelation(Result,P,figureOutDir);
plotPerFishModelComparison(Result,P,figureOutDir);

if P.classificationExport.enabled
    exportClassifiedCandidate(candidatePath,behaviorPath,allCandidateIDs,tunedCandidateOrder, ...
        classLabel,fitTable,classInfo,Tuning,P);
end
end

function [Y,tCa,phase,phaseOk,sourcePrefRad,candidateIDs,D] = loadCandidate(C,P)
required = {'calcium_traces','time_s','network_phase_rad','candidate_cell_ids'};
for k = 1:numel(required)
    assert(isfield(C,required{k}),'Candidate file missing %s.',required{k});
end
candidateIDs = double(C.candidate_cell_ids(:));
Y0 = double(C.calcium_traces);
if size(Y0,2)==numel(candidateIDs); Y=Y0;
elseif size(Y0,1)==numel(candidateIDs); Y=Y0';
else; error('calcium_traces dimensions do not match candidate_cell_ids.'); end
tCa = double(C.time_s(:));
phase = wrapToPiLocal(double(C.network_phase_rad(:)));
assert(numel(tCa)==size(Y,1) && numel(phase)==size(Y,1), ...
    'Candidate time/phase lengths do not match calcium traces.');
assert(all(diff(tCa(isfinite(tCa)))>0),'time_s must increase strictly.');
[phaseOk,phaseOkSource] = getPhaseOk(C,numel(phase),P);
[sourcePrefRad,prefSource] = getSourcePreferredPhase(C,size(Y,2));
D = struct('activitySource','calcium_traces','phaseSource','network_phase_rad', ...
    'phaseOkSource',phaseOkSource, ...
    'referencePreferredPhaseSource',prefSource,'nFrames',size(Y,1), ...
    'nNeurons',size(Y,2));
end

function [phaseOk,source] = getPhaseOk(C,nFrames,P)
phaseOk = true(nFrames,1);
source = 'all finite phase frames (phase_ok unavailable)';
if ~P.tuning.usePhaseOkIfAvailable
    source = 'all finite phase frames (phase_ok use disabled)';
    return;
end

candidate = [];
topLevelNames = {'phase_ok','network_phase_ok','phase_valid_mask', ...
    'network_phase_valid','manual_phase_ok'};
for k = 1:numel(topLevelNames)
    if isfield(C,topLevelNames{k}) && numel(C.(topLevelNames{k}))==nFrames
        candidate = C.(topLevelNames{k});
        source = topLevelNames{k};
        break;
    end
end

% Also support candidate files that retain the original ORI_V15 hdMetrics.
if isempty(candidate) && isfield(C,'hdMetrics') && isstruct(C.hdMetrics) && ...
        isfield(C.hdMetrics,'r1pi') && isfield(C.hdMetrics.r1pi,'phase')
    tags = {'r1pi','all'};
    for k = 1:numel(tags)
        tag = tags{k};
        if isfield(C.hdMetrics.r1pi.phase,tag) && ...
                isfield(C.hdMetrics.r1pi.phase.(tag),'phase_ok') && ...
                numel(C.hdMetrics.r1pi.phase.(tag).phase_ok)==nFrames
            candidate = C.hdMetrics.r1pi.phase.(tag).phase_ok;
            source = ['hdMetrics.r1pi.phase.' tag '.phase_ok'];
            break;
        end
    end
end

if isempty(candidate)
    if P.tuning.requirePhaseOk
        error('No phase_ok vector of length %d was found in the candidate file.',nFrames);
    end
    warning(['phase_ok was not found; matching ORI_V15 fallback by using all ' ...
        'finite network-phase frames.']);
    return;
end

phaseOk = logical(candidate(:));
phaseOk(~isfinite(double(candidate(:)))) = false;
end

function Tuning = selectPhaseTunedCellsORI(Y,phase,phaseOk,candidateIDs,sessionName,P)
% Reimplementation of ORI_V15 STEP 17 on candidate-file variables.
assert(size(Y,1)==numel(phase) && numel(phaseOk)==numel(phase), ...
    'Activity, phase and phase_ok must have matching frame counts.');
assert(size(Y,2)==numel(candidateIDs), ...
    'Activity columns must match candidate_cell_ids.');

oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(P.tuning.randomSeed,'twister');

% Exact ORI_V15 preprocessing: timewise z-score, framewise clipping and
% centering across candidate cells, followed by positive rectification.
dataZ = zscore(double(Y),0,1);
dataZ(~isfinite(dataZ)) = 0;
pLow = prctile(dataZ,P.tuning.percentileLow,2);
pHigh = prctile(dataZ,P.tuning.percentileHigh,2);
Fz = bsxfun(@min,dataZ,pHigh);
Fz = bsxfun(@max,Fz,pLow);
Fz = bsxfun(@minus,Fz,mean(Fz,2));

validFrame = isfinite(phase(:)) & logical(phaseOk(:));
phiUse = phase(validFrame);
FzUse = Fz(validFrame,:);
nValid = numel(phiUse);
assert(nValid>=P.tuning.minValidPhaseFrames, ...
    'Only %d valid phase frames; at least %d are required.', ...
    nValid,P.tuning.minValidPhaseFrames);

nb = P.tuning.nPhaseBins;
edges = linspace(-pi,pi,nb+1);
centers = ((edges(1:end-1)+edges(2:end))/2).';
phaseBin = discretize(phiUse,edges);
validBin = ~isnan(phaseBin);
occupancyCount = accumarray(phaseBin(validBin),1,[nb,1],@sum,0);
occupancy = occupancyCount/max(sum(occupancyCount),1);

nCells = size(FzUse,2);
infoBits = nan(nCells,1);
prefRad = nan(nCells,1);
vectorStrength = nan(nCells,1);
rateBins = nan(nb,nCells);

for k = 1:nCells
    activityPositive = max(FzUse(:,k),0);
    rb = accumarray(phaseBin(validBin),activityPositive(validBin), ...
        [nb,1],@mean,0);
    rb = circularSmoothBins(rb,P.tuning.smoothBins);
    [infoBits(k),prefRad(k),vectorStrength(k)] = ...
        phaseInfoFromBins(rb,occupancy,centers);
    rateBins(:,k) = rb;
end

nShuffles = P.tuning.nShuffles;
shuffleInfo = nan(nCells,nShuffles);
nSequence = nValid;
minShift = max(1,floor(P.tuning.minShiftFraction*nSequence));
maxShift = max(1,floor(P.tuning.maxShiftFraction*nSequence));
assert(maxShift>minShift, ...
    'Invalid circular-shift range [%d %d] for %d valid frames.', ...
    minShift,maxShift,nSequence);

fprintf('Running %d ORI_V15 circular shifts for each of %d candidate cells...\n', ...
    nShuffles,nCells);
for k = 1:nCells
    activityPositive = max(FzUse(:,k),0);
    for s = 1:nShuffles
        shiftFrames = randi([minShift maxShift]);
        shifted = circshift(activityPositive,shiftFrames);
        rb = accumarray(phaseBin(validBin),shifted(validBin), ...
            [nb,1],@mean,0);
        rb = circularSmoothBins(rb,P.tuning.smoothBins);
        shuffleInfo(k,s) = phaseInfoFromBins(rb,occupancy,centers);
    end
end

% This intentionally matches ORI_V15: no +1 pseudocount is added.
pShuffle = mean(shuffleInfo>=infoBits,2,'omitnan');
try
    qFDR = mafdr(pShuffle,'BHFDR',true);
catch
    qFDR = bhFdrLocal(pShuffle);
end
isPhaseTuned = isfinite(qFDR) & qFDR<P.tuning.fdrAlpha;
null95 = prctile(shuffleInfo,95,2);

Tuning = struct();
Tuning.method = 'ORI_V15_STEP17_information_circular_shift_BH_FDR';
Tuning.session = sessionName;
Tuning.candidateIDs = candidateIDs(:);
Tuning.infoBits = infoBits;
Tuning.prefRad = prefRad;
Tuning.vectorStrength = vectorStrength;
Tuning.rateBins = rateBins;
Tuning.occupancy = occupancy;
Tuning.pShuffle = pShuffle;
Tuning.qFDR = qFDR;
Tuning.null95 = null95;
Tuning.isPhaseTuned = isPhaseTuned;
Tuning.validFrameMask = validFrame;
Tuning.nValidPhaseFrames = nValid;
Tuning.nCandidates = nCells;
Tuning.nSelected = nnz(isPhaseTuned);
Tuning.selectedCandidateOrder = find(isPhaseTuned);
Tuning.selectedCandidateIDs = candidateIDs(isPhaseTuned);
Tuning.params = P.tuning;
Tuning.params.phaseEdges = edges;
Tuning.params.phaseCenters = centers;
Tuning.params.minShiftFrames = minShift;
Tuning.params.maxShiftFrames = maxShift;
Tuning.preprocessing = ['zscore over time; clip each frame across candidate ' ...
    'cells to percentileLow/percentileHigh; subtract framewise mean; ' ...
    'positive rectification'];
Tuning.circularityCaveat = ['Cells are tested against a population phase to ' ...
    'which they may contribute; circular shifts preserve temporal dynamics ' ...
    'but do not remove this self-inclusion.'];
if P.tuning.saveShuffleMatrix
    Tuning.shuffleInfo = shuffleInfo;
else
    Tuning.shuffleInfo = [];
end
end

function T = makeTuningTable(sessionName,candidateIDs,Q)
n = numel(candidateIDs);
session = repmat({sessionName},n,1);
candidateOrder = (1:n)';
nValidPhaseFrames = repmat(Q.nValidPhaseFrames,n,1);
T = table(session,candidateOrder,candidateIDs(:),Q.infoBits,Q.prefRad, ...
    rad2deg(Q.prefRad),Q.vectorStrength,Q.null95,Q.pShuffle,Q.qFDR, ...
    Q.isPhaseTuned,nValidPhaseFrames, ...
    'VariableNames',{'session','candidateOrder','neuronID','phaseInfoBits', ...
    'phaseTuningPrefRad','phaseTuningPrefDeg','phaseTuningVectorStrength', ...
    'phaseInfoShuffle95','phaseInfoShuffleP','phaseInfoFDRQ', ...
    'isPhaseTuned','nValidPhaseFrames'});
end

function T = makeTuningSummary(sessionName,Q)
T = table(string(sessionName),Q.nCandidates,Q.nSelected, ...
    Q.nSelected/max(Q.nCandidates,1),Q.nValidPhaseFrames, ...
    median(Q.infoBits,'omitnan'),median(Q.pShuffle,'omitnan'), ...
    median(Q.qFDR,'omitnan'),string(Q.phaseOkSource), ...
    Q.params.nShuffles,Q.params.fdrAlpha,Q.params.randomSeed, ...
    'VariableNames',{'session','nCandidatesPhaseTested','nPhaseTunedSelected', ...
    'phaseTunedFraction','nValidPhaseFrames','medianPhaseInfoBits', ...
    'medianShuffleP','medianFDRQ','phaseOkSource','nShuffles', ...
    'fdrAlpha','randomSeed'});
end

function y = circularSmoothBins(x,width)
x = x(:);
if width<=1; y=x; return; end
kernel = ones(round(width),1)/round(width);
n = numel(x);
triple = conv([x;x;x],kernel,'same');
y = triple((n+1):(2*n));
end

function [I,pref,vecR] = phaseInfoFromBins(rate,occupancy,centers)
rate = rate(:); occupancy = occupancy(:); centers = centers(:);
use = occupancy>0 & isfinite(rate) & isfinite(occupancy);
if ~any(use)
    I=0; pref=NaN; vecR=0; return;
end
rbar = sum(occupancy(use).*rate(use));
if ~isfinite(rbar) || rbar<=0
    I=0; pref=NaN; vecR=0; return;
end
ratio = rate(use)/rbar;
I = sum(occupancy(use).*ratio.*log2(ratio+eps));
v = sum(rate(use).*exp(1i*centers(use)));
if sum(rate(use))>0; v=v/sum(rate(use)); else; v=0; end
pref = angle(v);
vecR = abs(v);
end

function q = bhFdrLocal(p)
p = p(:); q = nan(size(p)); valid = isfinite(p); pv = p(valid);
[sortedP,order] = sort(pv,'ascend'); m = numel(sortedP);
if m==0; return; end
sortedQ = sortedP.*m./(1:m)';
for k=m-1:-1:1; sortedQ(k)=min(sortedQ(k),sortedQ(k+1)); end
sortedQ = min(sortedQ,1);
unsortedQ = nan(size(pv)); unsortedQ(order)=sortedQ; q(valid)=unsortedQ;
end

function [pref,source] = getSourcePreferredPhase(C,N)
pref = NaN(N,1); source = 'unavailable';
if isfield(C,'candidate_features') && istable(C.candidate_features)
    vn = C.candidate_features.Properties.VariableNames;
    radNames = {'prefRad','preferred_phase_rad','preferredPhaseRad','pref_phi_rad'};
    degNames = {'prefDeg','preferred_phase_deg','preferredPhaseDeg','pref_phi_deg'};
    for k=1:numel(radNames)
        q=find(strcmpi(vn,radNames{k}),1);
        if ~isempty(q) && height(C.candidate_features)==N
            pref=wrapToPiLocal(double(C.candidate_features{:,q})); source=vn{q}; return;
        end
    end
    for k=1:numel(degNames)
        q=find(strcmpi(vn,degNames{k}),1);
        if ~isempty(q) && height(C.candidate_features)==N
            pref=wrapToPiLocal(deg2rad(double(C.candidate_features{:,q}))); source=vn{q}; return;
        end
    end
end
if isfield(C,'prefRad') && numel(C.prefRad)==N
    pref=wrapToPiLocal(double(C.prefRad(:))); source='prefRad';
elseif isfield(C,'prefDeg') && numel(C.prefDeg)==N
    pref=wrapToPiLocal(deg2rad(double(C.prefDeg(:)))); source='prefDeg';
end
end

function [B,D] = buildBehaviorRegressors(beh,tCa,P)
% Construct all motor predictors at behavior resolution before sampling at
% calcium times. Event and continuous signals need different discrete
% convolution scaling; see comments below.
assert(isfield(beh,'fps') && isfinite(double(beh.fps)) && double(beh.fps)>0, ...
    'Behavior fps is missing or invalid.');
fps = double(beh.fps); fps = fps(1);
nBeh = inferBehaviorLength(beh);
assert(nBeh>=2,'Could not infer a valid behavior length.');

[tBeh,timeSource] = chooseBehaviorTime(beh,nBeh,fps,tCa,P.behaviorTimeMode);

if isfield(beh,'dtheta') && ~isempty(beh.dtheta)
    dtheta = cleanAngleUnitsToRad(double(beh.dtheta(:)));
    dthetaSource = 'dtheta';
elseif isfield(beh,'LI') && ~isempty(beh.LI)
    dtheta = cleanAngleUnitsToRad(double(beh.LI(:)));
    dthetaSource = 'LI fallback';
    warning('Using LI as the dtheta fallback.');
else
    error('Behavior file contains neither dtheta nor LI.');
end
assert(isfield(beh,'swimFrames') && ~isempty(beh.swimFrames), ...
    'swimFrames is required to construct event regressors.');
[boutStart,boutEnd] = parseSwimFrames(double(beh.swimFrames),numel(dtheta));
nB = min([numel(dtheta),numel(boutStart),numel(boutEnd)]);
dtheta=dtheta(1:nB); boutStart=boutStart(1:nB); boutEnd=boutEnd(1:nB);

switch lower(P.ahvBoutTime)
    case 'start';  boutFrame=boutStart;
    case 'center'; boutFrame=round((boutStart+boutEnd)./2);
    case 'end';    boutFrame=boutEnd;
    otherwise; error('P.ahvBoutTime must be start, center, or end.');
end
boutFrame = round(double(boutFrame(:)));
validBout = isfinite(boutFrame) & boutFrame>=1 & boutFrame<=nBeh & isfinite(dtheta);

% AHV impulses carry integrated angle (rad). Convolution with k(t), whose
% units are 1/s, therefore produces rad/s; no 1/fps factor is applied.
angleImpulse = accumarray(boutFrame(validBout),dtheta(validBout),[nBeh 1],@sum,0);
kernelT = (0:round(P.kernel.cutoffTau*P.kernel.tauSeconds*fps))'./fps;
kEvent = (1/P.kernel.tauSeconds).*exp(-kernelT/P.kernel.tauSeconds);
ahvBeh = causalConvolutionFFT(angleImpulse,kEvent);

% Determine forward bouts. This reproduces the behavioral summary's first
% choice (LI_threshold), with an explicit angle fallback.
[isForward,forwardDefinition] = identifyForwardBouts(beh,dtheta,nB,P);
isForward = isForward(:) & validBout;
[forwardWeight,forwardWeightSource,nMissingAmp] = ...
    getForwardWeights(beh,isForward,nB,P);
forwardImpulse = accumarray(boutFrame(isForward),forwardWeight(isForward), ...
    [nBeh 1],@sum,0);
forwardBeh = causalConvolutionFFT(forwardImpulse,kEvent);

% Vigor is a continuous signal. A Riemann-sum approximation to continuous
% convolution requires k/fps. Median centering avoids a spurious startup
% transient and redundancy with the intercept; it does not alter slopes in
% a model with an intercept away from the boundary.
if ~isfield(beh,'vigor') || isempty(beh.vigor)
    if P.vigor.requireField; error('Behavior file lacks continuous vigor.');
    else; vigorRaw=zeros(nBeh,1); vigorSource='zero fallback'; end
else
    vigorRaw=trimOrPad(double(beh.vigor(:)),nBeh);
    vigorSource='pass2FileResult.vigor';
end
missingVigor = mean(~isfinite(vigorRaw));
assert(missingVigor<=P.vigor.maximumMissingFraction, ...
    'Vigor missing fraction %.3f exceeds configured maximum %.3f.', ...
    missingVigor,P.vigor.maximumMissingFraction);
vigorRaw = fillFiniteLinear(vigorRaw);
if P.vigor.centerBeforeFiltering
    vigorCenter=median(vigorRaw,'omitnan');
    vigorForFilter=vigorRaw-vigorCenter;
else
    vigorCenter=0; vigorForFilter=vigorRaw;
end
vigorBeh = causalConvolutionFFT(vigorForFilter,kEvent./fps);

% Sample every predictor at the exact calcium timestamps.
B = struct();
B.ahv = interp1(tBeh,ahvBeh,tCa(:),'linear',NaN);
B.forward = interp1(tBeh,forwardBeh,tCa(:),'linear',NaN);
B.vigor = interp1(tBeh,vigorBeh,tCa(:),'linear',NaN);

D = struct();
D.fpsBehavior=fps; D.nBehaviorFrames=nBeh; D.nBouts=nB;
D.nValidBouts=sum(validBout); D.nForwardBouts=sum(isForward);
D.forwardBoutFraction=safeDivide(D.nForwardBouts,D.nValidBouts);
D.dthetaSource=dthetaSource; D.behaviorTimeSource=timeSource;
D.forwardDefinition=forwardDefinition;
D.forwardWeightSource=forwardWeightSource;
D.nForwardAmp70Replacements=nMissingAmp;
D.vigorSource=vigorSource; D.vigorMissingFraction=missingVigor;
D.vigorCenterBeforeFiltering=vigorCenter;
D.kernelTauSeconds=P.kernel.tauSeconds;
D.kernelCutoffTau=P.kernel.cutoffTau;
D.totalPositiveRotations=sum(max(dtheta(validBout),0),'omitnan')/(2*pi);
D.totalNegativeRotations=sum(max(-dtheta(validBout),0),'omitnan')/(2*pi);
D.minDirectionalRotations=min(D.totalPositiveRotations,D.totalNegativeRotations);
D.ahvRangeBehavior=[min(ahvBeh) max(ahvBeh)];
D.forwardRangeBehavior=[min(forwardBeh) max(forwardBeh)];
D.vigorFilteredRangeBehavior=[min(vigorBeh) max(vigorBeh)];
D.nCalciumSamplesAHV=sum(isfinite(B.ahv));
D.nCalciumSamplesForward=sum(isfinite(B.forward));
D.nCalciumSamplesVigor=sum(isfinite(B.vigor));
D.corrAbsAHV_AHV=safeCorr(abs(B.ahv),B.ahv);
D.corrForwardVigor=safeCorr(B.forward,B.vigor);
D.corrAbsAHVForward=safeCorr(abs(B.ahv),B.forward);
D.corrAbsAHVVigor=safeCorr(abs(B.ahv),B.vigor);
end

function n = inferBehaviorLength(beh)
n = NaN;
priority={'vigor','heading_est','tail_angle'};
for k=1:numel(priority)
    if isfield(beh,priority{k}) && ~isempty(beh.(priority{k}))
        n=numel(beh.(priority{k})); return;
    end
end
if isfield(beh,'nFrames') && isfinite(double(beh.nFrames))
    n=round(double(beh.nFrames)); n=n(1);
end
end

function [t,source] = chooseBehaviorTime(beh,n,fps,tCa,mode)
local=(0:n-1)'./fps;
startFrame=1;
if isfield(beh,'startFrame') && isfinite(double(beh.startFrame))
    startFrame=round(double(beh.startFrame)); startFrame=startFrame(1);
end
fromStart=((startFrame:startFrame+n-1)'-1)./fps;
explicit=[];
if isfield(beh,'behaviorTimeSec') && numel(beh.behaviorTimeSec)==n
    explicit=double(beh.behaviorTimeSec(:));
end

switch lower(mode)
    case 'local'; t=local; source='local: (frame-1)/fps';
    case 'startframe'; t=fromStart; source='startFrame-adjusted';
    case 'behaviortimesec'
        assert(~isempty(explicit),'behaviorTimeSec requested but unavailable.');
        t=explicit; source='pass2FileResult.behaviorTimeSec';
    case 'auto'
        candidates={local,fromStart}; names={'local: (frame-1)/fps','startFrame-adjusted'};
        if ~isempty(explicit); candidates{end+1}=explicit; names{end+1}='behaviorTimeSec'; end
        best=1; bestCoverage=-Inf; bestStart=Inf;
        for k=1:numel(candidates)
            q=candidates{k}; finiteQ=q(isfinite(q));
            if numel(finiteQ)<2; continue; end
            coverage=sum(isfinite(tCa) & tCa>=min(finiteQ) & tCa<=max(finiteQ));
            startDistance=abs(min(tCa(isfinite(tCa)))-min(finiteQ));
            if coverage>bestCoverage || (coverage==bestCoverage && startDistance<bestStart)
                best=k; bestCoverage=coverage; bestStart=startDistance;
            end
        end
        t=candidates{best}; source=['auto -> ' names{best}];
    otherwise; error('Unknown P.behaviorTimeMode: %s',mode);
end
assert(numel(t)==n && all(diff(t(isfinite(t)))>0),'Invalid behavior time vector.');
end

function [isForward,source] = identifyForwardBouts(beh,dtheta,nB,P)
isForward=false(nB,1); source='';
if P.forward.useLIThresholdWhenAvailable && isfield(beh,'LI') && ...
        numel(beh.LI)>=nB && isfield(beh,'LI_threshold') && ...
        ~isempty(beh.LI_threshold) && isfinite(double(beh.LI_threshold(1)))
    li=double(beh.LI(1:nB)); li=li(:);
    thr=abs(double(beh.LI_threshold(1)));
    isForward=isfinite(li) & abs(li)<=thr;
    source=sprintf('|LI| <= LI_threshold (%.6g)',thr);
else
    thr=deg2rad(P.forward.fallbackAbsDthetaDeg);
    isForward=isfinite(dtheta) & abs(dtheta)<=thr;
    source=sprintf('|dtheta| <= %.3g deg fallback',P.forward.fallbackAbsDthetaDeg);
end
end

function [w,source,nReplaced] = getForwardWeights(beh,isForward,nB,P)
w=ones(nB,1); nReplaced=0;
if ~P.forward.useAmp70Weights
    source='unit event weights (configured)'; return;
end
if isfield(beh,'amp70') && ~isempty(beh.amp70)
    amp=NaN(nB,1); m=min(nB,numel(beh.amp70));
    amp(1:m)=double(beh.amp70(1:m));
    finiteForward=isForward & isfinite(amp);
    if any(finiteForward)
        replacement=median(amp(finiteForward),'omitnan');
        missing=isForward & ~isfinite(amp);
        if any(missing)
            assert(P.forward.replaceMissingAmp70WithMedian, ...
                'Some forward bouts have missing amp70.');
            amp(missing)=replacement; nReplaced=sum(missing);
        end
        w=amp; w(~isfinite(w))=1;
        source='amp70 (missing forward values replaced by forward median)';
        if nReplaced==0; source='amp70'; end
        return;
    end
end
assert(P.forward.allowUnitWeightFallback, ...
    'amp70 unavailable for forward bouts and unit fallback is disabled.');
w=ones(nB,1); source='unit event weights: amp70 unavailable';
warning('amp70 unavailable for forward bouts; using unit event weights.');
end

function y = causalConvolutionFFT(x,k)
x=double(x(:)); k=double(k(:));
n=numel(x); m=numel(k); nFFT=2^nextpow2(n+m-1);
yFull=real(ifft(fft(x,nFFT).*fft(k,nFFT)));
y=yFull(1:n);
end

function [startFrame,endFrame] = parseSwimFrames(swimFrames,nBouts)
if size(swimFrames,1)==2
    startFrame=swimFrames(1,:)'; endFrame=swimFrames(2,:)';
elseif size(swimFrames,2)==2
    startFrame=swimFrames(:,1); endFrame=swimFrames(:,2);
else
    error('swimFrames must be 2 x nBouts or nBouts x 2.');
end
if nargin>1
    n=min([nBouts,numel(startFrame),numel(endFrame)]);
    startFrame=startFrame(1:n); endFrame=endFrame(1:n);
end
end

function x = trimOrPad(x,n)
x=x(:);
if numel(x)>n; x=x(1:n);
elseif numel(x)<n; x(end+1:n,1)=NaN; end
end

function x = fillFiniteLinear(x)
x=double(x(:)); ok=isfinite(x);
assert(any(ok),'Continuous predictor has no finite values.');
if all(ok); return; end
idx=(1:numel(x))';
if sum(ok)==1; x(~ok)=x(ok); return; end
x(~ok)=interp1(idx(ok),x(ok),idx(~ok),'linear','extrap');
end

function [C,VIF,conditionZ] = predictorDiagnostics(X)
ok=all(isfinite(X),2); X=X(ok,:);
if size(X,1)<size(X,2)+2
    C=NaN(size(X,2)); VIF=NaN(1,size(X,2)); conditionZ=NaN; return;
end
C=corrcoef(X);
Z=zscoreLocal(X);
conditionZ=cond([Z,ones(size(Z,1),1)]);
VIF=NaN(1,size(X,2));
for j=1:size(X,2)
    y=Z(:,j); other=setdiff(1:size(X,2),j);
    M=[Z(:,other),ones(size(Z,1),1)];
    beta=solveOLS(M,y); pred=M*beta;
    r2=localR2(y,pred);
    if isfinite(r2) && 1-r2>eps; VIF(j)=1/(1-r2);
    elseif isfinite(r2); VIF(j)=Inf; end
end
end

function Fit = fitAllNeurons(Y,phase,ahv,forward,vigor,sourcePrefRad,P)
% Augmented predictor order: cosPhi, sinPhi, absAHV, AHV, F, V, intercept.
[~,N]=size(Y); nanN=NaN(N,1);
fields={'a','d','b0','b1','b2','bF','bV','c0', ...
    'seA','seD','seB0','sePrefRad','sePrefDeg','seB1','seB2','seBF','seBV','seC0', ...
    'pA','pD','pB1','pB2','pBF','pBV','pC0','fHDJoint','pHDJoint', ...
    'fAddedBehavior','pAddedBehavior','partialR2HD','partialR2AddedBehavior', ...
    'prefRad','prefDeg','prefShiftRad','prefShiftDeg', ...
    'aStd','dStd','b0Std','b1Std','b2Std','bFStd','bVStd', ...
    'epsilonStd','epsilonRMS','r2','adjR2','aic','bic', ...
    'phaseOnlyA','phaseOnlyD','phaseOnlyB0','phaseOnlyC0', ...
    'phaseOnlySeA','phaseOnlySeD','phaseOnlySeB0','phaseOnlySePrefRad','phaseOnlySePrefDeg', ...
    'phaseOnlyPA','phaseOnlyPD','phaseOnlyFJoint','phaseOnlyPJoint','phaseOnlyPartialR2', ...
    'phaseOnlyPrefRad','phaseOnlyPrefDeg', ...
    'phaseOnlyAStd','phaseOnlyDStd','phaseOnlyB0Std', ...
    'phaseOnlyR2','phaseOnlyAdjR2','phaseOnlyAIC','phaseOnlyBIC', ...
    'originalA','originalD','originalB0','originalB1','originalB2','originalC0', ...
    'originalPrefRad','originalPrefDeg', ...
    'originalAStd','originalDStd','originalB0Std','originalB1Std','originalB2Std', ...
    'originalR2','originalAdjR2','originalAIC','originalBIC', ...
    'behaviorB1','behaviorB2','behaviorBF','behaviorBV','behaviorC0', ...
    'behaviorSeB1','behaviorSeB2','behaviorSeBF','behaviorSeBV','behaviorSeC0', ...
    'behaviorPB1','behaviorPB2','behaviorPBF','behaviorPBV','behaviorPC0', ...
    'behaviorB1Std','behaviorB2Std','behaviorBFStd','behaviorBVStd', ...
    'behaviorR2','behaviorAdjR2','behaviorAIC','behaviorBIC', ...
    'behaviorEpsilonStd','behaviorEpsilonRMS', ...
    'fAHVBeyondPhase','pAHVBeyondPhase','partialR2AHVBeyondPhase', ...
    'fAllBehaviorBeyondPhase','pAllBehaviorBeyondPhase','partialR2AllBehaviorBeyondPhase', ...
    'deltaB0FromControls','deltaB1FromControls','deltaB2FromControls', ...
    'cvR2PhaseOnly','cvR2Original','cvR2Behavior','cvR2Full','cvR2Direct', ...
    'cvUniquePhase','cvDeltaR2Phase','cvUniqueAddedBehavior', ...
    'cvDeltaR2AddedBehavior','cvUniqueAHVBeyondPhase','cvDeltaR2AHVBeyondPhase', ...
    'cvUniqueAllBehaviorBeyondPhase','cvDeltaR2AllBehaviorBeyondPhase', ...
    'cvSSEPhaseOnly','cvSSEOriginal','cvSSEBehavior','cvSSEFull', ...
    'designRank','conditionNumberRaw','conditionNumberStandardized', ...
    'vifCosPhi','vifSinPhi','vifAbsAHV','vifAHV','vifForward','vifVigor','maxVIF'};
Fit=struct();
for k=1:numel(fields); Fit.(fields{k})=nanN; end
% Cell arrays are necessary because the number of usable held-out frames
% can differ between neurons. single precision limits the MAT-file size.
Fit.cvObserved=cell(N,1);
Fit.cvPredPhaseOnly=cell(N,1);
Fit.cvPredOriginal=cell(N,1);
Fit.cvPredBehavior=cell(N,1);
Fit.cvPredFull=cell(N,1);
Fit.nFitSamples=zeros(N,1);
Fit.modelFormula=['zActivity = a*cos(phi) + d*sin(phi) + b1*abs(AHV) + ' ...
    'b2*AHV + bF*F + bV*V + c0 + epsilon'];
Fit.coefficientOrder={'a','d','b1','b2','bF','bV','c0'};
Fit.originalModelFormula=['zActivity = a*cos(phi) + d*sin(phi) + ' ...
    'b1*abs(AHV) + b2*AHV + c0 + epsilon'];
Fit.phaseOnlyModelFormula=['zActivity = a*cos(phi) + d*sin(phi) + ' ...
    'c0 + epsilon'];
Fit.phaseOnlyCoefficientOrder={'phaseOnlyA','phaseOnlyD','phaseOnlyC0'};
Fit.behaviorModelFormula=['zActivity = b1*abs(AHV) + b2*AHV + ' ...
    'bF*F + bV*V + c0 + epsilon'];
Fit.behaviorCoefficientOrder={'behaviorB1','behaviorB2','behaviorBF', ...
    'behaviorBV','behaviorC0'};
Fit.augmentedModelFormula=Fit.modelFormula;
Fit.augmentedCoefficientOrder=Fit.coefficientOrder;

phase=phase(:); ahv=ahv(:); forward=forward(:); vigor=vigor(:);
sourcePrefRad=sourcePrefRad(:);
base=[cos(phase),sin(phase),abs(ahv),ahv,forward,vigor];

for i=1:N
    y=Y(:,i);
    Xall=[base,ones(numel(y),1)];
    ok=all(isfinite(Xall),2) & isfinite(y);
    if sum(ok)<P.minSamplesForNeuronFit; continue; end
    X=Xall(ok,:); yy=y(ok); n=size(X,1); kFull=size(X,2); dof=n-kFull;
    if dof<=0; continue; end

    beta=solveOLS(X,yy); pred=X*beta; resid=yy-pred;
    yz=zscoreLocal(yy);
    sseFull=sum(resid.^2); mse=sseFull/dof;
    covBeta=mse*pinv(X'*X); se=sqrt(max(0,diag(covBeta)));
    tStat=beta./se; p=2*(1-tcdf(abs(tStat),dof));

    a=beta(1); d=beta(2); b0=hypot(a,d); pref=atan2(d,a);
    Fit.a(i)=a; Fit.d(i)=d; Fit.b0(i)=b0;
    Fit.b1(i)=beta(3); Fit.b2(i)=beta(4); Fit.bF(i)=beta(5);
    Fit.bV(i)=beta(6); Fit.c0(i)=beta(7);
    Fit.seA(i)=se(1); Fit.seD(i)=se(2); Fit.seB1(i)=se(3);
    Fit.seB2(i)=se(4); Fit.seBF(i)=se(5); Fit.seBV(i)=se(6); Fit.seC0(i)=se(7);
    Fit.pA(i)=p(1); Fit.pD(i)=p(2); Fit.pB1(i)=p(3); Fit.pB2(i)=p(4);
    Fit.pBF(i)=p(5); Fit.pBV(i)=p(6); Fit.pC0(i)=p(7);
    Fit.prefRad(i)=wrapToPiLocal(pref); Fit.prefDeg(i)=rad2deg(Fit.prefRad(i));
    if isfinite(sourcePrefRad(i))
        Fit.prefShiftRad(i)=wrapToPiLocal(pref-sourcePrefRad(i));
        Fit.prefShiftDeg(i)=rad2deg(Fit.prefShiftRad(i));
    end
    if isfinite(b0) && b0>sqrt(eps)
        covAD=covBeta(1:2,1:2);
        gB=[a;d]/b0; gP=[-d;a]/(b0^2);
        Fit.seB0(i)=sqrt(max(0,gB'*covAD*gB));
        Fit.sePrefRad(i)=sqrt(max(0,gP'*covAD*gP));
        Fit.sePrefDeg(i)=rad2deg(Fit.sePrefRad(i));
    end

    % Behavior-only coefficients fitted on the identical samples. This is a
    % complete saved fit, not merely an internal CV design matrix.
    Xbehavior=X(:,3:7); % |AHV|, AHV, F, V, intercept
    betaB=solveOLS(Xbehavior,yy); predB=Xbehavior*betaB;
    residB=yy-predB; sseBehavior=sum(residB.^2);
    kBehavior=size(Xbehavior,2); dofBehavior=n-kBehavior;
    Fit.behaviorB1(i)=betaB(1); Fit.behaviorB2(i)=betaB(2);
    Fit.behaviorBF(i)=betaB(3); Fit.behaviorBV(i)=betaB(4);
    Fit.behaviorC0(i)=betaB(5);
    Fit.behaviorR2(i)=localR2(yy,predB);
    Fit.behaviorAdjR2(i)=adjustedR2Local(Fit.behaviorR2(i),n,kBehavior);
    [Fit.behaviorAIC(i),Fit.behaviorBIC(i)]= ...
        informationCriteriaLocal(sseBehavior,n,kBehavior);
    Fit.behaviorEpsilonStd(i)=std(residB,0,'omitnan');
    Fit.behaviorEpsilonRMS(i)=sqrt(mean(residB.^2,'omitnan'));
    if dofBehavior>0
        mseB=sseBehavior/dofBehavior;
        covBetaB=mseB*pinv(Xbehavior'*Xbehavior);
        seB=sqrt(max(0,diag(covBetaB)));
        tB=betaB./seB;
        pB=2*(1-tcdf(abs(tB),dofBehavior));
        Fit.behaviorSeB1(i)=seB(1); Fit.behaviorSeB2(i)=seB(2);
        Fit.behaviorSeBF(i)=seB(3); Fit.behaviorSeBV(i)=seB(4);
        Fit.behaviorSeC0(i)=seB(5);
        Fit.behaviorPB1(i)=pB(1); Fit.behaviorPB2(i)=pB(2);
        Fit.behaviorPBF(i)=pB(3); Fit.behaviorPBV(i)=pB(4);
        Fit.behaviorPC0(i)=pB(5);
    end
    ZB=zscoreLocal(Xbehavior(:,1:4));
    betaBZ=solveOLS([ZB,ones(n,1)],yz);
    Fit.behaviorB1Std(i)=betaBZ(1); Fit.behaviorB2Std(i)=betaBZ(2);
    Fit.behaviorBFStd(i)=betaBZ(3); Fit.behaviorBVStd(i)=betaBZ(4);

    % Nested in-sample phase test: augmented versus behavior-only model.
    [Fit.fHDJoint(i),Fit.pHDJoint(i),Fit.partialR2HD(i)] = ...
        nestedFTest(yy,X,sseFull,dof,Xbehavior);

    % Added-behavior test: augmented versus the exactly matched original model.
    Xoriginal=X(:,[1:4 7]);
    [Fit.fAddedBehavior(i),Fit.pAddedBehavior(i),Fit.partialR2AddedBehavior(i)] = ...
        nestedFTest(yy,X,sseFull,dof,Xoriginal);

    % Standardized coefficients are diagnostics. b0/pref above always come
    % from raw cos/sin coefficients so preferred phase remains well-defined.
    Z=zscoreLocal(X(:,1:6));
    betaZ=solveOLS([Z,ones(n,1)],yz);
    Fit.aStd(i)=betaZ(1); Fit.dStd(i)=betaZ(2);
    Fit.b0Std(i)=hypot(betaZ(1),betaZ(2));
    Fit.b1Std(i)=betaZ(3); Fit.b2Std(i)=betaZ(4);
    Fit.bFStd(i)=betaZ(5); Fit.bVStd(i)=betaZ(6);

    % Phase-only coefficients fitted on the identical samples.
    Xphase=X(:,[1 2 7]);
    betaP=solveOLS(Xphase,yy); predP=Xphase*betaP;
    residP=yy-predP; ssePhase=sum(residP.^2);
    kPhase=size(Xphase,2); dofPhase=n-kPhase;
    Fit.phaseOnlyA(i)=betaP(1); Fit.phaseOnlyD(i)=betaP(2);
    Fit.phaseOnlyB0(i)=hypot(betaP(1),betaP(2));
    Fit.phaseOnlyC0(i)=betaP(3);
    Fit.phaseOnlyPrefRad(i)=wrapToPiLocal(atan2(betaP(2),betaP(1)));
    Fit.phaseOnlyPrefDeg(i)=rad2deg(Fit.phaseOnlyPrefRad(i));
    Fit.phaseOnlyR2(i)=localR2(yy,predP);
    Fit.phaseOnlyAdjR2(i)=adjustedR2Local(Fit.phaseOnlyR2(i),n,kPhase);
    [Fit.phaseOnlyAIC(i),Fit.phaseOnlyBIC(i)]= ...
        informationCriteriaLocal(ssePhase,n,kPhase);
    if dofPhase>0
        mseP=ssePhase/dofPhase;
        covBetaP=mseP*pinv(Xphase'*Xphase);
        seP=sqrt(max(0,diag(covBetaP)));
        tP=betaP./seP;
        pP=2*(1-tcdf(abs(tP),dofPhase));
        Fit.phaseOnlySeA(i)=seP(1); Fit.phaseOnlySeD(i)=seP(2);
        Fit.phaseOnlyPA(i)=pP(1); Fit.phaseOnlyPD(i)=pP(2);
        if Fit.phaseOnlyB0(i)>sqrt(eps)
            covADP=covBetaP(1:2,1:2);
            gBP=betaP(1:2)/Fit.phaseOnlyB0(i);
            gPP=[-betaP(2);betaP(1)]/(Fit.phaseOnlyB0(i)^2);
            Fit.phaseOnlySeB0(i)=sqrt(max(0,gBP'*covADP*gBP));
            Fit.phaseOnlySePrefRad(i)=sqrt(max(0,gPP'*covADP*gPP));
            Fit.phaseOnlySePrefDeg(i)=rad2deg(Fit.phaseOnlySePrefRad(i));
        end
        [Fit.phaseOnlyFJoint(i),Fit.phaseOnlyPJoint(i),Fit.phaseOnlyPartialR2(i)] = ...
            nestedFTest(yy,Xphase,ssePhase,dofPhase,ones(n,1));
    end
    ZP=zscoreLocal(Xphase(:,1:2));
    betaPZ=solveOLS([ZP,ones(n,1)],yz);
    Fit.phaseOnlyAStd(i)=betaPZ(1); Fit.phaseOnlyDStd(i)=betaPZ(2);
    Fit.phaseOnlyB0Std(i)=hypot(betaPZ(1),betaPZ(2));

    % Original coefficients fitted on identical samples.
    betaO=solveOLS(Xoriginal,yy); predO=Xoriginal*betaO;
    residO=yy-predO; sseOriginal=sum(residO.^2);
    dofOriginal=n-size(Xoriginal,2);
    Fit.originalA(i)=betaO(1); Fit.originalD(i)=betaO(2);
    Fit.originalB0(i)=hypot(betaO(1),betaO(2));
    Fit.originalB1(i)=betaO(3); Fit.originalB2(i)=betaO(4); Fit.originalC0(i)=betaO(5);
    Fit.originalPrefRad(i)=wrapToPiLocal(atan2(betaO(2),betaO(1)));
    Fit.originalPrefDeg(i)=rad2deg(Fit.originalPrefRad(i));
    Fit.originalR2(i)=localR2(yy,predO);
    Fit.originalAdjR2(i)=adjustedR2Local(Fit.originalR2(i),n,size(Xoriginal,2));
    [Fit.originalAIC(i),Fit.originalBIC(i)]= ...
        informationCriteriaLocal(sseOriginal,n,size(Xoriginal,2));
    ZO=zscoreLocal(Xoriginal(:,1:4));
    betaOZ=solveOLS([ZO,ones(n,1)],yz);
    Fit.originalAStd(i)=betaOZ(1); Fit.originalDStd(i)=betaOZ(2);
    Fit.originalB0Std(i)=hypot(betaOZ(1),betaOZ(2));
    Fit.originalB1Std(i)=betaOZ(3); Fit.originalB2Std(i)=betaOZ(4);
    [Fit.fAHVBeyondPhase(i),Fit.pAHVBeyondPhase(i),Fit.partialR2AHVBeyondPhase(i)] = ...
        nestedFTest(yy,Xoriginal,sseOriginal,dofOriginal,Xphase);
    [Fit.fAllBehaviorBeyondPhase(i),Fit.pAllBehaviorBeyondPhase(i), ...
        Fit.partialR2AllBehaviorBeyondPhase(i)] = ...
        nestedFTest(yy,X,sseFull,dof,Xphase);
    Fit.deltaB0FromControls(i)=Fit.b0(i)-Fit.originalB0(i);
    Fit.deltaB1FromControls(i)=Fit.b1(i)-Fit.originalB1(i);
    Fit.deltaB2FromControls(i)=Fit.b2(i)-Fit.originalB2(i);

    Fit.epsilonStd(i)=std(resid,0,'omitnan');
    Fit.epsilonRMS(i)=sqrt(mean(resid.^2,'omitnan'));
    Fit.r2(i)=localR2(yy,pred);
    Fit.adjR2(i)=adjustedR2Local(Fit.r2(i),n,kFull);
    [Fit.aic(i),Fit.bic(i)]=informationCriteriaLocal(sseFull,n,kFull);
    Fit.designRank(i)=rank(X); Fit.conditionNumberRaw(i)=cond(X);
    [~,vif,conditionZ]=predictorDiagnostics(X(:,1:6));
    Fit.conditionNumberStandardized(i)=conditionZ;
    Fit.vifCosPhi(i)=vif(1); Fit.vifSinPhi(i)=vif(2);
    Fit.vifAbsAHV(i)=vif(3); Fit.vifAHV(i)=vif(4);
    Fit.vifForward(i)=vif(5); Fit.vifVigor(i)=vif(6);
    Fit.maxVIF(i)=max(vif,[],'omitnan'); Fit.nFitSamples(i)=n;

    if P.cv.enabled
        CV=blockedCrossValidation(yy,Xphase,Xoriginal,Xbehavior,X,P);
        Fit.cvR2PhaseOnly(i)=CV.r2PhaseOnly;
        Fit.cvR2Original(i)=CV.r2Original;
        Fit.cvR2Behavior(i)=CV.r2Behavior;
        Fit.cvR2Full(i)=CV.r2Full;
        Fit.cvR2Direct(i)=CV.r2Full; % compatibility alias; this is augmented full
        Fit.cvSSEPhaseOnly(i)=CV.ssePhaseOnly;
        Fit.cvSSEOriginal(i)=CV.sseOriginal;
        Fit.cvSSEBehavior(i)=CV.sseBehavior;
        Fit.cvSSEFull(i)=CV.sseFull;
        Fit.cvUniquePhase(i)=safeOneMinusRatio(CV.sseFull,CV.sseBehavior);
        Fit.cvDeltaR2Phase(i)=CV.r2Full-CV.r2Behavior;
        Fit.cvUniqueAddedBehavior(i)=safeOneMinusRatio(CV.sseFull,CV.sseOriginal);
        Fit.cvDeltaR2AddedBehavior(i)=CV.r2Full-CV.r2Original;
        Fit.cvUniqueAHVBeyondPhase(i)=safeOneMinusRatio(CV.sseOriginal,CV.ssePhaseOnly);
        Fit.cvDeltaR2AHVBeyondPhase(i)=CV.r2Original-CV.r2PhaseOnly;
        Fit.cvUniqueAllBehaviorBeyondPhase(i)=safeOneMinusRatio(CV.sseFull,CV.ssePhaseOnly);
        Fit.cvDeltaR2AllBehaviorBeyondPhase(i)=CV.r2Full-CV.r2PhaseOnly;
        if P.cv.saveOOFTraces
            Fit.cvObserved{i}=single(CV.observed);
            Fit.cvPredPhaseOnly{i}=single(CV.predPhaseOnly);
            Fit.cvPredOriginal{i}=single(CV.predOriginal);
            Fit.cvPredBehavior{i}=single(CV.predBehavior);
            Fit.cvPredFull{i}=single(CV.predFull);
        end
    end
end
end

function [f,p,partialR2] = nestedFTest(y,Xfull,sseFull,dofFull,Xreduced)
f=NaN; p=NaN; partialR2=NaN;
betaR=solveOLS(Xreduced,y); residR=y-Xreduced*betaR; sseR=sum(residR.^2);
q=rank(Xfull)-rank(Xreduced);
if q>0 && dofFull>0 && sseFull>0 && sseR>=sseFull-1e-10
    f=max(0,((sseR-sseFull)/q)/(sseFull/dofFull));
    p=1-fcdf(f,q,dofFull);
    if sseR>0; partialR2=max(0,(sseR-sseFull)/sseR); end
end
end

function CV = blockedCrossValidation(y,Xphase,Xoriginal,Xbehavior,Xfull,P)
n=numel(y); nFolds=min(P.cv.nBlockedFolds,floor(n/P.cv.minTestSamplesPerFold));
CV=struct('r2PhaseOnly',NaN,'r2Original',NaN,'r2Behavior',NaN,'r2Full',NaN, ...
    'ssePhaseOnly',NaN,'sseOriginal',NaN,'sseBehavior',NaN,'sseFull',NaN, ...
    'observed',[],'predPhaseOnly',[],'predOriginal',[], ...
    'predBehavior',[],'predFull',[]);
if nFolds<2; return; end
edges=round(linspace(0,n,nFolds+1));
pP=NaN(n,1); pO=NaN(n,1); pB=NaN(n,1); pF=NaN(n,1);
for f=1:nFolds
    test=(edges(f)+1):edges(f+1); train=true(n,1); train(test)=false;
    if numel(test)<P.cv.minTestSamplesPerFold || sum(train)<=size(Xfull,2); continue; end
    pP(test)=Xphase(test,:)*solveOLS(Xphase(train,:),y(train));
    pO(test)=Xoriginal(test,:)*solveOLS(Xoriginal(train,:),y(train));
    pB(test)=Xbehavior(test,:)*solveOLS(Xbehavior(train,:),y(train));
    pF(test)=Xfull(test,:)*solveOLS(Xfull(train,:),y(train));
end
valid=isfinite(y) & isfinite(pP) & isfinite(pO) & isfinite(pB) & isfinite(pF);
if sum(valid)<P.minSamplesForNeuronFit; return; end
yy=y(valid); pP=pP(valid); pO=pO(valid); pB=pB(valid); pF=pF(valid);
% Preserve only the common, aligned OOF support. These are the same samples
% used below for all four CV-R2 values and therefore permit paired ROC/AUC.
CV.observed=yy;
CV.predPhaseOnly=pP;
CV.predOriginal=pO;
CV.predBehavior=pB;
CV.predFull=pF;
CV.ssePhaseOnly=sum((yy-pP).^2); CV.sseOriginal=sum((yy-pO).^2);
CV.sseBehavior=sum((yy-pB).^2); CV.sseFull=sum((yy-pF).^2);
sst=sum((yy-mean(yy)).^2);
if sst>0
    CV.r2PhaseOnly=1-CV.ssePhaseOnly/sst;
    CV.r2Original=1-CV.sseOriginal/sst;
    CV.r2Behavior=1-CV.sseBehavior/sst;
    CV.r2Full=1-CV.sseFull/sst;
end
end

function beta = solveOLS(X,y)
if rank(X)<size(X,2); beta=pinv(X)*y; else; beta=X\y; end
end

function value = safeOneMinusRatio(numer,denom)
if isfinite(numer) && isfinite(denom) && denom>0; value=1-numer/denom;
else; value=NaN; end
end

function [labels,ahs,info] = classifyAugmented(Fit,P)
N=numel(Fit.b2); labels=repmat({'NotFit'},N,1); ahs=repmat({'AHS_not_fit'},N,1);
valid=isfinite(Fit.b2) & isfinite(Fit.pB2) & Fit.nFitSamples>=P.minSamplesForNeuronFit;
nValid=sum(valid);
switch lower(P.class.bonferroniMode)
    case 'per_fish'; nTests=max(1,nValid);
    case 'fixed'
        nTests=round(P.class.fixedNTests);
        assert(nTests>=nValid,['P.class.fixedNTests=%d is smaller than %d valid neurons. ' ...
            'Use the common maximum candidate count.'],nTests,nValid);
    otherwise; error('Unknown Bonferroni mode.');
end
thr=P.class.alpha/nTests; sig=valid & Fit.pB2<thr;
labels(valid & ~sig)={'Symmetric'};
if P.ahvPositiveIsCCW
    labels(sig & Fit.b2>0)={'CCW'}; labels(sig & Fit.b2<0)={'CW'};
else
    labels(sig & Fit.b2>0)={'CW'}; labels(sig & Fit.b2<0)={'CCW'};
end
sigA=valid & isfinite(Fit.pB1) & Fit.pB1<thr;
ahs(valid & ~sigA)={'AHS_not_significant'};
ahs(sigA & Fit.b1>0)={'AHS_positive'}; ahs(sigA & Fit.b1<0)={'AHS_negative'};
info=struct('alpha',P.class.alpha,'alphaBonferroni',thr, ...
    'bonferroniMode',P.class.bonferroniMode,'nBonferroniTests',nTests, ...
    'nValidNeuronsThisFish',nValid, ...
    'classificationRule','Sign/significance of augmented-model b2', ...
    'pValueCaveat',['Conventional OLS p-values do not correct for temporal ' ...
    'autocorrelation; use continuous coefficients and blocked CV for primary inference.']);
end

function T = makeFitTable(session,candidateOrder,neuronID,F,classLabel,ahsLabel,info,fishVIF,fishConditionZ)
N=numel(neuronID);
T=table(session,candidateOrder,neuronID(:), ...
    F.prefRad,F.prefDeg,F.prefShiftRad,F.prefShiftDeg, ...
    F.a,F.d,F.b0,F.b1,F.b2,F.bF,F.bV,F.c0, ...
    F.seA,F.seD,F.seB0,F.sePrefRad,F.sePrefDeg,F.seB1,F.seB2,F.seBF,F.seBV,F.seC0, ...
    F.pA,F.pD,F.pHDJoint,F.fHDJoint,F.pB1,F.pB2,F.pBF,F.pBV,F.pC0, ...
    F.pAddedBehavior,F.fAddedBehavior, ...
    F.aStd,F.dStd,F.b0Std,F.b1Std,F.b2Std,F.bFStd,F.bVStd, ...
    F.partialR2HD,F.partialR2AddedBehavior,F.r2,F.adjR2,F.aic,F.bic, ...
    F.phaseOnlyA,F.phaseOnlyD,F.phaseOnlyB0,F.phaseOnlyC0, ...
    F.phaseOnlySeA,F.phaseOnlySeD,F.phaseOnlySeB0, ...
    F.phaseOnlySePrefRad,F.phaseOnlySePrefDeg, ...
    F.phaseOnlyPA,F.phaseOnlyPD,F.phaseOnlyFJoint,F.phaseOnlyPJoint,F.phaseOnlyPartialR2, ...
    F.phaseOnlyPrefRad,F.phaseOnlyPrefDeg, ...
    F.phaseOnlyAStd,F.phaseOnlyDStd,F.phaseOnlyB0Std, ...
    F.phaseOnlyR2,F.phaseOnlyAdjR2,F.phaseOnlyAIC,F.phaseOnlyBIC, ...
    F.originalA,F.originalD,F.originalB0,F.originalB1,F.originalB2,F.originalC0, ...
    F.originalPrefRad,F.originalPrefDeg, ...
    F.originalAStd,F.originalDStd,F.originalB0Std,F.originalB1Std,F.originalB2Std, ...
    F.originalR2,F.originalAdjR2,F.originalAIC,F.originalBIC, ...
    F.behaviorB1,F.behaviorB2,F.behaviorBF,F.behaviorBV,F.behaviorC0, ...
    F.behaviorSeB1,F.behaviorSeB2,F.behaviorSeBF,F.behaviorSeBV,F.behaviorSeC0, ...
    F.behaviorPB1,F.behaviorPB2,F.behaviorPBF,F.behaviorPBV,F.behaviorPC0, ...
    F.behaviorB1Std,F.behaviorB2Std,F.behaviorBFStd,F.behaviorBVStd, ...
    F.behaviorR2,F.behaviorAdjR2,F.behaviorAIC,F.behaviorBIC, ...
    F.behaviorEpsilonStd,F.behaviorEpsilonRMS, ...
    F.fAHVBeyondPhase,F.pAHVBeyondPhase,F.partialR2AHVBeyondPhase, ...
    F.fAllBehaviorBeyondPhase,F.pAllBehaviorBeyondPhase,F.partialR2AllBehaviorBeyondPhase, ...
    F.deltaB0FromControls,F.deltaB1FromControls,F.deltaB2FromControls, ...
    F.cvR2PhaseOnly,F.cvR2Original,F.cvR2Behavior,F.cvR2Full,F.cvR2Direct, ...
    F.cvUniquePhase,F.cvDeltaR2Phase,F.cvUniqueAddedBehavior,F.cvDeltaR2AddedBehavior, ...
    F.cvUniqueAHVBeyondPhase,F.cvDeltaR2AHVBeyondPhase, ...
    F.cvUniqueAllBehaviorBeyondPhase,F.cvDeltaR2AllBehaviorBeyondPhase, ...
    F.cvSSEPhaseOnly,F.cvSSEOriginal,F.cvSSEBehavior,F.cvSSEFull, ...
    F.epsilonStd,F.epsilonRMS,F.designRank,F.conditionNumberRaw, ...
    F.conditionNumberStandardized,F.vifCosPhi,F.vifSinPhi,F.vifAbsAHV,F.vifAHV, ...
    F.vifForward,F.vifVigor,F.maxVIF,F.nFitSamples,classLabel(:),ahsLabel(:), ...
    'VariableNames',{'session','candidateOrder','neuronID', ...
    'prefRad','prefDeg','prefShiftRad','prefShiftDeg', ...
    'a','d','b0','b1','b2','bF','bV','c0', ...
    'seA','seD','seB0','sePrefRad','sePrefDeg','seB1','seB2','seBF','seBV','seC0', ...
    'pA','pD','pHDJoint','fHDJoint','pB1','pB2','pBF','pBV','pC0', ...
    'pAddedBehavior','fAddedBehavior', ...
    'a_std','d_std','b0_std','b1_std','b2_std','bF_std','bV_std', ...
    'partialR2HD','partialR2AddedBehavior','r2','adjR2','AIC','BIC', ...
    'phaseOnlyA','phaseOnlyD','phaseOnlyB0','phaseOnlyC0', ...
    'phaseOnlySeA','phaseOnlySeD','phaseOnlySeB0', ...
    'phaseOnlySePrefRad','phaseOnlySePrefDeg', ...
    'phaseOnlyPA','phaseOnlyPD','phaseOnlyFJoint','phaseOnlyPJoint','phaseOnlyPartialR2', ...
    'phaseOnlyPrefRad','phaseOnlyPrefDeg', ...
    'phaseOnlyA_std','phaseOnlyD_std','phaseOnlyB0_std', ...
    'phaseOnlyR2','phaseOnlyAdjR2','phaseOnlyAIC','phaseOnlyBIC', ...
    'originalA','originalD','originalB0','originalB1','originalB2','originalC0', ...
    'originalPrefRad','originalPrefDeg', ...
    'originalA_std','originalD_std','originalB0_std','originalB1_std','originalB2_std', ...
    'originalR2','originalAdjR2','originalAIC','originalBIC', ...
    'behaviorB1','behaviorB2','behaviorBF','behaviorBV','behaviorC0', ...
    'behaviorSeB1','behaviorSeB2','behaviorSeBF','behaviorSeBV','behaviorSeC0', ...
    'behaviorPB1','behaviorPB2','behaviorPBF','behaviorPBV','behaviorPC0', ...
    'behaviorB1_std','behaviorB2_std','behaviorBF_std','behaviorBV_std', ...
    'behaviorR2','behaviorAdjR2','behaviorAIC','behaviorBIC', ...
    'behaviorEpsilonStd','behaviorEpsilonRMS', ...
    'fAHVBeyondPhase','pAHVBeyondPhase','partialR2AHVBeyondPhase', ...
    'fAllBehaviorBeyondPhase','pAllBehaviorBeyondPhase','partialR2AllBehaviorBeyondPhase', ...
    'deltaB0FromControls','deltaB1FromControls','deltaB2FromControls', ...
    'cvR2PhaseOnly','cvR2Original','cvR2Behavior','cvR2Full','cvR2Direct', ...
    'cvUniquePhase','cvDeltaR2Phase','cvUniqueAddedBehavior','cvDeltaR2AddedBehavior', ...
    'cvUniqueAHVBeyondPhase','cvDeltaR2AHVBeyondPhase', ...
    'cvUniqueAllBehaviorBeyondPhase','cvDeltaR2AllBehaviorBeyondPhase', ...
    'cvSSEPhaseOnly','cvSSEOriginal','cvSSEBehavior','cvSSEFull', ...
    'epsilonStd','epsilonRMS','designRank','conditionNumberRaw', ...
    'conditionNumberStandardized','vifCosPhi','vifSinPhi','vifAbsAHV','vifAHV', ...
    'vifForward','vifVigor','maxVIF','nFitSamples','class','ahsLabel'});

T.pThresholdAHV=repmat(info.alphaBonferroni,N,1);
T.bonferroniNTests=repmat(info.nBonferroniTests,N,1);
T.nValidNeuronsThisFish=repmat(info.nValidNeuronsThisFish,N,1);
T.fishVIFCosPhi=repmat(fishVIF(1),N,1);
T.fishVIFSinPhi=repmat(fishVIF(2),N,1);
T.fishVIFAbsAHV=repmat(fishVIF(3),N,1);
T.fishVIFAHV=repmat(fishVIF(4),N,1);
T.fishVIFForward=repmat(fishVIF(5),N,1);
T.fishVIFVigor=repmat(fishVIF(6),N,1);
T.fishConditionNumberStandardized=repmat(fishConditionZ,N,1);
% Exact names expected by compare_HD_AHV_models_across_morphs_updated.m.
% Each table row contains one numeric column vector of aligned OOF samples.
T.cvObserved=F.cvObserved;
T.cvPredPhaseOnly=F.cvPredPhaseOnly;
T.cvPredOriginal=F.cvPredOriginal;
T.cvPredBehavior=F.cvPredBehavior;
T.cvPredFull=F.cvPredFull;
% Explicit augmented aliases improve readability while the *Full names are
% retained for compatibility with the existing cross-morph analysis.
T.augmentedR2=T.r2;
T.augmentedAdjR2=T.adjR2;
T.augmentedAIC=T.AIC;
T.augmentedBIC=T.BIC;
T.cvR2Augmented=T.cvR2Full;
T.cvSSEAugmented=T.cvSSEFull;
T.cvPredAugmented=F.cvPredFull;
end

function T = makeFishSummary(session,fitTable,D,C,VIF,conditionZ)
ok=isfinite(fitTable.b0);
cvPair=ok & isfinite(fitTable.cvR2Original) & isfinite(fitTable.cvR2Full);
if any(cvPair)
    fractionFullBetter=mean(fitTable.cvR2Full(cvPair)>fitTable.cvR2Original(cvPair));
else
    fractionFullBetter=NaN;
end
cvPhaseOriginalPair=ok & isfinite(fitTable.cvR2PhaseOnly) & isfinite(fitTable.cvR2Original);
if any(cvPhaseOriginalPair)
    fractionOriginalBetterThanPhaseOnly=mean( ...
        fitTable.cvR2Original(cvPhaseOriginalPair)>fitTable.cvR2PhaseOnly(cvPhaseOriginalPair));
else
    fractionOriginalBetterThanPhaseOnly=NaN;
end
cvPhaseFullPair=ok & isfinite(fitTable.cvR2PhaseOnly) & isfinite(fitTable.cvR2Full);
if any(cvPhaseFullPair)
    fractionFullBetterThanPhaseOnly=mean( ...
        fitTable.cvR2Full(cvPhaseFullPair)>fitTable.cvR2PhaseOnly(cvPhaseFullPair));
else
    fractionFullBetterThanPhaseOnly=NaN;
end
T=table(string(session),sum(ok),D.nValidBouts,D.nForwardBouts,D.forwardBoutFraction, ...
    D.totalPositiveRotations,D.totalNegativeRotations,D.minDirectionalRotations, ...
    string(D.behaviorTimeSource),string(D.forwardDefinition),string(D.forwardWeightSource), ...
    D.vigorMissingFraction,D.corrAbsAHV_AHV,D.corrForwardVigor, ...
    D.corrAbsAHVForward,D.corrAbsAHVVigor,conditionZ,max(VIF,[],'omitnan'), ...
    median(fitTable.phaseOnlyB0(ok),'omitnan'), ...
    median(fitTable.originalB0(ok),'omitnan'),median(fitTable.b0(ok),'omitnan'), ...
    median(fitTable.originalB1(ok),'omitnan'),median(fitTable.b1(ok),'omitnan'), ...
    median(abs(fitTable.originalB2(ok)),'omitnan'),median(abs(fitTable.b2(ok)),'omitnan'), ...
    median(fitTable.b2(ok),'omitnan'),mean(fitTable.b2(ok)>0,'omitnan'), ...
    median(fitTable.bF_std(ok),'omitnan'),median(fitTable.bV_std(ok),'omitnan'), ...
    median(fitTable.cvR2PhaseOnly(ok),'omitnan'), ...
    median(fitTable.cvR2Original(ok),'omitnan'),median(fitTable.cvR2Behavior(ok),'omitnan'), ...
    median(fitTable.cvR2Full(ok),'omitnan'),median(fitTable.cvUniquePhase(ok),'omitnan'), ...
    median(fitTable.cvUniqueAddedBehavior(ok),'omitnan'), ...
    median(fitTable.cvDeltaR2AddedBehavior(ok),'omitnan'), ...
    median(fitTable.cvUniqueAHVBeyondPhase(ok),'omitnan'), ...
    median(fitTable.cvDeltaR2AHVBeyondPhase(ok),'omitnan'), ...
    median(fitTable.cvUniqueAllBehaviorBeyondPhase(ok),'omitnan'), ...
    median(fitTable.cvDeltaR2AllBehaviorBeyondPhase(ok),'omitnan'), ...
    fractionOriginalBetterThanPhaseOnly,fractionFullBetterThanPhaseOnly,fractionFullBetter, ...
    median(fitTable.phaseOnlyPartialR2(ok),'omitnan'), ...
    median(fitTable.partialR2AHVBeyondPhase(ok),'omitnan'), ...
    median(fitTable.partialR2AllBehaviorBeyondPhase(ok),'omitnan'), ...
    median(fitTable.partialR2HD(ok),'omitnan'),median(fitTable.partialR2AddedBehavior(ok),'omitnan'), ...
    'VariableNames',{'session','nFitNeurons','nBouts','nForwardBouts','forwardBoutFraction', ...
    'positiveRotations','negativeRotations','minDirectionalRotations', ...
    'behaviorTimeSource','forwardDefinition','forwardWeightSource','vigorMissingFraction', ...
    'corrAbsAHV_AHV','corrForwardVigor','corrAbsAHVForward','corrAbsAHVVigor', ...
    'conditionNumberStandardized','maxVIF', ...
    'medianPhaseOnlyB0','medianOriginalB0','medianAugmentedB0', ...
    'medianOriginalB1','medianAugmentedB1', ...
    'medianAbsOriginalB2','medianAbsAugmentedB2','medianSignedAugmentedB2', ...
    'fractionAugmentedB2Positive','medianBFStd','medianBVStd', ...
    'medianCvR2PhaseOnly','medianCvR2Original','medianCvR2Behavior','medianCvR2Full', ...
    'medianCvUniquePhase','medianCvUniqueAddedBehavior','medianCvDeltaR2AddedBehavior', ...
    'medianCvUniqueAHVBeyondPhase','medianCvDeltaR2AHVBeyondPhase', ...
    'medianCvUniqueAllBehaviorBeyondPhase','medianCvDeltaR2AllBehaviorBeyondPhase', ...
    'fractionNeuronsOriginalCVBetterThanPhaseOnly', ...
    'fractionNeuronsFullCVBetterThanPhaseOnly', ...
    'fractionNeuronsFullCVBetterThanOriginal', ...
    'medianPhaseOnlyPartialR2','medianPartialR2AHVBeyondPhase', ...
    'medianPartialR2AllBehaviorBeyondPhase','medianPartialR2HD', ...
    'medianPartialR2AddedBehavior'});

% Complete four-model fish-level summaries. These explicit names are used
% by the plotting functions below; the legacy columns above are preserved.
T.medianR2PhaseOnly=median(fitTable.phaseOnlyR2(ok),'omitnan');
T.medianR2Original=median(fitTable.originalR2(ok),'omitnan');
T.medianR2Behavior=median(fitTable.behaviorR2(ok),'omitnan');
T.medianR2Augmented=median(fitTable.augmentedR2(ok),'omitnan');
T.medianAICPhaseOnly=median(fitTable.phaseOnlyAIC(ok),'omitnan');
T.medianAICOriginal=median(fitTable.originalAIC(ok),'omitnan');
T.medianAICBehavior=median(fitTable.behaviorAIC(ok),'omitnan');
T.medianAICAugmented=median(fitTable.augmentedAIC(ok),'omitnan');
T.medianBICPhaseOnly=median(fitTable.phaseOnlyBIC(ok),'omitnan');
T.medianBICOriginal=median(fitTable.originalBIC(ok),'omitnan');
T.medianBICBehavior=median(fitTable.behaviorBIC(ok),'omitnan');
T.medianBICAugmented=median(fitTable.augmentedBIC(ok),'omitnan');
T.medianBehaviorB1Std=median(fitTable.behaviorB1_std(ok),'omitnan');
T.medianBehaviorB2Std=median(fitTable.behaviorB2_std(ok),'omitnan');
T.medianBehaviorBFStd=median(fitTable.behaviorBF_std(ok),'omitnan');
T.medianBehaviorBVStd=median(fitTable.behaviorBV_std(ok),'omitnan');
end

function plotPhaseTuningSelection(sessionName,Q,P,outDir)
fig=figure('Visible',P.figureVisible,'Color','w','Position',[100 100 1350 820]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl);
histogram(ax,Q.infoBits(isfinite(Q.infoBits)),40,'FaceColor',[.55 .55 .55], ...
    'EdgeColor','none'); hold(ax,'on');
histogram(ax,Q.infoBits(Q.isPhaseTuned),40,'FaceColor',[.85 .25 .20], ...
    'EdgeColor','none');
xlabel(ax,'Observed phase information (bits)'); ylabel(ax,'Candidate cells');
legend(ax,{'All candidates','BH-FDR selected'},'Location','best'); grid(ax,'on');

ax=nexttile(tl);
scatter(ax,Q.null95,Q.infoBits,28,Q.isPhaseTuned,'filled'); hold(ax,'on');
finitePair=isfinite(Q.null95)&isfinite(Q.infoBits);
if any(finitePair)
    lo=min([Q.null95(finitePair);Q.infoBits(finitePair)]);
    hi=max([Q.null95(finitePair);Q.infoBits(finitePair)]);
    if hi<=lo; hi=lo+1; end
    plot(ax,[lo hi],[lo hi],'k--','LineWidth',1);
end
xlabel(ax,'95th percentile shuffled information');
ylabel(ax,'Observed information'); grid(ax,'on');

ax=nexttile(tl);
scatter(ax,Q.pShuffle,Q.qFDR,28,Q.isPhaseTuned,'filled');
xline(ax,P.tuning.fdrAlpha,'k:'); yline(ax,P.tuning.fdrAlpha,'k:');
xlabel(ax,'Raw circular-shift p'); ylabel(ax,'BH-FDR q');
xlim(ax,[0 1]); ylim(ax,[0 1]); grid(ax,'on');

ax=nexttile(tl);
selectedPref=Q.prefRad(Q.isPhaseTuned & isfinite(Q.prefRad));
if isempty(selectedPref)
    text(ax,.5,.5,'No phase-tuned cells','HorizontalAlignment','center');
    axis(ax,'off');
else
    histogram(ax,rad2deg(selectedPref),18,'BinLimits',[-180 180], ...
        'FaceColor',[.25 .50 .85],'EdgeColor','none');
    xlim(ax,[-180 180]); xticks(ax,-180:90:180);
    xlabel(ax,'Selected-cell preferred phase (deg)'); ylabel(ax,'Cells'); grid(ax,'on');
end

title(tl,sprintf('%s | %s: ORI\\_V15 phase-tuning selection | %d/%d cells', ...
    P.morphDisplayName,sessionName,Q.nSelected,Q.nCandidates), ...
    'Interpreter','tex');
saveFigureLocal(fig,outDir,[sessionName '_ORI_V15_phase_tuning_selection'],P);
end

function plotRegressorQC(R,P,outDir)
idx=evenlySpacedIndices(numel(R.tCa),P.qc.maxTimePoints);
t=R.tCa(idx)/60;
fig=figure('Visible',P.figureVisible,'Color','w','Position',[100 100 1450 850]);
tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); plot(ax,t,wrapToPiLocal(R.phaseCa(idx)),'k');
ylabel(ax,'phase (rad)'); title(ax,'Network phase'); grid(ax,'on');
ax=nexttile(tl); plot(ax,t,R.ahvCa(idx)*180/pi,'Color',[0.15 0.35 0.85]);
ylabel(ax,'AHV (deg/s)'); title(ax,'Signed AHV'); grid(ax,'on');
ax=nexttile(tl); plot(ax,t,zscoreLocal(R.forwardCa(idx)),'Color',[0.85 0.35 0.10]);
ylabel(ax,'F (z for display)'); title(ax,'Filtered amp70-weighted forward events'); grid(ax,'on');
ax=nexttile(tl); plot(ax,t,zscoreLocal(R.vigorCa(idx)),'Color',[0.20 0.60 0.25]);
ylabel(ax,'V (z for display)'); xlabel(ax,'Time (min)'); title(ax,'Filtered continuous vigor'); grid(ax,'on');
linkaxes(findall(fig,'Type','axes'),'x');
title(tl,sprintf('%s | %s: augmented behavior regressors', ...
    P.morphDisplayName,R.name),'Interpreter','none');
saveFigureLocal(fig,outDir,[R.name '_augmented_behavior_regressor_QC'],P);
end

function plotPredictorCorrelation(R,P,outDir)
fig=figure('Visible',P.figureVisible,'Color','w','Position',[100 100 800 700]);
ax=axes(fig); imagesc(ax,R.predictorCorrelation,[-1 1]); axis(ax,'square');
colormap(ax,blueWhiteRedMap(257)); colorbar(ax); labels=R.predictorNames;
xticks(ax,1:numel(labels)); yticks(ax,1:numel(labels));
xticklabels(ax,labels); yticklabels(ax,labels); xtickangle(ax,35);
for i=1:numel(labels)
    for j=1:numel(labels)
        v=R.predictorCorrelation(i,j);
        if isfinite(v); text(ax,j,i,sprintf('%.2f',v),'HorizontalAlignment','center', ...
                'Color',double(abs(v)>0.55)*[1 1 1]+double(abs(v)<=0.55)*[0 0 0]); end
    end
end
title(ax,sprintf('%s | %s predictor correlations | max VIF %.2f', ...
    P.morphDisplayName,R.name,max(R.predictorVIF,[],'omitnan')), ...
    'Interpreter','none');
saveFigureLocal(fig,outDir,[R.name '_augmented_predictor_correlation'],P);
end

function plotPerFishModelComparison(R,P,outDir)
T=R.fitTable;
fig=figure('Visible',P.figureVisible,'Color','w','Position',[40 60 1850 950]);
tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');

cvR2=[T.cvR2PhaseOnly,T.cvR2Original,T.cvR2Full,T.cvR2Behavior];
inSampleR2=[T.phaseOnlyR2,T.originalR2,T.augmentedR2,T.behaviorR2];
aic=[T.phaseOnlyAIC,T.originalAIC,T.augmentedAIC,T.behaviorAIC];
bic=[T.phaseOnlyBIC,T.originalBIC,T.augmentedBIC,T.behaviorBIC];

ax=nexttile(tl); fourModelPlot(ax,cvR2,'Blocked-CV R^2',true);
ax=nexttile(tl); fourModelPlot(ax,inSampleR2,'In-sample R^2',true);
ax=nexttile(tl); fourModelPlot(ax,aic-aic(:,1),'Delta AIC vs phase-only',true); yline(ax,0,'k:');
ax=nexttile(tl); fourModelPlot(ax,bic-bic(:,1),'Delta BIC vs phase-only',true); yline(ax,0,'k:');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAHVBeyondPhase,40, ...
    [0.55 0.30 0.75],'1 - SSE_{original}/SSE_{phase-only}','AHV beyond phase');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAllBehaviorBeyondPhase,40, ...
    [0.85 0.40 0.15],'1 - SSE_{augmented}/SSE_{phase-only}','All behavior beyond phase');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniquePhase,40, ...
    [0.25 0.25 0.25],'1 - SSE_{augmented}/SSE_{behavior}','Phase beyond all behavior');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAddedBehavior,40, ...
    [0.15 0.55 0.90],'1 - SSE_{augmented}/SSE_{original}','Added F+V beyond original');
title(tl,[P.morphDisplayName ' | ' R.name ...
    ': matched four-model comparison'],'Interpreter','none');
saveFigureLocal(fig,outDir,[R.name '_phase_only_and_augmented_model_comparison'],P);
end

function plotAllFishModelComparison(T,P,outDir)
fig=figure('Visible',P.figureVisible,'Color','w','Position',[40 60 1850 950]);
tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');

cvR2=[T.cvR2PhaseOnly,T.cvR2Original,T.cvR2Full,T.cvR2Behavior];
inSampleR2=[T.phaseOnlyR2,T.originalR2,T.augmentedR2,T.behaviorR2];
aic=[T.phaseOnlyAIC,T.originalAIC,T.augmentedAIC,T.behaviorAIC];
bic=[T.phaseOnlyBIC,T.originalBIC,T.augmentedBIC,T.behaviorBIC];

ax=nexttile(tl); fourModelPlot(ax,cvR2,'Blocked-CV R^2',false);
ax=nexttile(tl); fourModelPlot(ax,inSampleR2,'In-sample R^2',false);
ax=nexttile(tl); fourModelPlot(ax,aic-aic(:,1),'Delta AIC vs phase-only',false); yline(ax,0,'k:');
ax=nexttile(tl); fourModelPlot(ax,bic-bic(:,1),'Delta BIC vs phase-only',false); yline(ax,0,'k:');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAHVBeyondPhase,60, ...
    [0.55 0.30 0.75],'1 - SSE_{original}/SSE_{phase-only}','AHV beyond phase');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAllBehaviorBeyondPhase,60, ...
    [0.85 0.40 0.15],'1 - SSE_{augmented}/SSE_{phase-only}','All behavior beyond phase');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniquePhase,60, ...
    [0.25 0.25 0.25],'1 - SSE_{augmented}/SSE_{behavior}','Phase beyond all behavior');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAddedBehavior,60, ...
    [0.15 0.55 0.90],'1 - SSE_{augmented}/SSE_{original}','Added F+V beyond original');
title(tl,[P.morphDisplayName ...
    ' | all fish: phase-only, original, behavior-only and augmented models']);
saveFigureLocal(fig,outDir,'ALL_FISH_phase_only_and_augmented_model_comparison',P);
end

function plotAllFishCoefficientSummary(T,P,outDir)
vars={'phaseOnlyB0_std','originalB0_std','b0_std', ...
    'behaviorB1_std','behaviorB2_std','behaviorBF_std','behaviorBV_std', ...
    'b1_std','b2_std','bF_std','bV_std', ...
    'cvR2PhaseOnly','cvR2Original','cvR2Behavior','cvR2Full'};
labels={'Phase-only b0 std','Original b0 std','Augmented b0 std', ...
    'Behavior-only |AHV| std','Behavior-only AHV std', ...
    'Behavior-only forward std','Behavior-only vigor std', ...
    'Augmented |AHV| std','Augmented AHV std', ...
    'Augmented forward std','Augmented vigor std', ...
    'Phase-only blocked-CV R^2','Original blocked-CV R^2', ...
    'Behavior-only blocked-CV R^2','Augmented blocked-CV R^2'};
fig=figure('Visible',P.figureVisible,'Color','w','Position',[70 70 1750 1000]);
tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
for k=1:numel(vars)
    ax=nexttile(tl); x=T.(vars{k}); x=x(isfinite(x));
    histogram(ax,x,50,'FaceColor',[0.25 0.25 0.25],'EdgeColor','none');
    xline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,labels{k}); ylabel(ax,'Neurons'); grid(ax,'on');
end
title(tl,[P.morphDisplayName ...
    ' | all fish: four-model standardized coefficients and prediction']);
saveFigureLocal(fig,outDir,'ALL_FISH_phase_only_and_augmented_distributions',P);
end

function plotFishLevelSummary(F,P,outDir)
n=height(F); x=(1:n)';
fig=figure('Visible',P.figureVisible,'Color','w','Position',[80 80 1650 900]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl); pairedFishFour(ax,x,F.medianR2PhaseOnly,F.medianR2Original, ...
    F.medianR2Augmented,F.medianR2Behavior,'Median in-sample R^2');
ax=nexttile(tl); pairedFishFour(ax,x,F.medianCvR2PhaseOnly,F.medianCvR2Original, ...
    F.medianCvR2Full,F.medianCvR2Behavior,'Median blocked-CV R^2');
ax=nexttile(tl); plot(ax,x,F.medianCvUniqueAHVBeyondPhase,'o-','LineWidth',1.2, ...
    'MarkerFaceColor',[0.55 0.30 0.75]); yline(ax,0,'k:');
ylabel(ax,'Median unique AHV CV'); grid(ax,'on'); title(ax,'AHV beyond phase');
ax=nexttile(tl); plot(ax,x,F.medianCvUniqueAllBehaviorBeyondPhase,'o-','LineWidth',1.2, ...
    'MarkerFaceColor',[0.90 0.40 0.15]); yline(ax,0,'k:');
ylabel(ax,'Median unique behavior CV'); grid(ax,'on'); title(ax,'All behavior beyond phase');
ax=nexttile(tl); plot(ax,x,F.medianCvUniquePhase,'o-','LineWidth',1.2, ...
    'MarkerFaceColor',[0.20 0.55 0.85]); yline(ax,0,'k:');
ylabel(ax,'Median unique phase CV'); grid(ax,'on'); title(ax,'Phase beyond all behavior');
ax=nexttile(tl); plot(ax,x,F.medianCvUniqueAddedBehavior,'o-','LineWidth',1.2, ...
    'MarkerFaceColor',[0.90 0.40 0.15]); yline(ax,0,'k:');
ylabel(ax,'Median unique F+V CV'); grid(ax,'on'); title(ax,'Added behavior beyond original');

axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    xlim(axs(k),[0.5 n+0.5]); xticks(axs(k),x); xticklabels(axs(k),F.session); xtickangle(axs(k),45);
end
title(tl,[P.morphDisplayName ' | fish-level matched four-model comparison']);
saveFigureLocal(fig,outDir,'ALL_FISH_phase_only_and_augmented_fish_level_summary',P);
end

function pairedFishFour(ax,x,phaseOnly,original,augmented,behavior,yLabel)
offset=[-0.27 -0.09 0.09 0.27];
values=[phaseOnly(:),original(:),augmented(:),behavior(:)];
colors=modelColors();
hold(ax,'on');
for i=1:numel(x)
    vals=values(i,:);
    if all(isfinite(vals))
        plot(ax,x(i)+offset,vals,'-','Color',[.70 .70 .70], ...
            'HandleVisibility','off');
    end
end
names={'Phase-only','Original','Augmented','Behavior-only'};
for j=1:4
    plot(ax,x+offset(j),values(:,j),'o','MarkerFaceColor',colors(j,:), ...
        'MarkerEdgeColor','k','DisplayName',names{j});
end
ylabel(ax,yLabel); grid(ax,'on'); title(ax,[yLabel ': matched frames']);
legend(ax,'Location','best');
end

function fourModelPlot(ax,M,yLabel,showMatchedLines)
names={'Phase-only','Original','Augmented','Behavior-only'};
colors=modelColors();
ok=all(isfinite(M),2);
M=M(ok,:);
hold(ax,'on');
if isempty(M)
    xlim(ax,[0.5 4.5]); xticks(ax,1:4); xticklabels(ax,names);
    ylabel(ax,yLabel); grid(ax,'on'); title(ax,'No complete matched fits');
    return;
end
if showMatchedLines
    for i=1:size(M,1)
        plot(ax,1:4,M(i,:),'-','Color',[.78 .78 .78], ...
            'LineWidth',0.4,'HandleVisibility','off');
    end
end
h=gobjects(4,1);
for j=1:4
    if showMatchedLines
        xx=repmat(j,size(M,1),1);
    else
        xx=j+0.14*sin((1:size(M,1))'*sqrt(17+j));
    end
    h(j)=scatter(ax,xx,M(:,j),12,colors(j,:),'filled', ...
        'MarkerFaceAlpha',0.28,'DisplayName',names{j});
    med=median(M(:,j),'omitnan');
    plot(ax,[j-.22 j+.22],[med med],'-','Color',colors(j,:), ...
        'LineWidth',3,'HandleVisibility','off');
end
xlim(ax,[0.5 4.5]); xticks(ax,1:4); xticklabels(ax,names); xtickangle(ax,25);
ylabel(ax,yLabel); grid(ax,'on');
title(ax,sprintf('%s; n = %d matched neurons',yLabel,size(M,1)));
legend(ax,h,'Location','best');
end

function contributionHistogram(ax,x,nBins,color,xLabelText,titleText)
x=x(isfinite(x));
histogram(ax,x,nBins,'FaceColor',color,'EdgeColor','none');
xline(ax,0,'k:','HandleVisibility','off');
xlabel(ax,xLabelText); ylabel(ax,'Neurons'); grid(ax,'on');
title(ax,[titleText ', blocked CV']);
end

function colors = modelColors()
% Phase-only purple; original gray; augmented blue; behavior-only orange.
colors=[0.55 0.30 0.75; 0.55 0.55 0.55; 0.15 0.55 0.90; 0.90 0.40 0.15];
end

function idx = evenlySpacedIndices(n,maxN)
if n<=maxN; idx=(1:n)'; else; idx=unique(round(linspace(1,n,maxN)))'; end
end

function saveFigureLocal(fig,outDir,name,P)
if ~P.saveFigures; return; end
name = [P.figureNamePrefix name];
path=fullfile(outDir,[name '.' lower(P.figureFormat)]);
try
    exportgraphics(fig,path,'Resolution',200);
catch
    saveas(fig,path);
end
svgPath=fullfile(outDir,[name '.svg']);
if ~strcmpi(path,svgPath)
    try
        exportgraphics(fig,svgPath,'ContentType','vector');
    catch
        print(fig,svgPath,'-dsvg');
    end
end
close(fig);
end

function map = blueWhiteRedMap(n)
if nargin<1; n=257; end
x=linspace(0,1,n)';
map=[min(1,2*x),min(1,2*min(x,1-x)+0.02),min(1,2*(1-x))];
map(round((n+1)/2),:)=[1 1 1];
end

function exportClassifiedCandidate(candidatePath,behaviorPath,allCandidateIDs,tunedOrder,labels,fitTable,info,Tuning,P)
[~,name,ext]=fileparts(candidatePath);
ensureDirectory(P.classificationExport.outputDir);
out=fullfile(P.classificationExport.outputDir, ...
    [name P.classificationExport.suffix ext]);
if exist(out,'file')==2
    if P.classificationExport.overwrite; delete(out);
    else; error('Classified candidate already exists: %s',out); end
end
[ok,msg]=copyfile(candidatePath,out); assert(ok,'Candidate copy failed: %s',msg);
assert(numel(tunedOrder)==numel(labels),'Tuned-cell labels do not match tuned indices.');
neuron_class_label=repmat({'NotPhaseTuned'},numel(allCandidateIDs),1); %#ok<NASGU>
neuron_class_label(tunedOrder)=labels(:);
class_masks=struct('CW',strcmp(neuron_class_label,'CW'), ...
    'CCW',strcmp(neuron_class_label,'CCW'), ...
    'Symmetric',strcmp(neuron_class_label,'Symmetric'), ...
    'NotFit',strcmp(neuron_class_label,'NotFit'), ...
    'NotPhaseTuned',strcmp(neuron_class_label,'NotPhaseTuned')); %#ok<NASGU>
class_candidate_indices=struct('CW',find(class_masks.CW), ...
    'CCW',find(class_masks.CCW),'Symmetric',find(class_masks.Symmetric), ...
    'NotFit',find(class_masks.NotFit), ...
    'NotPhaseTuned',find(class_masks.NotPhaseTuned)); %#ok<NASGU>
class_original_roi_ids=struct('CW',allCandidateIDs(class_masks.CW), ...
    'CCW',allCandidateIDs(class_masks.CCW), ...
    'Symmetric',allCandidateIDs(class_masks.Symmetric), ...
    'NotFit',allCandidateIDs(class_masks.NotFit), ...
    'NotPhaseTuned',allCandidateIDs(class_masks.NotPhaseTuned)); %#ok<NASGU>
source_candidate_path=candidatePath; %#ok<NASGU>
source_behavior_path=behaviorPath; %#ok<NASGU>
behavior_augmented_model_fit_table=fitTable; %#ok<NASGU>
phase_tuning_selection=Tuning; %#ok<NASGU>
behavior_augmented_model_definition=struct( ...
    'phaseOnlyFormula','a*cos(phi)+d*sin(phi)+c0', ...
    'originalFormula','a*cos(phi)+d*sin(phi)+b1*abs(AHV)+b2*AHV+c0', ...
    'behaviorOnlyFormula','b1*abs(AHV)+b2*AHV+bF*F+bV*V+c0', ...
    'augmentedFormula', ...
    'a*cos(phi)+d*sin(phi)+b1*abs(AHV)+b2*AHV+bF*F+bV*V+c0', ...
    'formula','a*cos(phi)+d*sin(phi)+b1*abs(AHV)+b2*AHV+bF*F+bV*V+c0', ...
    'derivedAmplitude','b0=hypot(a,d)','derivedPreferredPhase','atan2(d,a)', ...
    'comparisonPolicy','Identical neurons, frames and blocked-CV folds for all four models', ...
    'classInfo',info); %#ok<NASGU>
save(out,'neuron_class_label','class_masks','class_candidate_indices', ...
    'class_original_roi_ids','source_candidate_path','source_behavior_path', ...
    'behavior_augmented_model_fit_table', ...
    'phase_tuning_selection','behavior_augmented_model_definition','-append');
fprintf('Saved augmented classified copy: %s\n',out);
end

function Z = zscoreColumnsFinite(X)
Z=NaN(size(X));
for j=1:size(X,2); Z(:,j)=zscoreLocal(X(:,j)); end
end

function Z = zscoreLocal(X)
mu=mean(X,1,'omitnan'); sd=std(X,0,1,'omitnan'); sd(~isfinite(sd)|sd<=0)=1;
Z=(X-mu)./sd;
end

function r2 = localR2(y,pred)
ok=isfinite(y)&isfinite(pred); y=y(ok); pred=pred(ok);
if numel(y)<2; r2=NaN; return; end
sst=sum((y-mean(y)).^2); if sst<=0; r2=NaN; else; r2=1-sum((y-pred).^2)/sst; end
end

function r2 = adjustedR2Local(r2,n,k)
if ~isfinite(r2)||n<=k; r2=NaN; else; r2=1-(1-r2)*(n-1)/(n-k); end
end

function [aic,bic] = informationCriteriaLocal(sse,n,k)
if ~isfinite(sse)||sse<=0||n<=0; aic=NaN; bic=NaN; return; end
aic=n*log(sse/n)+2*k; bic=n*log(sse/n)+log(n)*k;
end

function r = safeCorr(x,y)
ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
if numel(x)<3 || std(x)==0 || std(y)==0; r=NaN; return; end
q=corrcoef(x,y); r=q(1,2);
end

function v = safeDivide(a,b)
if isfinite(a)&&isfinite(b)&&b~=0; v=a/b; else; v=NaN; end
end

function ang = cleanAngleUnitsToRad(ang)
ang=double(ang(:)); f=ang(isfinite(ang));
if ~isempty(f) && max(abs(f))>2*pi+0.5; ang=deg2rad(ang); end
ang=wrapToPiLocal(ang);
end

function ang = wrapToPiLocal(ang)
ang=mod(ang+pi,2*pi)-pi;
end

function name = getFolderName(path)
[~,name]=fileparts(path);
end
