function pipelineResults = fit_neuron_HD_AHV_three_models_behavior_tuned_cells_only_ROC(cfg)
%FIT_NEURON_HD_AHV_THREE_MODELS_BEHAVIOR_TUNED_CELLS_ONLY_ROC
% Fit Surface, Molino and Pachon in one reproducible pipeline run.
% Pass the struct returned by pipeline_config; outputs are centralized.

% First select network-phase-tuned candidate neurons using the exact
% ORI_V15 STEP 17 information/shuffle logic, then fit three HD/AHV models
% only to the selected neurons.
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
% For every phase-tuned candidate neuron i, the combined model is:
%   y_i(t) = betaTheta_i*cos(sourcePrefRad_i - phi(t)) ...
%          + b1_i*abs(AHV(t)) + b2_i*AHV(t) + c0_i + epsilon_i(t)
%
% sourcePrefRad_i is fixed before selection and behavior/AHV loading from
% the first circular harmonic of the stored tuning_curves_phi_all_cells. It
% is not re-estimated by any regression below.
%
% The script fits three models on identical samples:
%   phase + AHV: cos(sourcePrefRad_i-phi)+|AHV|+AHV
%   phase only:  cos(sourcePrefRad_i-phi)
%   AHV only:    |AHV|+AHV
%
% Legacy a/d/b0 fields remain in the output for downstream compatibility:
% a stores betaTheta, d is zero, b0 is abs(betaTheta), and prefRad is the
% fixed sourcePrefRad.
%
% Blocked cross-validation uses identical contiguous folds and common OOF
% support for all three models, as required by the downstream comparisons.
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

% Causal calcium kernel used for the AHV event impulse train.
P.kernel.tauSeconds = 5.0;
P.ahvBoutTime = 'start';       % 'start', 'center', or 'end'

P.kernel.ahvDurationSeconds = 20;

% Published convention: positive AHV = CCW.
% In p4_build_behavior_features, native positive angles are right/CW.
P.ahvNativeToPaperSign = -1;

% Mei TURN_BIAS: exclude forward/small bouts from the AHV impulse train.
P.ahvMinAbsBoutAngleRad = 0.239;

P.model.version = 'three_models_fixed_source_pref_turn_bias_v2';

% Retain rotation counts as descriptive QC only. Recording inclusion must not
% depend on how many positive or negative rotations a fish performed.
P.lowTurnQC.enabled = false;
P.lowTurnQC.minRotationsEachDirection = 0;

% Blocked time-series cross-validation. All nested models use identical folds.
P.cv.enabled = true;
P.cv.nBlockedFolds = 5;
P.cv.minTestSamplesPerFold = 20;
P.cv.saveOOFTraces = true; % required by the downstream ROC/AUC analysis

% CW/CCW/Symmetric labels from the phase+AHV model. All valid coefficient
% p-values are pooled across every included fish and morph. Bonferroni is
% then applied separately to betaTheta, b1 and b2 before labels are assigned.
P.class.alpha = 0.05;
P.class.bonferroniScope = 'all_included_cells_across_morphs';
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

% Labels must be assigned only after every included cell has been fitted.
% This pass pools coefficient p-values across all morphs, applies a separate
% Bonferroni family to betaTheta, b1 and b2, updates every saved result, and
% only then writes the classified candidate copies.
pipelineResults.pooledBonferroni = applyPooledBonferroniAcrossPipeline( ...
    pipelineResults.morphs, cfg);

save(fullfile(cfg.ModelDataDir, 'three_model_fit_pipeline_summary.mat'), ...
    'pipelineResults', '-v7.3');
disp(' ');
fprintf('All morphs complete. Processed data: %s', cfg.ModelDataDir);
fprintf('%s', newline);
end

function audit = applyPooledBonferroniAcrossPipeline(morphOutputs,cfg)
% Pool raw coefficient p-values across all successfully fitted cells from
% every included fish and morph. Each coefficient has its own family.
pooledPA=[]; pooledPB1=[]; pooledPB2=[]; nIncludedFish=0; pooledAlpha=NaN;
for m=1:numel(morphOutputs)
    if isempty(morphOutputs(m).resultFile) || exist(morphOutputs(m).resultFile,'file')~=2
        continue;
    end
    L=load(morphOutputs(m).resultFile,'Results','P');
    if isnan(pooledAlpha); pooledAlpha=L.P.class.alpha;
    else; assert(L.P.class.alpha==pooledAlpha, ...
            'P.class.alpha must be identical across morphs.'); end
    for s=1:numel(L.Results)
        R=L.Results{s};
        if ~isstruct(R) || ~isfield(R,'fit') || isempty(R.fit) || ...
                ~isfield(R.fit,'nFitSamples')
            continue;
        end
        eligible=R.fit.nFitSamples>=L.P.minSamplesForNeuronFit;
        pooledPA=[pooledPA;R.fit.pA(eligible & isfinite(R.fit.pA))]; %#ok<AGROW>
        pooledPB1=[pooledPB1;R.fit.pB1(eligible & isfinite(R.fit.pB1))]; %#ok<AGROW>
        pooledPB2=[pooledPB2;R.fit.pB2(eligible & isfinite(R.fit.pB2))]; %#ok<AGROW>
        nIncludedFish=nIncludedFish+1;
    end
end

familyN=struct('betaTheta',numel(pooledPA),'b1',numel(pooledPB1), ...
    'b2',numel(pooledPB2));
assert(familyN.b2>0,'No valid b2 p-values were available for pooled classification.');

audit=struct();
audit.scope='all included fitted cells across all fish and morphs';
audit.method='Bonferroni: adjustedP=min(1,Nfamily*rawP)';
audit.alpha=pooledAlpha;
audit.nIncludedFish=nIncludedFish;
audit.nTestsBetaTheta=familyN.betaTheta;
audit.nTestsB1=familyN.b1;
audit.nTestsB2=familyN.b2;
audit.alphaThresholdBetaTheta=safeBonferroniThreshold(audit.alpha,familyN.betaTheta);
audit.alphaThresholdB1=safeBonferroniThreshold(audit.alpha,familyN.b1);
audit.alphaThresholdB2=safeBonferroniThreshold(audit.alpha,familyN.b2);

fprintf(['\nApplying pooled Bonferroni across %d fitted fish: ' ...
    'betaTheta N=%d, b1 N=%d, b2 N=%d.\n'],nIncludedFish, ...
    familyN.betaTheta,familyN.b1,familyN.b2);

for m=1:numel(morphOutputs)
    resultFile=morphOutputs(m).resultFile;
    if isempty(resultFile) || exist(resultFile,'file')~=2; continue; end
    L=load(resultFile,'Results','AllNeurons','FishSummary', ...
        'AllCandidateTuning','TuningFishSummary','P','sessions');
    Results=L.Results; P=L.P;

    for s=1:numel(Results)
        R=Results{s};
        if ~isstruct(R) || ~isfield(R,'fit') || isempty(R.fit) || ...
                ~isfield(R,'fitTable') || isempty(R.fitTable)
            continue;
        end

        [R.fit,labels,ahs,info]=classifyPhaseAHVPooled(R.fit,P,familyN);
        R.classLabel=labels;
        R.classInfo=info;

        T=R.fitTable;
        T.pAAdjusted=R.fit.pAAdjusted;
        T.pBetaThetaAdjusted=R.fit.pAAdjusted;
        T.pB1Adjusted=R.fit.pB1Adjusted;
        T.pB2Adjusted=R.fit.pB2Adjusted;
        T.class=labels(:);
        T.ahsLabel=ahs(:);
        T.pThresholdBetaTheta=repmat(info.alphaThresholdBetaTheta,height(T),1);
        T.pThresholdB1=repmat(info.alphaThresholdB1,height(T),1);
        T.pThresholdAHV=repmat(info.alphaThresholdB2,height(T),1);
        T.bonferroniNTestsBetaTheta=repmat(familyN.betaTheta,height(T),1);
        T.bonferroniNTestsB1=repmat(familyN.b1,height(T),1);
        T.bonferroniNTestsB2=repmat(familyN.b2,height(T),1);
        T.bonferroniNTests=repmat(familyN.b2,height(T),1);
        T.nValidNeuronsThisFish=repmat(info.nValidNeuronsThisFish,height(T),1);
        R.fitTable=T;

        Results{s}=R;
        perFishFile=fullfile(morphOutputs(m).sessionDataDir, ...
            [R.name '_three_models_phase_tuned_only.mat']);
        Result=R; %#ok<NASGU>
        save(perFishFile,'Result','P','-v7.3');

        if P.classificationExport.enabled
            exportClassifiedCandidate(R.paths.candidatePath,R.paths.behaviorPath, ...
                R.tuning.candidateIDs,R.tunedCandidateOrder,labels,T,info,R.tuning,P);
        end
    end

    fitTables=cellfun(@resultFitTableOrEmpty,Results,'UniformOutput',false);
    useT=cellfun(@(x) istable(x) && height(x)>0,fitTables);
    if any(useT); AllNeurons=vertcat(fitTables{useT}); else; AllNeurons=table(); end
    fishTables=cellfun(@resultFishTableOrEmpty,Results,'UniformOutput',false);
    useF=cellfun(@(x) istable(x) && height(x)>0,fishTables);
    if any(useF); FishSummary=vertcat(fishTables{useF}); else; FishSummary=table(); end
    AllCandidateTuning=L.AllCandidateTuning;
    TuningFishSummary=L.TuningFishSummary;
    sessions=L.sessions;
    save(resultFile,'Results','AllNeurons','FishSummary','AllCandidateTuning', ...
        'TuningFishSummary','P','sessions','-v7.3');

    nameStem='HD_AHV_three_models_phase_tuned_only';
    if ~isempty(AllNeurons)
        writeNeuronScalarCSV(AllNeurons,fullfile(cfg.ModelDataDir, ...
            sprintf('%s_neurons_%s.csv',nameStem,P.morphName)));
    end
    if ~isempty(FishSummary)
        writetable(FishSummary,fullfile(cfg.ModelDataDir, ...
            sprintf('%s_fish_summary_%s.csv',nameStem,P.morphName)));
    end
end
end

function [Fit,labels,ahs,info] = classifyPhaseAHVPooled(Fit,P,familyN)
N=numel(Fit.b2);
eligible=Fit.nFitSamples>=P.minSamplesForNeuronFit;
validA=eligible & isfinite(Fit.a) & isfinite(Fit.pA);
validB1=eligible & isfinite(Fit.b1) & isfinite(Fit.pB1);
validB2=eligible & isfinite(Fit.b2) & isfinite(Fit.pB2);

Fit.pAAdjusted=nan(N,1);
Fit.pB1Adjusted=nan(N,1);
Fit.pB2Adjusted=nan(N,1);
if familyN.betaTheta>0
    Fit.pAAdjusted(validA)=min(1,familyN.betaTheta.*Fit.pA(validA));
end
if familyN.b1>0
    Fit.pB1Adjusted(validB1)=min(1,familyN.b1.*Fit.pB1(validB1));
end
Fit.pB2Adjusted(validB2)=min(1,familyN.b2.*Fit.pB2(validB2));

labels=repmat({'NotFit'},N,1);
[sharedLabels,sharedAdjusted]=hdahv.classifyB2(Fit.b2,Fit.pB2, ...
    familyN.b2,P.class.alpha,P.ahvPositiveIsCCW,validB2);
labels=cellstr(sharedLabels); Fit.pB2Adjusted=sharedAdjusted;

ahs=repmat({'AHS_not_fit'},N,1);
ahs(validB1)={'AHS_not_significant'};
sigB1=validB1 & Fit.pB1Adjusted<P.class.alpha;
ahs(sigB1 & Fit.b1>0)={'AHS_positive'};
ahs(sigB1 & Fit.b1<0)={'AHS_negative'};

info=struct('alpha',P.class.alpha, ...
    'alphaBonferroni',safeBonferroniThreshold(P.class.alpha,familyN.b2), ...
    'alphaThresholdBetaTheta',safeBonferroniThreshold(P.class.alpha,familyN.betaTheta), ...
    'alphaThresholdB1',safeBonferroniThreshold(P.class.alpha,familyN.b1), ...
    'alphaThresholdB2',safeBonferroniThreshold(P.class.alpha,familyN.b2), ...
    'bonferroniMode','pooled_across_all_included_cells_and_morphs', ...
    'nBonferroniTests',familyN.b2, ...
    'nBonferroniTestsBetaTheta',familyN.betaTheta, ...
    'nBonferroniTestsB1',familyN.b1, ...
    'nBonferroniTestsB2',familyN.b2, ...
    'nValidNeuronsThisFish',sum(validB2), ...
    'adjustedPFormula','min(1,Nfamily*rawP)', ...
    'classificationRule','Sign of b2 when pooled-Bonferroni pB2Adjusted < alpha', ...
    'pValueCaveat',['Conventional OLS p-values do not correct for temporal ' ...
    'autocorrelation; use continuous coefficients and blocked CV for primary inference.']);
end

function value = safeBonferroniThreshold(alpha,nTests)
if nTests>0; value=alpha/nTests; else; value=NaN; end
end

function T = resultFitTableOrEmpty(R)
if isstruct(R) && isfield(R,'fitTable') && istable(R.fitTable); T=R.fitTable;
else; T=table(); end
end

function T = resultFishTableOrEmpty(R)
if isstruct(R) && isfield(R,'fishSummary') && istable(R.fishSummary); T=R.fishSummary;
else; T=table(); end
end

function writeNeuronScalarCSV(T,path)
traceVariables=intersect({'cvObserved','cvPredPhaseOnly', ...
    'cvPredBehavior','cvPredFull','cvPredOriginal', ...
    'cvPredPhaseAHV','cvPredAHVOnly'},T.Properties.VariableNames,'stable');
if ~isempty(traceVariables); T=removevars(T,traceVariables); end
writetable(T,path);
end

function morphOutput = runOneMorph(P, modelDataDir, sessionDataDir, figureDir)
%% ======================== DISCOVER / RUN ========================
[sessions, P.rootDataDir] = discoverSessions(P.rootDataDir, P);

nameStem = "HD_AHV_three_models_phase_tuned_only";
resultFile = fullfile(modelDataDir, sprintf("%s_results_%s.mat",nameStem,P.morphName));
cachedResults={}; cachedNames=strings(0,1);
if isfile(resultFile)
    Q=load(resultFile,"Results","P");
    compatible=isfield(Q,"Results")&&iscell(Q.Results)&&isfield(Q,"P")&& ...
        isfield(Q.P,"model")&&isfield(Q.P.model,"version")&&strcmp(Q.P.model.version,P.model.version)&& ...
        isfield(Q.P.model,"zscoreActivityWithinFitWindow")&& ...
        logical(Q.P.model.zscoreActivityWithinFitWindow)==logical(P.model.zscoreActivityWithinFitWindow);
    if compatible
        cachedResults=Q.Results; cachedNames=strings(numel(cachedResults),1);
        for i=1:numel(cachedResults)
            if isstruct(cachedResults{i})&&isfield(cachedResults{i},"name"), cachedNames(i)=string(cachedResults{i}.name); end
        end
    end
end

Results = cell(numel(sessions),1);
allTables = cell(numel(sessions),1);
fishRows = cell(numel(sessions),1);
allTuningTables = cell(numel(sessions),1);
tuningFishRows = cell(numel(sessions),1);

for s = 1:numel(sessions)
    fprintf('\n================ %s (%d/%d) ================\n', ...
        sessions(s).name, s, numel(sessions));
    cachedIndex=find(cachedNames==string(sessions(s).name),1);
    if ~isempty(cachedIndex)&&isReusableSessionResult(cachedResults{cachedIndex})
        Results{s}=cachedResults{cachedIndex};
        allTables{s}=Results{s}.fitTable; fishRows{s}=Results{s}.fishSummary;
        allTuningTables{s}=Results{s}.tuningTable; tuningFishRows{s}=Results{s}.tuningSummary;
        fprintf("Reusing successful cached session %s.\n",sessions(s).name);
        continue
    end
    fprintf("Rerunning failed or missing session %s.\n",sessions(s).name);
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

save(resultFile, ...
    'Results','AllNeurons','FishSummary','AllCandidateTuning', ...
    'TuningFishSummary','P','sessions','-v7.3');
if ~isempty(AllNeurons)
    % Variable-length OOF traces are retained in the MAT file but omitted
    % from the flat CSV, which is intended for scalar neuron summaries.
    AllNeuronsCSV = AllNeurons;
    traceVariables = intersect({'cvObserved','cvPredPhaseOnly', ...
        'cvPredBehavior','cvPredFull','cvPredOriginal', ...
        'cvPredPhaseAHV','cvPredAHVOnly'}, ...
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

% The legacy comparison plots require reduced-model fits and are therefore
% intentionally not generated by this single-model pipeline.

morphOutput = struct('morph', P.morphDisplayName, ...
    'resultFile', resultFile, 'sessionDataDir', sessionDataDir, ...
    'figureDir', figureDir, 'nSessions', numel(sessions), ...
    'nFittedNeurons', height(AllNeurons));
fprintf('Morph %s complete: %s', P.morphDisplayName, resultFile);
fprintf('%s', newline);
end

%% ======================== LOCAL FUNCTIONS ========================

function tf=isReusableSessionResult(R)
tf=isstruct(R)&&~isfield(R,"error")&&isfield(R,"name")&& ...
    isfield(R,"fitTable")&&istable(R.fitTable)&&isfield(R,"fishSummary")&&istable(R.fishSummary)&& ...
    isfield(R,"tuningTable")&&istable(R.tuningTable)&&isfield(R,"tuningSummary")&&istable(R.tuningSummary);
end

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
d = d(startsWith(lower(string({d.name})), "rec"));
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
[Y,tCa,phase,phaseOk,sourcePrefRad,candidateIDs,Harmonic,loadDiag] = loadCandidate(C,P);

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
        'tuningHarmonicDiagnostics',Harmonic, ...
        'tuningTable',tuningTable,'tuningSummary',tuningSummary, ...
        'fitTable',table(),'fishSummary',table());
    save(fullfile(dataOutDir,[sess.name '_three_models_phase_tuned_only.mat']), ...
        'Result','P','-v7.3');
    return;
end

Y = Y(:,tunedCandidateOrder);
sourcePrefRad = sourcePrefRad(tunedCandidateOrder);
Harmonic = subsetHarmonicDiagnostics(Harmonic,tunedCandidateOrder);
candidateIDs = candidateIDs(tunedCandidateOrder);
assert(numel(sourcePrefRad)==numel(candidateIDs), ...
    'Preferred-phase vector lost alignment with selected candidate IDs.');
missingSourcePref = ~isfinite(sourcePrefRad);
if any(missingSourcePref)
    error(['First-circular-harmonic sourcePrefRad is missing or invalid for ' ...
        'selected candidate ID(s): %s. No fallback preferred phase is allowed.'], ...
        strjoin(compose('%.15g',candidateIDs(missingSourcePref)),', '));
end
assert(all(sourcePrefRad > -pi & sourcePrefRad <= pi), ...
    'First-harmonic sourcePrefRad must always be in (-pi, pi].');

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
        'tuningHarmonicDiagnostics',Harmonic, ...
        'tuning',Tuning,'tuningTable',tuningTable, ...
        'tuningSummary',tuningSummary,'fitTable',table(),'fishSummary',table());
    save(fullfile(dataOutDir,[sess.name '_three_models_phase_tuned_only.mat']), ...
        'Result','P','-v7.3');
    return;
end

common = isfinite(tCa) & isfinite(phase) & isfinite(B.ahv);
if P.useTimeOverlapOnly
    Y = Y(common,:); tUsed = tCa(common); phase = phase(common);
    ahv = B.ahv(common);
else
    tUsed = tCa; ahv = B.ahv;
end
assert(sum(all(isfinite([phase,ahv]),2)) >= P.minSamplesForNeuronFit, ...
    'Too few common finite samples after behavior alignment.');

if P.model.zscoreActivityWithinFitWindow; Y = hdahv.preprocessActivity(Y,true); end

% Fish-level diagnostics use the same fixed-preference phase regressor as
% the neuron fits. Because preferred direction differs across cells, pool
% all neuron-by-time design rows for these summary diagnostics.
phaseFixed = cos(sourcePrefRad(:)' - phase(:));
predictors = [phaseFixed(:),repmat(abs(ahv(:)),size(Y,2),1), ...
    repmat(ahv(:),size(Y,2),1)];
predictorNames = {'cosPreferredMinusPhi','absAHV','AHV'};
[predictorCorr,predictorVIF,conditionNumberStandardized] = ...
    predictorDiagnostics(predictors);

fprintf('Fitting phase + AHV, AHV-only and phase-only models: %d neurons, %d frames.\n', ...
    size(Y,2),size(Y,1));
Fit = fitAllNeurons(Y,phase,ahv,sourcePrefRad,P);
% Classification is deliberately deferred until every fish and morph has
% been fitted, so coefficient p-values can be corrected in pooled families.
[classLabel,ahsLabel,classInfo] = makePendingClassification(Fit,P);

N = size(Y,2);
session = repmat({sess.name},N,1);
fitTable = makeFitTable(session,tunedCandidateOrder,candidateIDs,Fit, ...
    classLabel,ahsLabel,classInfo,predictorVIF,conditionNumberStandardized);
fitTable.phaseInfoBits = Tuning.infoBits(tunedCandidateOrder);
fitTable.phaseInfoShuffleP = Tuning.pShuffle(tunedCandidateOrder);
fitTable.phaseInfoFDRQ = Tuning.qFDR(tunedCandidateOrder);
fitTable.phaseTuningPrefRad = Tuning.prefRad(tunedCandidateOrder);
fitTable.phaseTuningPrefDeg = rad2deg(Tuning.prefRad(tunedCandidateOrder));
fitTable.sourcePrefRad = sourcePrefRad;
fitTable.sourcePrefDeg = rad2deg(sourcePrefRad);
fitTable.sourcePrefMethod = repmat( ...
    {'first circular harmonic of stored tuning_curves_phi_all_cells'},N,1);
fitTable.tuningHarmonicIntercept = Harmonic.intercept;
fitTable.tuningHarmonicCosCoeff = Harmonic.cosCoeff;
fitTable.tuningHarmonicSinCoeff = Harmonic.sinCoeff;
fitTable.tuningHarmonicAmplitude = Harmonic.amplitude;
fitTable.tuningHarmonicR2 = Harmonic.r2;
fitTable.harmonicVsInternalPrefDifferenceRad = atan2( ...
    sin(sourcePrefRad-fitTable.phaseTuningPrefRad), ...
    cos(sourcePrefRad-fitTable.phaseTuningPrefRad));
fitTable.harmonicVsInternalPrefDifferenceDeg = ...
    rad2deg(fitTable.harmonicVsInternalPrefDifferenceRad);
fitTable.phaseTuningVectorStrength = Tuning.vectorStrength(tunedCandidateOrder);
fitTable.selectedByPhaseTuning = true(N,1);

fishSummary = makeFishSummary(sess.name,fitTable,behaviorDiag, ...
    predictorVIF,conditionNumberStandardized);
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
Result.tuningHarmonicDiagnostics = Harmonic;
Result.fit = Fit;
Result.classLabel = classLabel;
Result.classInfo = classInfo;
Result.fitTable = fitTable;
Result.fishSummary = fishSummary;

perFishFile = fullfile(dataOutDir,[sess.name '_three_models_phase_tuned_only.mat']);
save(perFishFile,'Result','P','-v7.3');
plotRegressorQC(Result,P,figureOutDir);
plotPredictorCorrelation(Result,P,figureOutDir);

end

function [Y,tCa,phase,phaseOk,sourcePrefRad,candidateIDs,H,D] = loadCandidate(C,P)
required = {'calcium_traces','time_s','network_phase_rad','candidate_cell_ids', ...
    'phi_bin_centers_rad','tuning_curves_phi_all_cells'};
for k = 1:numel(required)
    assert(isfield(C,required{k}),'Candidate file missing %s.',required{k});
end
candidateIDs = double(C.candidate_cell_ids(:));
assert(numel(unique(candidateIDs))==numel(candidateIDs), ...
    'candidate_cell_ids must be unique within each candidate file.');
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
[sourcePrefRad,H] = getTuningHarmonicPreferredPhase(C,candidateIDs);
assert(numel(sourcePrefRad)==numel(candidateIDs) && ...
    isequaln(H.candidateIDs,candidateIDs), ...
    'First-harmonic preferred phases are not aligned with candidate_cell_ids.');
D = struct('activitySource','calcium_traces','phaseSource','network_phase_rad', ...
    'phaseOkSource',phaseOkSource, ...
    'referencePreferredPhaseSource', ...
    'first circular harmonic of tuning_curves_phi_all_cells', ...
    'preferredPhaseMethod', ...
    'atan2(sine coefficient, cosine coefficient) from stored tuning curve', ...
    'nFrames',size(Y,1), ...
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

function [sourcePrefRad,H] = getTuningHarmonicPreferredPhase(C,candidateIDs)
% Derive the sole model reference phase from each stored tuning curve.
phiBins = double(C.phi_bin_centers_rad(:));
Q = double(C.tuning_curves_phi_all_cells);
candidateIDs = double(candidateIDs(:));
nNeurons = numel(candidateIDs);

% Convert to: rows = phase bins, columns = neurons.
if size(Q,1) == numel(phiBins) && size(Q,2) == nNeurons
    % Already correctly oriented.
elseif size(Q,2) == numel(phiBins) && size(Q,1) == nNeurons
    Q = Q';
else
    error(['tuning_curves_phi_all_cells dimensions [%d %d] are incompatible ' ...
        'with %d phase bins and %d candidate neurons.'], ...
        size(Q,1),size(Q,2),numel(phiBins),nNeurons);
end

X = [ones(numel(phiBins),1),cos(phiBins),sin(phiBins)];
sourcePrefRad = nan(nNeurons,1);
tuningHarmonicIntercept = nan(nNeurons,1);
tuningHarmonicCosCoeff = nan(nNeurons,1);
tuningHarmonicSinCoeff = nan(nNeurons,1);
tuningAmplitude = nan(nNeurons,1);
tuningHarmonicR2 = nan(nNeurons,1);

for i = 1:nNeurons
    q = Q(:,i);
    valid = isfinite(phiBins) & isfinite(q);
    if nnz(valid) < 3
        continue;
    end
    beta = X(valid,:) \ q(valid);
    tuningHarmonicIntercept(i) = beta(1);
    tuningHarmonicCosCoeff(i) = beta(2);
    tuningHarmonicSinCoeff(i) = beta(3);
    sourcePrefRad(i) = atan2(beta(3),beta(2));
    sourcePrefRad(i) = atan2(sin(sourcePrefRad(i)),cos(sourcePrefRad(i)));
    % Enforce the requested half-open circular interval (-pi, pi].
    if sourcePrefRad(i) <= -pi
        sourcePrefRad(i) = pi;
    end
    tuningAmplitude(i) = hypot(beta(2),beta(3));
    qhat = X(valid,:) * beta;
    ssRes = sum((q(valid)-qhat).^2);
    ssTot = sum((q(valid)-mean(q(valid))).^2);
    if isfinite(ssTot) && ssTot > 0
        tuningHarmonicR2(i) = 1-ssRes/ssTot;
    end
end

finitePref = isfinite(sourcePrefRad);
assert(all(sourcePrefRad(finitePref)>-pi & sourcePrefRad(finitePref)<=pi), ...
    'First-harmonic sourcePrefRad escaped the required (-pi, pi] interval.');
H = struct('candidateIDs',candidateIDs, ...
    'intercept',tuningHarmonicIntercept, ...
    'cosCoeff',tuningHarmonicCosCoeff, ...
    'sinCoeff',tuningHarmonicSinCoeff, ...
    'amplitude',tuningAmplitude,'r2',tuningHarmonicR2, ...
    'method','first circular harmonic of stored tuning_curves_phi_all_cells');
end

function H = subsetHarmonicDiagnostics(H,tunedCandidateOrder)
fields = {'candidateIDs','intercept','cosCoeff','sinCoeff','amplitude','r2'};
for k = 1:numel(fields)
    assert(numel(H.(fields{k})) >= max(tunedCandidateOrder), ...
        'Harmonic diagnostic %s is not aligned with candidate order.',fields{k});
    H.(fields{k}) = H.(fields{k})(tunedCandidateOrder);
end
end

function [B,D] = buildBehaviorRegressors(beh,tCa,P)
% Construct only the Mei-style AHV predictor at calcium resolution.
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


% Convert repository convention to Mei convention: positive = CCW.
dthetaPaper = P.ahvNativeToPaperSign .* dtheta;

% Match Mei's TURN_BIAS definition: forward/small bouts do not enter AHV.
validAHVBout = validBout & ...
    abs(dthetaPaper) > P.ahvMinAbsBoutAngleRad;

boutRows = find(validAHVBout);
boutTimes = tBeh(boutFrame(boutRows));

% Map bout starts to their nearest imaging frames.
assert(all(isfinite(tCa)) && all(diff(tCa)>0), ...
    'Calcium timestamps must be finite and strictly increasing.');

imagingIndex = round(interp1( ...
    tCa(:),(1:numel(tCa))',boutTimes,'linear',NaN));

keep = isfinite(imagingIndex) & ...
    imagingIndex >= 1 & imagingIndex <= numel(tCa);

angleImpulseCa = accumarray( ...
    imagingIndex(keep), ...
    dthetaPaper(boutRows(keep)), ...
    [numel(tCa),1],@sum,0);

% Three-second exponential kernel normalized to unit sum.
dtCa = median(diff(tCa));
kernelTime = (0:floor(P.kernel.ahvDurationSeconds/dtCa))' .* dtCa;

kAHV = exp(-kernelTime/P.kernel.tauSeconds);
kAHV = kAHV ./ sum(kAHV);

assert(abs(sum(kAHV)-1) < 1e-12, ...
    'Mei AHV kernel must sum to one.');

ahvCa = causalConvolutionFFT(angleImpulseCa,kAHV);

finiteBehTime = tBeh(isfinite(tBeh));
inBehaviorSupport = tCa >= min(finiteBehTime) & ...
    tCa <= max(finiteBehTime);
ahvCa(~inBehaviorSupport) = NaN;

B = struct('ahv',ahvCa);

D = struct();
D.fpsBehavior=fps; D.nBehaviorFrames=nBeh; D.nBouts=nB;
D.nValidBouts=sum(validBout); D.nAHVBouts=sum(validAHVBout);
D.dthetaSource=dthetaSource; D.behaviorTimeSource=timeSource;
D.kernelTauSeconds=P.kernel.tauSeconds;
D.kernelDurationSeconds=P.kernel.ahvDurationSeconds;
D.turnBiasRad=P.ahvMinAbsBoutAngleRad;
D.totalPositiveRotations=sum(max(dtheta(validBout),0),'omitnan')/(2*pi);
D.totalNegativeRotations=sum(max(-dtheta(validBout),0),'omitnan')/(2*pi);
D.minDirectionalRotations=min(D.totalPositiveRotations,D.totalNegativeRotations);
D.ahvRangeCalcium = [ ...
    min(ahvCa,[],'omitnan'), ...
    max(ahvCa,[],'omitnan')];
D.nCalciumSamplesAHV=sum(isfinite(B.ahv));
D.corrAbsAHV_AHV=safeCorr(abs(B.ahv),B.ahv);
end

function n = inferBehaviorLength(beh)
n = NaN;
priority={'behaviorTimeSec','heading_est','tail_angle'};
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

function Fit = fitAllNeurons(Y,phase,ahv,sourcePrefRad,P)
% Full predictor order: cos(sourcePrefRad-phi), absAHV, AHV, intercept.
[~,N]=size(Y); nanN=NaN(N,1);
fields={'a','d','b0','b1','b2','c0', ...
    'seA','seD','seB0','sePrefRad','sePrefDeg','seB1','seB2','seC0', ...
    'pA','pD','pB1','pB2','pC0','fHDJoint','pHDJoint','partialR2HD', ...
    'prefRad','prefDeg','prefShiftRad','prefShiftDeg', ...
    'aStd','dStd','b0Std','b1Std','b2Std', ...
    'epsilonStd','epsilonRMS','r2','adjR2','aic','bic', ...
    'phaseOnlyA','phaseOnlyD','phaseOnlyB0','phaseOnlyC0', ...
    'phaseOnlySeA','phaseOnlySeD','phaseOnlySeB0','phaseOnlySePrefRad','phaseOnlySePrefDeg', ...
    'phaseOnlyPA','phaseOnlyPD','phaseOnlyFJoint','phaseOnlyPJoint','phaseOnlyPartialR2', ...
    'phaseOnlyPrefRad','phaseOnlyPrefDeg', ...
    'phaseOnlyAStd','phaseOnlyDStd','phaseOnlyB0Std', ...
    'phaseOnlyR2','phaseOnlyAdjR2','phaseOnlyAIC','phaseOnlyBIC', ...
    'behaviorB1','behaviorB2','behaviorC0', ...
    'behaviorSeB1','behaviorSeB2','behaviorSeC0', ...
    'behaviorPB1','behaviorPB2','behaviorPC0', ...
    'behaviorB1Std','behaviorB2Std', ...
    'behaviorR2','behaviorAdjR2','behaviorAIC','behaviorBIC', ...
    'behaviorEpsilonStd','behaviorEpsilonRMS', ...
    'fAHVBeyondPhase','pAHVBeyondPhase','partialR2AHVBeyondPhase', ...
    'cvR2PhaseOnly','cvR2Behavior','cvR2Full', ...
    'cvUniquePhase','cvDeltaR2Phase','cvUniqueAHV','cvDeltaR2AHV', ...
    'cvSSEPhaseOnly','cvSSEBehavior','cvSSEFull', ...
    'designRank','conditionNumberRaw','conditionNumberStandardized', ...
    'vifCosPhi','vifSinPhi','vifAbsAHV','vifAHV','maxVIF'};
Fit=struct();
for k=1:numel(fields); Fit.(fields{k})=nanN; end
Fit.cvObserved=cell(N,1);
Fit.cvPredPhaseOnly=cell(N,1);
Fit.cvPredBehavior=cell(N,1);
Fit.cvPredFull=cell(N,1);
Fit.nFitSamples=zeros(N,1);
Fit.responseScale = ternaryActivityScale(P.model.zscoreActivityWithinFitWindow);
Fit.modelFormula=[Fit.responseScale ' = betaTheta*cos(sourcePrefRad-phi) + ' ...
    'b1*abs(AHV) + b2*AHV + c0 + epsilon; sourcePrefRad from first ' ...
    'circular harmonic of stored tuning_curves_phi_all_cells'];
Fit.preferredPhaseSource = ...
    'first circular harmonic of stored tuning_curves_phi_all_cells';
Fit.coefficientOrder={'betaTheta','b1','b2','c0'};
Fit.phaseOnlyModelFormula=[Fit.responseScale ' = betaTheta*cos(sourcePrefRad-phi) + ' ...
    'c0 + epsilon; sourcePrefRad from first circular harmonic of stored ' ...
    'tuning_curves_phi_all_cells'];
Fit.phaseOnlyCoefficientOrder={'phaseOnlyBetaTheta','phaseOnlyC0'};
Fit.behaviorModelFormula=[Fit.responseScale ' = b1*abs(AHV) + b2*AHV + c0 + epsilon'];
Fit.behaviorCoefficientOrder={'behaviorB1','behaviorB2','behaviorC0'};

phase=phase(:); ahv=ahv(:); sourcePrefRad=sourcePrefRad(:);
sharedFull=hdahv.fitPhaseAHV(Y,phase,ahv,sourcePrefRad,P.minSamplesForNeuronFit);

for i=1:N
    y=Y(:,i);
    if ~isfinite(sourcePrefRad(i)); continue; end
    phaseRegressor=cos(sourcePrefRad(i)-phase);
    Xall=[phaseRegressor,abs(ahv),ahv,ones(numel(y),1)];
    ok=all(isfinite(Xall),2) & isfinite(y);
    if sum(ok)<P.minSamplesForNeuronFit; continue; end
    Xfull=Xall(ok,:); yy=y(ok); n=size(Xfull,1);
    Xphase=Xfull(:,[1 4]);
    Xbehavior=Xfull(:,2:4);
    kFull=size(Xfull,2); dofFull=n-kFull;
    if dofFull<=0; continue; end
    yz=zscoreLocal(yy);

    % Model 1: phase + AHV.
    beta=[sharedFull.betaTheta(i);sharedFull.b1(i);sharedFull.b2(i);sharedFull.c0(i)];
    assert(all(isfinite(beta)),'Shared phase+AHV fit unexpectedly invalid.'); pred=Xfull*beta; resid=yy-pred;
    sseFull=sum(resid.^2); mse=sseFull/dofFull;
    se=[sharedFull.seBetaTheta(i);sharedFull.seB1(i);sharedFull.seB2(i);sharedFull.seC0(i)];
    p=[sharedFull.pBetaTheta(i);sharedFull.pB1(i);sharedFull.pB2(i);sharedFull.pC0(i)];
    a=beta(1); d=0; b0=abs(a); pref=sourcePrefRad(i);
    Fit.a(i)=a; Fit.d(i)=d; Fit.b0(i)=b0;
    Fit.b1(i)=beta(2); Fit.b2(i)=beta(3); Fit.c0(i)=beta(4);
    Fit.seA(i)=se(1); Fit.seD(i)=0; Fit.seB1(i)=se(2);
    Fit.seB2(i)=se(3); Fit.seC0(i)=se(4);
    Fit.pA(i)=p(1); Fit.pD(i)=NaN; Fit.pB1(i)=p(2);
    Fit.pB2(i)=p(3); Fit.pC0(i)=p(4);
    Fit.prefRad(i)=wrapToPiLocal(pref); Fit.prefDeg(i)=rad2deg(Fit.prefRad(i));
    Fit.prefShiftRad(i)=0; Fit.prefShiftDeg(i)=0;
    Fit.seB0(i)=se(1); Fit.sePrefRad(i)=0; Fit.sePrefDeg(i)=0;
    Z=zscoreLocal(Xfull(:,1:3));
    betaZ=solveOLS([Z,ones(n,1)],yz);
    Fit.aStd(i)=betaZ(1); Fit.dStd(i)=0;
    Fit.b0Std(i)=abs(betaZ(1));
    Fit.b1Std(i)=betaZ(2); Fit.b2Std(i)=betaZ(3);
    Fit.epsilonStd(i)=std(resid,0,'omitnan');
    Fit.epsilonRMS(i)=sqrt(mean(resid.^2,'omitnan'));
    Fit.r2(i)=localR2(yy,pred);
    Fit.adjR2(i)=adjustedR2Local(Fit.r2(i),n,kFull);
    [Fit.aic(i),Fit.bic(i)]=informationCriteriaLocal(sseFull,n,kFull);

    % Model 2: phase only, fitted on exactly the same samples.
    betaP=solveOLS(Xphase,yy); predP=Xphase*betaP;
    residP=yy-predP; ssePhase=sum(residP.^2);
    kPhase=size(Xphase,2); dofPhase=n-kPhase;
    Fit.phaseOnlyA(i)=betaP(1); Fit.phaseOnlyD(i)=0;
    Fit.phaseOnlyB0(i)=abs(betaP(1));
    Fit.phaseOnlyC0(i)=betaP(2);
    Fit.phaseOnlyPrefRad(i)=wrapToPiLocal(sourcePrefRad(i));
    Fit.phaseOnlyPrefDeg(i)=rad2deg(Fit.phaseOnlyPrefRad(i));
    Fit.phaseOnlyR2(i)=localR2(yy,predP);
    Fit.phaseOnlyAdjR2(i)=adjustedR2Local(Fit.phaseOnlyR2(i),n,kPhase);
    [Fit.phaseOnlyAIC(i),Fit.phaseOnlyBIC(i)]= ...
        informationCriteriaLocal(ssePhase,n,kPhase);
    if dofPhase>0
        mseP=ssePhase/dofPhase;
        covBetaP=mseP*pinv(Xphase'*Xphase);
        seP=sqrt(max(0,diag(covBetaP)));
        tP=betaP./seP; pP=2*(1-tcdf(abs(tP),dofPhase));
        Fit.phaseOnlySeA(i)=seP(1); Fit.phaseOnlySeD(i)=0;
        Fit.phaseOnlyPA(i)=pP(1); Fit.phaseOnlyPD(i)=NaN;
        Fit.phaseOnlySeB0(i)=seP(1);
        Fit.phaseOnlySePrefRad(i)=0; Fit.phaseOnlySePrefDeg(i)=0;
        [Fit.phaseOnlyFJoint(i),Fit.phaseOnlyPJoint(i),Fit.phaseOnlyPartialR2(i)] = ...
            nestedFTest(yy,Xphase,ssePhase,dofPhase,ones(n,1));
    end
    ZP=zscoreLocal(Xphase(:,1));
    betaPZ=solveOLS([ZP,ones(n,1)],yz);
    Fit.phaseOnlyAStd(i)=betaPZ(1); Fit.phaseOnlyDStd(i)=0;
    Fit.phaseOnlyB0Std(i)=abs(betaPZ(1));

    % Model 3: AHV only, fitted on exactly the same samples.
    betaB=solveOLS(Xbehavior,yy); predB=Xbehavior*betaB;
    residB=yy-predB; sseBehavior=sum(residB.^2);
    kBehavior=size(Xbehavior,2); dofBehavior=n-kBehavior;
    Fit.behaviorB1(i)=betaB(1); Fit.behaviorB2(i)=betaB(2);
    Fit.behaviorC0(i)=betaB(3);
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
        tB=betaB./seB; pB=2*(1-tcdf(abs(tB),dofBehavior));
        Fit.behaviorSeB1(i)=seB(1); Fit.behaviorSeB2(i)=seB(2);
        Fit.behaviorSeC0(i)=seB(3);
        Fit.behaviorPB1(i)=pB(1); Fit.behaviorPB2(i)=pB(2);
        Fit.behaviorPC0(i)=pB(3);
    end
    ZB=zscoreLocal(Xbehavior(:,1:2));
    betaBZ=solveOLS([ZB,ones(n,1)],yz);
    Fit.behaviorB1Std(i)=betaBZ(1); Fit.behaviorB2Std(i)=betaBZ(2);

    % Nested in-sample contributions of phase and AHV to the combined model.
    [Fit.fHDJoint(i),Fit.pHDJoint(i),Fit.partialR2HD(i)] = ...
        nestedFTest(yy,Xfull,sseFull,dofFull,Xbehavior);
    [Fit.fAHVBeyondPhase(i),Fit.pAHVBeyondPhase(i),Fit.partialR2AHVBeyondPhase(i)] = ...
        nestedFTest(yy,Xfull,sseFull,dofFull,Xphase);

    Fit.designRank(i)=rank(Xfull); Fit.conditionNumberRaw(i)=cond(Xfull);
    [~,vif,conditionZ]=predictorDiagnostics(Xfull(:,1:3));
    Fit.conditionNumberStandardized(i)=conditionZ;
    Fit.vifCosPhi(i)=vif(1); Fit.vifSinPhi(i)=NaN;
    Fit.vifAbsAHV(i)=vif(2); Fit.vifAHV(i)=vif(3);
    Fit.maxVIF(i)=max(vif,[],'omitnan'); Fit.nFitSamples(i)=n;

    if P.cv.enabled
        CV=blockedCrossValidation(yy,Xphase,Xbehavior,Xfull,P);
        Fit.cvR2PhaseOnly(i)=CV.r2PhaseOnly;
        Fit.cvR2Behavior(i)=CV.r2Behavior;
        Fit.cvR2Full(i)=CV.r2Full;
        Fit.cvSSEPhaseOnly(i)=CV.ssePhaseOnly;
        Fit.cvSSEBehavior(i)=CV.sseBehavior;
        Fit.cvSSEFull(i)=CV.sseFull;
        Fit.cvUniquePhase(i)=safeOneMinusRatio(CV.sseFull,CV.sseBehavior);
        Fit.cvDeltaR2Phase(i)=CV.r2Full-CV.r2Behavior;
        Fit.cvUniqueAHV(i)=safeOneMinusRatio(CV.sseFull,CV.ssePhaseOnly);
        Fit.cvDeltaR2AHV(i)=CV.r2Full-CV.r2PhaseOnly;
        if P.cv.saveOOFTraces
            Fit.cvObserved{i}=single(CV.observed);
            Fit.cvPredPhaseOnly{i}=single(CV.predPhaseOnly);
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

function CV = blockedCrossValidation(y,Xphase,Xbehavior,Xfull,P)
n=numel(y); nFolds=min(P.cv.nBlockedFolds,floor(n/P.cv.minTestSamplesPerFold));
CV=struct('r2PhaseOnly',NaN,'r2Behavior',NaN,'r2Full',NaN, ...
    'ssePhaseOnly',NaN,'sseBehavior',NaN,'sseFull',NaN, ...
    'observed',[],'predPhaseOnly',[],'predBehavior',[],'predFull',[]);
if nFolds<2; return; end
edges=round(linspace(0,n,nFolds+1));
pP=NaN(n,1); pB=NaN(n,1); pF=NaN(n,1);
for f=1:nFolds
    test=(edges(f)+1):edges(f+1); train=true(n,1); train(test)=false;
    if numel(test)<P.cv.minTestSamplesPerFold || sum(train)<=size(Xfull,2); continue; end
    pP(test)=Xphase(test,:)*solveOLS(Xphase(train,:),y(train));
    pB(test)=Xbehavior(test,:)*solveOLS(Xbehavior(train,:),y(train));
    pF(test)=Xfull(test,:)*solveOLS(Xfull(train,:),y(train));
end
valid=isfinite(y) & isfinite(pP) & isfinite(pB) & isfinite(pF);
if sum(valid)<P.minSamplesForNeuronFit; return; end
yy=y(valid); pP=pP(valid); pB=pB(valid); pF=pF(valid);
% Preserve only the common, aligned OOF support. These are the same samples
% used below for all three CV-R2 values and therefore permit paired ROC/AUC.
CV.observed=yy;
CV.predPhaseOnly=pP;
CV.predBehavior=pB;
CV.predFull=pF;
CV.ssePhaseOnly=sum((yy-pP).^2);
CV.sseBehavior=sum((yy-pB).^2); CV.sseFull=sum((yy-pF).^2);
sst=sum((yy-mean(yy)).^2);
if sst>0
    CV.r2PhaseOnly=1-CV.ssePhaseOnly/sst;
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

function [labels,ahs,info] = makePendingClassification(Fit,P)
N=numel(Fit.b2);
labels=repmat({'PendingPooledBonferroni'},N,1);
ahs=repmat({'AHS_pending_pooled_Bonferroni'},N,1);
valid=isfinite(Fit.b2) & isfinite(Fit.pB2) & ...
    Fit.nFitSamples>=P.minSamplesForNeuronFit;
labels(~valid)={'NotFit'};
ahs(~valid)={'AHS_not_fit'};
info=struct('alpha',P.class.alpha,'alphaBonferroni',NaN, ...
    'bonferroniMode','pending_pooled_across_all_included_cells', ...
    'nBonferroniTests',NaN,'nValidNeuronsThisFish',sum(valid), ...
    'classificationRule',['Pending: sign of b2 after pooled Bonferroni ' ...
    'correction across all included cells and morphs'], ...
    'pValueCaveat',['Conventional OLS p-values do not correct for temporal ' ...
    'autocorrelation; use continuous coefficients and blocked CV for primary inference.']);
end

function T = makeFitTable(session,candidateOrder,neuronID,F,classLabel,ahsLabel,info,fishVIF,fishConditionZ)
N=numel(neuronID);
T=table(session,candidateOrder,neuronID(:), ...
    F.prefRad,F.prefDeg,F.prefShiftRad,F.prefShiftDeg, ...
    F.a,F.d,F.b0,F.b1,F.b2,F.c0, ...
    F.seA,F.seD,F.seB0,F.sePrefRad,F.sePrefDeg,F.seB1,F.seB2,F.seC0, ...
    F.pA,F.pD,F.pHDJoint,F.fHDJoint,F.pB1,F.pB2,F.pC0, ...
    F.aStd,F.dStd,F.b0Std,F.b1Std,F.b2Std, ...
    F.partialR2HD,F.r2,F.adjR2,F.aic,F.bic, ...
    F.phaseOnlyA,F.phaseOnlyD,F.phaseOnlyB0,F.phaseOnlyC0, ...
    F.phaseOnlySeA,F.phaseOnlySeD,F.phaseOnlySeB0, ...
    F.phaseOnlySePrefRad,F.phaseOnlySePrefDeg, ...
    F.phaseOnlyPA,F.phaseOnlyPD,F.phaseOnlyFJoint,F.phaseOnlyPJoint,F.phaseOnlyPartialR2, ...
    F.phaseOnlyPrefRad,F.phaseOnlyPrefDeg, ...
    F.phaseOnlyAStd,F.phaseOnlyDStd,F.phaseOnlyB0Std, ...
    F.phaseOnlyR2,F.phaseOnlyAdjR2,F.phaseOnlyAIC,F.phaseOnlyBIC, ...
    F.behaviorB1,F.behaviorB2,F.behaviorC0, ...
    F.behaviorSeB1,F.behaviorSeB2,F.behaviorSeC0, ...
    F.behaviorPB1,F.behaviorPB2,F.behaviorPC0, ...
    F.behaviorB1Std,F.behaviorB2Std, ...
    F.behaviorR2,F.behaviorAdjR2,F.behaviorAIC,F.behaviorBIC, ...
    F.behaviorEpsilonStd,F.behaviorEpsilonRMS, ...
    F.fAHVBeyondPhase,F.pAHVBeyondPhase,F.partialR2AHVBeyondPhase, ...
    F.cvR2PhaseOnly,F.cvR2Behavior,F.cvR2Full, ...
    F.cvUniquePhase,F.cvDeltaR2Phase,F.cvUniqueAHV,F.cvDeltaR2AHV, ...
    F.cvSSEPhaseOnly,F.cvSSEBehavior,F.cvSSEFull, ...
    F.epsilonStd,F.epsilonRMS,F.designRank,F.conditionNumberRaw, ...
    F.conditionNumberStandardized,F.vifCosPhi,F.vifSinPhi,F.vifAbsAHV,F.vifAHV, ...
    F.maxVIF,F.nFitSamples,classLabel(:),ahsLabel(:), ...
    'VariableNames',{'session','candidateOrder','neuronID', ...
    'prefRad','prefDeg','prefShiftRad','prefShiftDeg', ...
    'a','d','b0','b1','b2','c0', ...
    'seA','seD','seB0','sePrefRad','sePrefDeg','seB1','seB2','seC0', ...
    'pA','pD','pHDJoint','fHDJoint','pB1','pB2','pC0', ...
    'a_std','d_std','b0_std','b1_std','b2_std', ...
    'partialR2HD','r2','adjR2','AIC','BIC', ...
    'phaseOnlyA','phaseOnlyD','phaseOnlyB0','phaseOnlyC0', ...
    'phaseOnlySeA','phaseOnlySeD','phaseOnlySeB0', ...
    'phaseOnlySePrefRad','phaseOnlySePrefDeg', ...
    'phaseOnlyPA','phaseOnlyPD','phaseOnlyFJoint','phaseOnlyPJoint','phaseOnlyPartialR2', ...
    'phaseOnlyPrefRad','phaseOnlyPrefDeg', ...
    'phaseOnlyA_std','phaseOnlyD_std','phaseOnlyB0_std', ...
    'phaseOnlyR2','phaseOnlyAdjR2','phaseOnlyAIC','phaseOnlyBIC', ...
    'behaviorB1','behaviorB2','behaviorC0', ...
    'behaviorSeB1','behaviorSeB2','behaviorSeC0', ...
    'behaviorPB1','behaviorPB2','behaviorPC0', ...
    'behaviorB1_std','behaviorB2_std', ...
    'behaviorR2','behaviorAdjR2','behaviorAIC','behaviorBIC', ...
    'behaviorEpsilonStd','behaviorEpsilonRMS', ...
    'fAHVBeyondPhase','pAHVBeyondPhase','partialR2AHVBeyondPhase', ...
    'cvR2PhaseOnly','cvR2Behavior','cvR2Full', ...
    'cvUniquePhase','cvDeltaR2Phase','cvUniqueAHV','cvDeltaR2AHV', ...
    'cvSSEPhaseOnly','cvSSEBehavior','cvSSEFull', ...
    'epsilonStd','epsilonRMS','designRank','conditionNumberRaw', ...
    'conditionNumberStandardized','vifCosPhi','vifSinPhi','vifAbsAHV','vifAHV', ...
    'maxVIF','nFitSamples','class','ahsLabel'});

T.pThresholdAHV=repmat(info.alphaBonferroni,N,1);
T.bonferroniNTests=repmat(info.nBonferroniTests,N,1);
T.nValidNeuronsThisFish=repmat(info.nValidNeuronsThisFish,N,1);
T.fishVIFCosPhi=repmat(fishVIF(1),N,1);
T.fishVIFSinPhi=nan(N,1);
T.fishVIFAbsAHV=repmat(fishVIF(2),N,1);
T.fishVIFAHV=repmat(fishVIF(3),N,1);
T.fishConditionNumberStandardized=repmat(fishConditionZ,N,1);
% Explicit fixed-preference names; legacy a/d/b0 columns are retained.
T.sourcePrefRad=T.prefRad;
T.sourcePrefDeg=T.prefDeg;
T.betaTheta=T.a;
T.seBetaTheta=T.seA;
T.pBetaTheta=T.pA;
T.betaTheta_std=T.a_std;
% Each table row contains one numeric column vector of aligned OOF samples.
T.cvObserved=F.cvObserved;
T.cvPredPhaseOnly=F.cvPredPhaseOnly;
T.cvPredBehavior=F.cvPredBehavior;
T.cvPredFull=F.cvPredFull;
% Explicit aliases make the three model identities unambiguous downstream.
T.phaseAHVR2=T.r2;
T.phaseAHVAdjR2=T.adjR2;
T.phaseAHVAIC=T.AIC;
T.phaseAHVBIC=T.BIC;
T.cvR2PhaseAHV=T.cvR2Full;
T.cvSSEPhaseAHV=T.cvSSEFull;
T.cvPredPhaseAHV=F.cvPredFull;
T.ahvOnlyR2=T.behaviorR2;
T.ahvOnlyAdjR2=T.behaviorAdjR2;
T.cvR2AHVOnly=T.cvR2Behavior;
T.cvSSEAHVOnly=T.cvSSEBehavior;
T.cvPredAHVOnly=F.cvPredBehavior;
T.cvUniqueAHVBeyondPhase=T.cvUniqueAHV;
T.cvDeltaR2AHVBeyondPhase=T.cvDeltaR2AHV;
% Backward-compatible aliases for scripts that call phase+AHV "original".
T.originalA=T.a;
T.originalD=T.d;
T.originalB0=T.b0;
T.originalB1=T.b1;
T.originalB2=T.b2;
T.originalC0=T.c0;
T.originalPrefRad=T.prefRad;
T.originalPrefDeg=T.prefDeg;
T.originalA_std=T.a_std;
T.originalD_std=T.d_std;
T.originalB0_std=T.b0_std;
T.originalB1_std=T.b1_std;
T.originalB2_std=T.b2_std;
T.originalR2=T.r2;
T.originalAdjR2=T.adjR2;
T.originalAIC=T.AIC;
T.originalBIC=T.BIC;
T.cvR2Original=T.cvR2Full;
T.cvSSEOriginal=T.cvSSEFull;
T.cvPredOriginal=F.cvPredFull;
end

function T = makeFishSummary(session,fitTable,D,VIF,conditionZ)
ok=isfinite(fitTable.b0);
cvPhasePair=ok & isfinite(fitTable.cvR2PhaseOnly) & isfinite(fitTable.cvR2Full);
if any(cvPhasePair)
    fractionFullBetterThanPhase=mean( ...
        fitTable.cvR2Full(cvPhasePair)>fitTable.cvR2PhaseOnly(cvPhasePair));
else
    fractionFullBetterThanPhase=NaN;
end
cvBehaviorPair=ok & isfinite(fitTable.cvR2Behavior) & isfinite(fitTable.cvR2Full);
if any(cvBehaviorPair)
    fractionFullBetterThanAHV=mean( ...
        fitTable.cvR2Full(cvBehaviorPair)>fitTable.cvR2Behavior(cvBehaviorPair));
else
    fractionFullBetterThanAHV=NaN;
end
T=table(string(session),sum(ok),D.nValidBouts,D.nAHVBouts, ...
    D.totalPositiveRotations,D.totalNegativeRotations,D.minDirectionalRotations, ...
    string(D.behaviorTimeSource),D.corrAbsAHV_AHV,conditionZ,max(VIF,[],'omitnan'), ...
    median(fitTable.phaseOnlyB0(ok),'omitnan'), ...
    median(fitTable.b0(ok),'omitnan'),median(fitTable.b1(ok),'omitnan'), ...
    median(fitTable.b2(ok),'omitnan'),median(hypot(fitTable.b1(ok),fitTable.b2(ok)),'omitnan'), ...
    median(fitTable.cvR2PhaseOnly(ok),'omitnan'), ...
    median(fitTable.cvR2Behavior(ok),'omitnan'),median(fitTable.cvR2Full(ok),'omitnan'), ...
    median(fitTable.cvUniquePhase(ok),'omitnan'),median(fitTable.cvDeltaR2Phase(ok),'omitnan'), ...
    median(fitTable.cvUniqueAHV(ok),'omitnan'),median(fitTable.cvDeltaR2AHV(ok),'omitnan'), ...
    fractionFullBetterThanPhase,fractionFullBetterThanAHV, ...
    median(fitTable.phaseOnlyPartialR2(ok),'omitnan'), ...
    median(fitTable.partialR2AHVBeyondPhase(ok),'omitnan'),median(fitTable.partialR2HD(ok),'omitnan'), ...
    'VariableNames',{'session','nFitNeurons','nBouts','nAHVBouts', ...
    'positiveRotations','negativeRotations','minDirectionalRotations', ...
    'behaviorTimeSource','corrAbsAHV_AHV','conditionNumberStandardized','maxVIF', ...
    'medianPhaseOnlyB0','medianPhaseAHVB0','medianPhaseAHVB1','medianPhaseAHVB2', ...
    'medianPhaseAHVModulationStrength','medianCvR2PhaseOnly','medianCvR2Behavior', ...
    'medianCvR2Full','medianCvUniquePhase','medianCvDeltaR2Phase', ...
    'medianCvUniqueAHV','medianCvDeltaR2AHV', ...
    'fractionNeuronsFullCVBetterThanPhaseOnly', ...
    'fractionNeuronsFullCVBetterThanAHVOnly', ...
    'medianPhaseOnlyPartialR2','medianPartialR2AHVBeyondPhase','medianPartialR2PhaseBeyondAHV'});

T.medianR2PhaseOnly=median(fitTable.phaseOnlyR2(ok),'omitnan');
T.medianR2Behavior=median(fitTable.behaviorR2(ok),'omitnan');
T.medianR2Full=median(fitTable.r2(ok),'omitnan');
T.medianAICPhaseOnly=median(fitTable.phaseOnlyAIC(ok),'omitnan');
T.medianAICBehavior=median(fitTable.behaviorAIC(ok),'omitnan');
T.medianAICFull=median(fitTable.AIC(ok),'omitnan');
T.medianBICPhaseOnly=median(fitTable.phaseOnlyBIC(ok),'omitnan');
T.medianBICBehavior=median(fitTable.behaviorBIC(ok),'omitnan');
T.medianBICFull=median(fitTable.BIC(ok),'omitnan');
T.medianBehaviorB1Std=median(fitTable.behaviorB1_std(ok),'omitnan');
T.medianBehaviorB2Std=median(fitTable.behaviorB2_std(ok),'omitnan');
T.medianR2PhaseAHV=T.medianR2Full;
T.medianR2AHVOnly=T.medianR2Behavior;
T.medianCvR2PhaseAHV=T.medianCvR2Full;
T.medianCvR2AHVOnly=T.medianCvR2Behavior;
T.medianR2Original=T.medianR2Full;
T.medianAICOriginal=T.medianAICFull;
T.medianBICOriginal=T.medianBICFull;
T.medianCvR2Original=T.medianCvR2Full;
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
    xlabel(ax,'Selected-cell phase-tuning preference (deg)'); ylabel(ax,'Cells'); grid(ax,'on');
end

title(tl,sprintf('%s | %s: ORI\\_V15 phase-tuning selection | %d/%d cells', ...
    P.morphDisplayName,sessionName,Q.nSelected,Q.nCandidates), ...
    'Interpreter','tex');
saveFigureLocal(fig,outDir,[sessionName '_ORI_V15_phase_tuning_selection'],P);
end

function plotRegressorQC(R,P,outDir)
idx=evenlySpacedIndices(numel(R.tCa),P.qc.maxTimePoints);
t=R.tCa(idx)/60;
fig=figure('Visible',P.figureVisible,'Color','w','Position',[100 100 1450 650]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); plot(ax,t,wrapToPiLocal(R.phaseCa(idx)),'k');
ylabel(ax,'phase (rad)'); title(ax,'Network phase'); grid(ax,'on');
ax=nexttile(tl); plot(ax,t,R.ahvCa(idx)*180/pi,'Color',[0.15 0.35 0.85]);
ylabel(ax,'Mei AHV proxy (deg)'); xlabel(ax,'Time (min)'); title(ax,'Signed AHV'); grid(ax,'on');
linkaxes(findall(fig,'Type','axes'),'x');
title(tl,sprintf('%s | %s: phase and AHV regressors', ...
    P.morphDisplayName,R.name),'Interpreter','none');
saveFigureLocal(fig,outDir,[R.name '_phase_AHV_regressor_QC'],P);
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
saveFigureLocal(fig,outDir,[R.name '_phase_AHV_predictor_correlation'],P);
end

function plotPerFishModelComparison(R,P,outDir)
T=R.fitTable;
fig=figure('Visible',P.figureVisible,'Color','w','Position',[40 60 1850 950]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

cvR2=[T.cvR2PhaseOnly,T.cvR2Full,T.cvR2Behavior];
inSampleR2=[T.phaseOnlyR2,T.r2,T.behaviorR2];
aic=[T.phaseOnlyAIC,T.AIC,T.behaviorAIC];
bic=[T.phaseOnlyBIC,T.BIC,T.behaviorBIC];

ax=nexttile(tl); threeModelPlot(ax,cvR2,'Blocked-CV R^2',true);
ax=nexttile(tl); threeModelPlot(ax,inSampleR2,'In-sample R^2',true);
ax=nexttile(tl); threeModelPlot(ax,aic-aic(:,1),'Delta AIC vs phase-only',true); yline(ax,0,'k:');
ax=nexttile(tl); threeModelPlot(ax,bic-bic(:,1),'Delta BIC vs phase-only',true); yline(ax,0,'k:');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAHV,40, ...
    [0.55 0.30 0.75],'1 - SSE_{phase+AHV}/SSE_{phase-only}','AHV beyond phase');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniquePhase,40, ...
    [0.25 0.25 0.25],'1 - SSE_{phase+AHV}/SSE_{AHV-only}','Phase beyond AHV');
title(tl,[P.morphDisplayName ' | ' R.name ...
    ': matched three-model comparison'],'Interpreter','none');
saveFigureLocal(fig,outDir,[R.name '_three_model_comparison'],P);
end

function plotAllFishModelComparison(T,P,outDir)
fig=figure('Visible',P.figureVisible,'Color','w','Position',[40 60 1850 950]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

cvR2=[T.cvR2PhaseOnly,T.cvR2Full,T.cvR2Behavior];
inSampleR2=[T.phaseOnlyR2,T.r2,T.behaviorR2];
aic=[T.phaseOnlyAIC,T.AIC,T.behaviorAIC];
bic=[T.phaseOnlyBIC,T.BIC,T.behaviorBIC];

ax=nexttile(tl); threeModelPlot(ax,cvR2,'Blocked-CV R^2',false);
ax=nexttile(tl); threeModelPlot(ax,inSampleR2,'In-sample R^2',false);
ax=nexttile(tl); threeModelPlot(ax,aic-aic(:,1),'Delta AIC vs phase-only',false); yline(ax,0,'k:');
ax=nexttile(tl); threeModelPlot(ax,bic-bic(:,1),'Delta BIC vs phase-only',false); yline(ax,0,'k:');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniqueAHV,60, ...
    [0.55 0.30 0.75],'1 - SSE_{phase+AHV}/SSE_{phase-only}','AHV beyond phase');
ax=nexttile(tl); contributionHistogram(ax,T.cvUniquePhase,60, ...
    [0.25 0.25 0.25],'1 - SSE_{phase+AHV}/SSE_{AHV-only}','Phase beyond AHV');
title(tl,[P.morphDisplayName ...
    ' | all fish: phase, phase + AHV, and AHV models']);
saveFigureLocal(fig,outDir,'ALL_FISH_three_model_comparison',P);
end

function plotAllFishCoefficientSummary(T,P,outDir)
vars={'phaseOnlyB0_std','b0_std','behaviorB1_std', ...
    'behaviorB2_std','cvR2PhaseOnly','cvR2Full','cvR2Behavior'};
labels={'Phase b0 std','Phase + AHV b0 std','AHV-only |AHV| std', ...
    'AHV-only signed AHV std','Phase blocked-CV R^2', ...
    'Phase + AHV blocked-CV R^2','AHV blocked-CV R^2'};
fig=figure('Visible',P.figureVisible,'Color','w','Position',[70 70 1750 1000]);
tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
for k=1:numel(vars)
    ax=nexttile(tl); x=T.(vars{k}); x=x(isfinite(x));
    histogram(ax,x,50,'FaceColor',[0.25 0.25 0.25],'EdgeColor','none');
    xline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,labels{k}); ylabel(ax,'Neurons'); grid(ax,'on');
end
title(tl,[P.morphDisplayName ...
    ' | all fish: three-model standardized coefficients and prediction']);
saveFigureLocal(fig,outDir,'ALL_FISH_three_model_distributions',P);
end

function plotFishLevelSummary(F,P,outDir)
n=height(F); x=(1:n)';
fig=figure('Visible',P.figureVisible,'Color','w','Position',[80 80 1450 850]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl); pairedFishThree(ax,x,F.medianR2PhaseOnly,F.medianR2Full, ...
    F.medianR2Behavior,'Median in-sample R^2');
ax=nexttile(tl); pairedFishThree(ax,x,F.medianCvR2PhaseOnly,F.medianCvR2Full, ...
    F.medianCvR2Behavior,'Median blocked-CV R^2');
ax=nexttile(tl); plot(ax,x,F.medianCvUniqueAHV,'o-','LineWidth',1.2, ...
    'MarkerFaceColor',[0.55 0.30 0.75]); yline(ax,0,'k:');
ylabel(ax,'Median unique AHV CV'); grid(ax,'on'); title(ax,'AHV beyond phase');
ax=nexttile(tl); plot(ax,x,F.medianCvUniquePhase,'o-','LineWidth',1.2, ...
    'MarkerFaceColor',[0.20 0.55 0.85]); yline(ax,0,'k:');
ylabel(ax,'Median unique phase CV'); grid(ax,'on'); title(ax,'Phase beyond AHV');

axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    xlim(axs(k),[0.5 n+0.5]); xticks(axs(k),x); xticklabels(axs(k),F.session); xtickangle(axs(k),45);
end
title(tl,[P.morphDisplayName ' | fish-level matched three-model comparison']);
saveFigureLocal(fig,outDir,'ALL_FISH_three_model_fish_level_summary',P);
end

function pairedFishThree(ax,x,phaseOnly,phaseAHV,ahv,yLabel)
offset=[-0.18 0 0.18];
values=[phaseOnly(:),phaseAHV(:),ahv(:)];
colors=modelColors();
hold(ax,'on');
for i=1:numel(x)
    vals=values(i,:);
    if all(isfinite(vals))
        plot(ax,x(i)+offset,vals,'-','Color',[.70 .70 .70], ...
            'HandleVisibility','off');
    end
end
names={'Phase','Phase + AHV','AHV'};
for j=1:size(values,2)
    plot(ax,x+offset(j),values(:,j),'o','MarkerFaceColor',colors(j,:), ...
        'MarkerEdgeColor','k','DisplayName',names{j});
end
ylabel(ax,yLabel); grid(ax,'on'); title(ax,[yLabel ': matched frames']);
legend(ax,'Location','best');
end

function threeModelPlot(ax,M,yLabel,showMatchedLines)
names={'Phase','Phase + AHV','AHV'};
colors=modelColors();
ok=all(isfinite(M),2);
M=M(ok,:);
hold(ax,'on');
if isempty(M)
    xlim(ax,[0.5 size(M,2)+0.5]); xticks(ax,1:size(M,2)); xticklabels(ax,names);
    ylabel(ax,yLabel); grid(ax,'on'); title(ax,'No complete matched fits');
    return;
end
if showMatchedLines
    for i=1:size(M,1)
        plot(ax,1:size(M,2),M(i,:),'-','Color',[.78 .78 .78], ...
            'LineWidth',0.4,'HandleVisibility','off');
    end
end
h=gobjects(size(M,2),1);
for j=1:size(M,2)
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
xlim(ax,[0.5 size(M,2)+0.5]); xticks(ax,1:size(M,2)); xticklabels(ax,names); xtickangle(ax,25);
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
% Phase-only purple; phase+AHV blue; AHV-only orange.
colors=[0.55 0.30 0.75; 0.15 0.55 0.90; 0.90 0.40 0.15];
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
three_model_fit_table=fitTable; %#ok<NASGU>
phase_tuning_selection=Tuning; %#ok<NASGU>
three_model_definition=struct( ...
    'phaseOnlyFormula',['betaTheta*cos(sourcePrefRad-phi)+c0; sourcePrefRad ' ...
    'from first circular harmonic of stored tuning_curves_phi_all_cells'], ...
    'ahvOnlyFormula','b1*abs(AHV)+b2*AHV+c0', ...
    'phaseAHVFormula',['betaTheta*cos(sourcePrefRad-phi)+' ...
    'b1*abs(AHV)+b2*AHV+c0; sourcePrefRad from first circular harmonic ' ...
    'of stored tuning_curves_phi_all_cells'], ...
    'fixedPreferredPhase', ...
    'sourcePrefRad from first circular harmonic of stored tuning_curves_phi_all_cells', ...
    'preferredPhaseSource', ...
    'first circular harmonic of tuning_curves_phi_all_cells', ...
    'legacyFieldMapping','a=betaTheta; d=0; b0=abs(betaTheta); prefRad=sourcePrefRad', ...
    'comparisonPolicy','Identical neurons, frames and blocked-CV folds for all three models', ...
    'classInfo',info); %#ok<NASGU>
save(out,'neuron_class_label','class_masks','class_candidate_indices', ...
    'class_original_roi_ids','source_candidate_path','source_behavior_path', ...
    'three_model_fit_table', ...
    'phase_tuning_selection','three_model_definition','-append');
fprintf('Saved three-model classified copy: %s\n',out);
end

function label=ternaryActivityScale(doZscore)
if doZscore, label='zActivity'; else, label='F_over_F0'; end
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
