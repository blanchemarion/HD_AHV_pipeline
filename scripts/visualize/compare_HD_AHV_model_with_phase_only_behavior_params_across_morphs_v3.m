function results = compare_HD_AHV_model_with_phase_only_behavior_params_across_morphs_v3(cfg)
%COMPARE_HD_AHV_MODEL_WITH_PHASE_ONLY_BEHAVIOR_PARAMS_ACROSS_MORPHS_V3
% Compare three explicitly named neural-activity models across morphs.
%
% This script reads the morph-level output produced by:
%   fit_neuron_HD_AHV_three_models_behavior_tuned_cells_only_ROC.m
%
% It expects one result file per morph:
%   HD_AHV_behavior_augmented_phase_tuned_only_with_phase_only_results.mat
%
% THREE PRIMARY MODELS USE THE SAME NEURONS, FRAMES AND CV FOLDS
% ----------------------------------------------------------------
%   phase + AHV: cos(phi)+sin(phi)+|AHV|+AHV
%   phase:       cos(phi)+sin(phi)
%   AHV:         |AHV|+AHV
%
% The saved variable names cvR2Original and cvR2Full are retained by the
% fitting script for backward compatibility. In every user-facing label,
% "Original" is called "phase + AHV" and "Full" is called
% "phase + AHV".
%
% The fitting script uses five contiguous blocked folds. This comparison
% script checks the saved fold count and reports aggregate held-out R^2.
% Folds are not treated as independent replicates.
%
% PRIMARY INFERENTIAL UNIT
% ------------------------
% Fish (true animal ID) is the independent biological replicate. Neurons are nested
% within fish and are therefore NOT entered as independent observations in
% cross-morph statistical tests. Each requested outcome is first summarized
% within fish.
%
% PRIMARY QUESTIONS
% -----------------
%   1. Do phase + AHV coefficients differ by morph?
%   2. How do the three primary models compare in held-out R^2?
%   3. Is phase prediction stronger than AHV-only prediction?
%   4. Does AHV add held-out information beyond phase?
%   5. Is phase predictive beyond the AHV-only model?
%   6. How does unique AHV contribution relate to phase-tuning consistency?
%
% IMPORTANT INTERPRETATION
% ------------------------
% * Because neural activity was z-scored identically within each neuron's fit
%   window, raw b0 is in activity-SD units.
% * Raw b1 and b2 are slopes per unit AHV (rad/s in the fitting script).
%   They are comparable as physical gain only if AHV units, filtering and
%   tau are identical for all morphs. This script checks saved settings.
% * Raw bF and bV depend on the scaling of amp70-weighted F and vigor.
%   Therefore standardized bF and bV are the primary cross-morph controls;
%   raw values are retained as sensitivity outputs.
% * Median signed b1 can be close to zero when positive and negative speed
%   effects coexist. Therefore median abs(b1) is added as a complementary
%   speed-encoding-strength metric.
% * Negative cross-validated R^2 and unique contributions are valid: they mean
%   worse than the test-set mean baseline. Values are not clipped at zero.
% * All finite, adequately fitted neurons selected as phase-tuned upstream
%   are included. No additional coefficient-significance filter is applied.
%
% STATISTICS
% ----------
% * Three pre-specified two-group Kruskal-Wallis comparisons for every plotted
%   fish-level metric:
%       Molino - Surface, Pachon - Surface.
% * Pairwise effect size: difference between fish-level medians, with a
%   fish-level bootstrap 95% CI, plus Cliff's delta.
% * No three-group omnibus test is run, so Molino-Pachon differences do not
%   contribute to the requested Surface-referenced tests.
% * The three planned pairwise contrasts are reported without multiplicity
%   correction by default, as explicitly configured below.
% * Within each morph, model scores are paired by fish and compared with
%   sign-rank tests. Cross-morph tests on matched unique-contribution
%   summaries provide a non-parametric model-by-morph sensitivity analysis.
%
% OUTPUTS
% -------
%   cross_morph_HD_AHV_phase_only_behavior_augmented_analysis.mat
%   README_analysis_notes.txt
%   fish_summary_raw_coefficients.csv
%   figures as PNG and SVG (no PDF files)

if nargin < 1 || isempty(cfg)
    cfg = pipeline_config();
end
assert(isstruct(cfg) && isfield(cfg, 'ModelDataDir'), ...
    'Pass the struct returned by pipeline_config.');

%% ======================== USER PARAMETERS ========================
P = struct();

% Source-output contract. These labels are used in figures, legends,
% statistical tables and the analysis notes. Keep this order everywhere.
P.sourceFittingScript = ...
    'fit_neuron_HD_AHV_three_models_behavior_tuned_cells_only_ROC.m';
P.modelNames = {'phase + AHV', 'AHV', 'phase'};
P.modelMetricNames = {'median_cvR2Original', 'median_cvR2Behavior', ...
    'median_cvR2PhaseOnly'};

% Input folders supplied by the user. Folder order does not determine plot
% order; P.morphOrder below does.
P.inputs = struct( ...
    'morph', {'Surface','Molino','Pachon'}, ...
    'folder', {cfg.ModelDataDir, cfg.ModelDataDir, cfg.ModelDataDir});

P.resultsFileName = 'HD_AHV_three_models_phase_tuned_only_results.mat';
for iInput = 1:numel(P.inputs)
    P.inputs(iInput).folder = cfg.ModelDataDir;
    P.inputs(iInput).fileName = sprintf( ...
        'HD_AHV_three_models_phase_tuned_only_results_%s.mat', ...
        lower(P.inputs(iInput).morph));
end

% Optional recording exclusions. Recording numbers are extracted from the
% saved session labels, so entries such as rec44, rec044 and rec44_behXX all
% match the value 44. Exclusion is applied before every fish summary, plot,
% ROC analysis and statistical test.
%
% IMPORTANT: excluding recordings because their behavior-only R^2 is high is
% outcome-dependent. Use this run as a clearly labelled sensitivity analysis,
% while retaining the all-fish run as the primary analysis unless there is an
% independent technical QC reason for exclusion.
P.exclusion.enabled = false;
P.exclusion.recNumbers = [44 70 52 80 81];
P.exclusion.reason = ...
    'Sensitivity analysis: unusually high behavior-only cross-validated R2';
P.exclusion.warnIfRequestedRecNotFound = true;

% Animal identity. Prefer an explicit animal/fish ID saved upstream. A mapping
% table can merge repeated recording IDs from the same animal. If neither is
% available, session is retained as recording ID with an explicit warning.
P.animalIdentity.sessionToAnimal = table();

% Output. A timestamp prevents accidental overwriting of an earlier run.
P.outputParent = cfg.ModelComparisonDir;
P.addTimestampToOutput = false;

% Morph order and colors: Surface red, Molino green, Pachon yellow.
P.morphOrder = {'Surface','Molino','Pachon'};
P.morphColors = [ ...
    0.85 0.20 0.20; ...
    0.20 0.65 0.30; ...
    0.95 0.72 0.10];

% Quality control. These filters concern numerical fit validity only, not
% coefficient significance or functional class.
P.quality.minFitSamples = 100;
P.quality.requireFullDesignRank = true;
P.quality.expectedDesignRank = 4; % 3 predictors plus intercept
P.quality.maxConditionNumber = 1e8;
P.quality.minNeuronsPerFish = 5;
P.quality.excludeByMaxVIF = false; % diagnostic by default, not an exclusion
P.quality.maxVIF = 10;

% Require the saved fitting run to use five contiguous blocked folds. The
% comparison script reads already-computed CV outputs; it does not refit.
P.cv.expectedBlockedFolds = 5;
P.cv.requireExpectedFoldCount = true;

% Fish-level model-performance exclusion.
% A fish is excluded from every downstream analysis when its median
% held-out R2 is below this threshold for at least one of the three models.
P.fishR2Exclusion.enabled = true;
P.fishR2Exclusion.threshold = -1;

% Statistics.
P.stats.alpha = 0.05;
P.stats.minFishPerMorph = 3;
P.stats.nBootstrap = 10000;
P.stats.bootstrapCI = [2.5 97.5];
P.stats.randomSeed = 7;
P.stats.pairwiseComparisons = { ...
    'Surface','Molino'; ...
    'Surface','Pachon'; ...
    'Molino','Pachon'};
% The three contrasts are explicit planned comparisons and are left
% uncorrected. Set this to 'bh' only if a corrected sensitivity analysis is
% desired; Holm correction is intentionally not implemented in this script.
P.stats.multipleComparisonMethod = 'none';

% Retained for saved-configuration compatibility. Three-group omnibus tests
% are disabled; only the three planned pairwise KW tests are run.
P.stats.correctOmnibusAcrossMetrics = false;

% Phase-dominance inference below uses held-out blocked-CV R^2 values.
P.stats.withinModelComparisons = { ...
    'median_cvR2Original','median_cvR2Behavior','phase + AHV','AHV'; ...
    'median_cvR2Original','median_cvR2PhaseOnly','phase + AHV','phase'; ...
    'median_cvR2Behavior','median_cvR2PhaseOnly','AHV','phase'};
P.stats.withinCoefficientComparisons = { ...
    'median_phaseOnlyB0_std','median_b0_std', ...
        'phase: standardized b0', ...
        'phase + AHV: standardized b0'; ...
    'median_phaseOnlyB0','median_originalB0', ...
        'phase: b0','phase + AHV: b0'; ...
    'median_phaseOnlyB0','median_b0_aug', ...
        'phase: b0','phase + AHV: b0'; ...
    'median_originalB0','median_b0_aug', ...
        'phase + AHV: b0','phase + AHV: b0'; ...
    'median_originalB1','median_b1_aug', ...
        'phase + AHV: signed b1', ...
        'phase + AHV: signed b1'; ...
    'median_abs_originalB1','median_abs_b1_aug', ...
        'phase + AHV: |b1|','phase + AHV: |b1|'; ...
    'median_originalB2','median_b2_aug', ...
        'phase + AHV: signed b2', ...
        'phase + AHV: signed b2'; ...
    'median_abs_originalB2','median_abs_b2_aug', ...
        'phase + AHV: |b2|','phase + AHV: |b2|'};

% Secondary class-stratified description. This is useful for asking whether
% an overall morph effect is due to class composition, but interpret it
% cautiously because class membership was derived from these coefficients.
P.classAnalysis.enabled = true;
P.classAnalysis.classes = {'CW','CCW','Symmetric'};
P.classAnalysis.minNeuronsPerFishClass = 5;

% Figure export.
P.figure.visible = cfg.FigureVisible;       % 'on' for interactive inspection
P.figure.savePNG = true;
P.figure.dpi = 300;
P.figure.pointSize = 46;
P.figure.jitterWidth = 0.11;
P.figure.morphBarWidth = 0.46;
P.figure.morphXLimits = [0.72 numel(P.morphOrder)+0.28];
P.ahvDistribution.nBins = 81;
P.ahvDistribution.tailPercentiles = [0.5 99.5];
P.ahvDistribution.nPermutations = 10000;

% Sample across the good-fit population and show distinct recording parts.
P.examples.numberPerMorph = 10;
P.examples.windowSeconds = 180;
P.examples.scoreQuantiles = linspace(0.92,0.38,P.examples.numberPerMorph);

rng(P.stats.randomSeed, 'twister');

%% ======================== METRIC DEFINITIONS ========================
% transform: identity | abs
% summary:   median | fraction_positive | mean
Metrics = [ ...
    ... % Primary raw coefficients from the fitted phase + AHV model
    makeMetric('median_b0_raw', 'b0_raw', 'identity','median', ...
        'Raw coefficient b_0 (fish median)', true, NaN), ...
    makeMetric('median_b1_raw', 'b1_raw', 'identity','median', ...
        'Raw coefficient b_1 (fish median)', true, 0), ...
    makeMetric('median_b2_raw', 'b2_raw', 'identity','median', ...
        'Raw coefficient b_2 (fish median)', true, 0), ...
    ... % phase model
    makeMetric('median_phaseOnlyB0',       'phaseOnlyB0', 'identity','median', ...
        'Median b_0: phase', true, NaN), ...
    makeMetric('median_phaseOnlyB0_std',   'phaseOnlyB0_std','identity','median', ...
        'Median standardized b_0: phase', true, NaN), ...
    ... % phase + AHV coefficients
    makeMetric('median_b0_aug',            'b0',          'identity','median', ...
        'Median b_0: phase + AHV', true, NaN), ...
    makeMetric('median_b1_aug',            'b1',          'identity','median', ...
        'Median signed b_1: phase + AHV', true, 0), ...
    makeMetric('median_abs_b1_aug',        'b1',          'abs','median', ...
        'Median |b_1|: phase + AHV', true, NaN), ...
    makeMetric('median_b2_aug',            'b2',          'identity','median', ...
        'Median signed b_2: phase + AHV', true, 0), ...
    makeMetric('median_abs_b2_aug',        'b2',          'abs','median', ...
        'Median |b_2|: phase + AHV', true, NaN), ...
    makeMetric('fraction_b2_positive_aug', 'b2',          'identity','fraction_positive', ...
        'Fraction b_2 > 0: phase + AHV', true, 0.5), ...
    ... % Exactly matched three-model five-fold contiguous blocked CV
    makeMetric('median_cvR2PhaseOnly',     'cvR2PhaseOnly','identity','median', ...
        'Median blocked-CV R^2: phase', true, 0), ...
    makeMetric('median_cvR2Original',      'cvR2Original','identity','median', ...
        'Median blocked-CV R^2: phase + AHV', true, 0), ...
    makeMetric('median_cvR2Behavior',      'cvR2Behavior','identity','median', ...
        'Median blocked-CV R^2: AHV-only', true, 0), ...
    makeMetric('median_cvR2Full',          'cvR2Full',    'identity','median', ...
        'Median blocked-CV R^2: phase + AHV', true, 0), ...
    makeMetric('median_deltaR2Phase','deltaR2Phase','identity','median', ...
        'Median Delta R^2 phase: phase + AHV minus AHV', true, 0), ...
    makeMetric('median_deltaR2Behavior','deltaR2Behavior','identity','median', ...
        'Median Delta R^2 AHV: phase + AHV minus phase', true, 0), ...
    makeMetric('median_phaseDominanceR2','phaseDominanceR2','identity','median', ...
        'Blocked-CV phase dominance (fish median; Delta R^2_phase minus Delta R^2_AHV)', true, 0), ...
    makeMetric('median_cvDeltaR2AddedBehavior','cvDeltaR2AddedBehavior','identity','median', ...
        'Median \DeltaCV R^2: adding forward + vigor to phase + AHV', true, 0), ...
    makeMetric('median_cvUniqueAddedBehavior','cvUniqueAddedBehavior','identity','median', ...
        'Median unique forward + vigor beyond phase + AHV', true, 0), ...
    makeMetric('median_cvDeltaR2Phase',    'cvDeltaR2Phase','identity','median', ...
        'Median ΔCV ^2: adding phase to AHV-only', true, 0), ...
    makeMetric('median_cvUniquePhase',     'cvUniquePhase','identity','median', ...
        'Median unique phase beyond AHV-only', true, 0), ...
    makeMetric('median_cvDeltaR2AHVBeyondPhase','cvDeltaR2AHVBeyondPhase', ...
        'identity','median', ...
        'Median DeltaCV-R^2: adding AHV to phase', true, 0), ...
    makeMetric('median_cvUniqueAHVBeyondPhase','cvUniqueAHVBeyondPhase', ...
        'identity','median', ...
        'Median unique AHV beyond phase', true, 0), ...
    makeMetric('median_cvDeltaR2AllBehaviorBeyondPhase', ...
        'cvDeltaR2AllBehaviorBeyondPhase','identity','median', ...
        'Median DeltaCV-R^2: adding AHV + forward + vigor to phase', true, 0), ...
    makeMetric('median_cvUniqueAllBehaviorBeyondPhase', ...
        'cvUniqueAllBehaviorBeyondPhase','identity','median', ...
        'Median unique AHV + forward + vigor beyond phase', true, 0), ...
    makeMetric('fraction_original_better_phaseOnly','originalBetterPhaseOnly', ...
        'identity','mean','Fraction neurons: phase + AHV > phase', false, 0.5), ...
    makeMetric('fraction_full_better_phaseOnly','fullBetterPhaseOnly', ...
        'identity','mean', ...
        'Fraction neurons: phase + AHV > phase', false, 0.5), ...
    makeMetric('fraction_full_better_original','fullBetterOriginal','identity','mean', ...
        'Fraction neurons: phase + AHV > phase + AHV', ...
        false, 0.5), ...
    makeMetric('fraction_full_better_behavior','fullBetterBehavior','identity','mean', ...
        'Fraction neurons: phase + AHV > AHV + forward + vigor', ...
        false, 0.5), ...
    makeMetric('fraction_cvR2PhaseOnly_positive','cvR2PhaseOnly','identity', ...
        'fraction_positive','Fraction phase CV-R^2 > 0', false, 0.5), ...
    makeMetric('fraction_cvR2Original_positive','cvR2Original','identity','fraction_positive', ...
        'Fraction phase + AHV CV-R^2 > 0', false, 0.5), ...
    makeMetric('fraction_cvR2Behavior_positive','cvR2Behavior','identity','fraction_positive', ...
        'Fraction AHV-only CV-R^2 > 0', false, 0.5), ...
    makeMetric('fraction_cvR2Full_positive','cvR2Full','identity','fraction_positive', ...
        'Fraction phase + AHV CV-R^2 > 0', false, 0.5), ...
    ... % Added behavioral controls; standardized values are primary
    makeMetric('median_bF_std',            'bF_std',      'identity','median', ...
        'Median standardized b_F', true, 0), ...
    makeMetric('median_abs_bF_std',        'bF_std',      'abs','median', ...
        'Median |standardized b_F|', true, NaN), ...
    makeMetric('median_bV_std',            'bV_std',      'identity','median', ...
        'Median standardized b_V', true, 0), ...
    makeMetric('median_abs_bV_std',        'bV_std',      'abs','median', ...
        'Median |standardized b_V|', true, NaN), ...
    makeMetric('median_bF_raw',            'bF',          'identity','median', ...
        'Median raw b_F', false, 0), ...
    makeMetric('median_abs_bF_raw',        'bF',          'abs','median', ...
        'Median |raw b_F|', false, NaN), ...
    makeMetric('median_bV_raw',            'bV',          'identity','median', ...
        'Median raw b_V', false, 0), ...
    makeMetric('median_abs_bV_raw',        'bV',          'abs','median', ...
        'Median |raw b_V|', false, NaN), ...
    ... % Paired coefficient sensitivity: adding forward + vigor
    makeMetric('median_originalB0',        'originalB0',  'identity','median', ...
        'Median b_0: phase + AHV', false, NaN), ...
    makeMetric('median_originalB1',        'originalB1',  'identity','median', ...
        'Median signed b_1: phase + AHV', false, 0), ...
    makeMetric('median_abs_originalB1',    'originalB1',  'abs','median', ...
        'Median |b_1|: phase + AHV', false, NaN), ...
    makeMetric('median_originalB2',        'originalB2',  'identity','median', ...
        'Median signed b_2: phase + AHV', false, 0), ...
    makeMetric('median_abs_originalB2',    'originalB2',  'abs','median', ...
        'Median |b_2|: phase + AHV', false, NaN), ...
    makeMetric('median_deltaB0_controls',  'deltaB0FromControls','identity','median', ...
        'Median Delta b_0 after adding forward + vigor', true, 0), ...
    makeMetric('median_deltaB1_controls',  'deltaB1FromControls','identity','median', ...
        'Median Delta b_1 after adding forward + vigor', true, 0), ...
    makeMetric('median_deltaAbsB1_controls','deltaAbsB1FromControls','identity','median', ...
        'Median Delta|b_1| after adding forward + vigor', true, 0), ...
    makeMetric('median_deltaB2_controls',  'deltaB2FromControls','identity','median', ...
        'Median Delta b_2 after adding forward + vigor', true, 0), ...
    makeMetric('median_deltaAbsB2_controls','deltaAbsB2FromControls','identity','median', ...
        'Median Delta|b_2| after adding forward + vigor', true, 0), ...
    makeMetric('median_abs_prefShift_original','absPrefShiftFromOriginalDeg','identity','median', ...
        'Median |preferred-phase shift| (deg)', false, 0), ...
    makeMetric('median_deltaB0_AHVBeyondPhase','deltaB0FromPhaseOnlyToOriginal', ...
        'identity','median','Median Delta b_0: phase + AHV minus phase', false, 0), ...
    makeMetric('median_deltaB0_AllBehaviorBeyondPhase','deltaB0FromPhaseOnlyToFull', ...
        'identity','median', ...
        'Median Delta b_0: phase + AHV minus phase', false, 0), ...
    makeMetric('median_abs_prefShift_AHVBeyondPhase', ...
        'absPrefShiftOriginalFromPhaseOnlyDeg','identity','median', ...
        'Median |phase shift|: phase + AHV vs phase (deg)', false, 0), ...
    makeMetric('median_abs_prefShift_AllBehaviorBeyondPhase', ...
        'absPrefShiftFullFromPhaseOnlyDeg','identity','median', ...
        'Median |phase shift|: phase + AHV vs phase (deg)', ...
        false, 0), ...
    ... % Cross-morph standardized phase/AHV sensitivity controls
    makeMetric('median_b0_std',            'b0_std',      'identity','median', ...
        'Median standardized b_0: phase + AHV', false, NaN), ...
    makeMetric('median_b1_std',            'b1_std',      'identity','median', ...
        'Median standardized b_1: phase + AHV', false, 0), ...
    makeMetric('median_abs_b1_std',        'b1_std',      'abs','median', ...
        'Median |standardized b_1|', false, NaN), ...
    makeMetric('median_b2_std',            'b2_std',      'identity','median', ...
        'Median standardized b_2: phase + AHV', false, 0), ...
    makeMetric('median_abs_b2_std',        'b2_std',      'abs','median', ...
        'Median |standardized b_2|', false, NaN), ...
    ... % Variance partitioning and diagnostics
    makeMetric('median_partialR2HD',       'partialR2HD','identity','median', ...
        'Median in-sample partial R^2 of phase', false, 0), ...
    makeMetric('median_phaseOnlyPartialR2','phaseOnlyPartialR2','identity','median', ...
        'Median in-sample phase R^2 vs intercept', false, 0), ...
    makeMetric('median_partialR2AHVBeyondPhase','partialR2AHVBeyondPhase', ...
        'identity','median','Median in-sample partial R^2 of AHV beyond phase', ...
        false, 0), ...
    makeMetric('median_partialR2AllBehaviorBeyondPhase', ...
        'partialR2AllBehaviorBeyondPhase','identity','median', ...
        'Median in-sample partial R^2 of all behavior beyond phase', false, 0), ...
    makeMetric('median_partialR2AddedBehavior','partialR2AddedBehavior','identity','median', ...
        'Median in-sample partial R^2 of F+V', false, 0), ...
    makeMetric('median_maxVIF',            'maxVIF','identity','median', ...
        'Median maximum predictor VIF', false, NaN), ...
    makeMetric('median_vifForward',        'vifForward','identity','median', ...
        'Median VIF of F', false, NaN), ...
    makeMetric('median_vifVigor',          'vifVigor','identity','median', ...
        'Median VIF of V', false, NaN), ...
    makeMetric('median_conditionNumber',   'conditionNumberStandardized','identity','median', ...
        'Median standardized condition number', false, NaN), ...
    ... % Interpretable positive/negative AHV gains
    makeMetric('median_gPositive',         'gPositive','identity','median', ...
        'Raw coefficient gain for positive AHV (fish median)', false, 0), ...
    makeMetric('median_gNegative',         'gNegative','identity','median', ...
        'Raw coefficient gain for negative AHV (fish median)', false, 0), ...
    makeMetric('median_overallRawAHVStrength', ...
        'overallRawAHVStrength','identity','median', ...
        'Overall raw AHV coefficient strength (fish median)', false, NaN), ...
    makeMetric('median_c0_raw', 'c0','identity','median', ...
        'Raw intercept c_0 (fish median)', false, 0), ...
    makeMetric('median_directionalAsymmetry','directionalAsymmetry','identity','median', ...
        'Median normalized directional asymmetry', false, 0)];

modelR2MetricNames = P.modelMetricNames;

% Only these outcomes are plotted and included in the revised inferential
% families. Keeping the full Metrics array above preserves all requested CSV
% summaries without generating inferential rows for unused diagnostics.
primaryRawCoefficientMetricNames = {'median_b0_raw','median_b1_raw', ...
    'median_b2_raw'};
incrementalR2MetricNames = {'median_deltaR2Phase','median_deltaR2Behavior'};
phaseDominanceMetricNames = {'median_phaseDominanceR2'};
classCoefficientMetricNames = primaryRawCoefficientMetricNames;
classPhaseDominanceMetricNames = phaseDominanceMetricNames;
classMetricNames = unique([classCoefficientMetricNames, ...
    classPhaseDominanceMetricNames],'stable');

StatFamilies = [ ...
    makeFamily('model_cv', modelR2MetricNames), ...
    makeFamily('raw_phase_AHV_coefficients', primaryRawCoefficientMetricNames), ...
    makeFamily('incremental_R2', incrementalR2MetricNames), ...
    makeFamily('phase_dominance_R2', phaseDominanceMetricNames)];
statsMetricNames = unique([StatFamilies.metrics], 'stable');
StatsMetrics = selectMetrics(Metrics, statsMetricNames);

%% ======================== LOAD ALL MORPHS ========================
fprintf('\nLoading matched three-model outputs from %s...\n', ...
    P.sourceFittingScript);
[AllNeurons, InputReport, ExclusionReport] = loadAllMorphResults(P);
assert(~isempty(AllNeurons), 'No neuron tables were loaded. Check P.inputs.');

AllNeurons = addDerivedNeuronVariables(AllNeurons);
requiredVariables = unique([{Metrics.source}, ...
    {'session','animalID','morph','fishKey'}]);
missing = setdiff(requiredVariables, AllNeurons.Properties.VariableNames);
assert(isempty(missing), 'AllNeurons is missing required variables: %s', ...
    strjoin(missing, ', '));

AllNeurons.baseQualityOK = makeBaseQualityMask(AllNeurons, P);

fprintf('Loaded %d neurons from %d fish.\n', height(AllNeurons), ...
    numel(unique(AllNeurons.fishKey)));
fprintf('Matched finite three-model CV outputs: %d/%d quality-passing neurons.\n', ...
    sum(AllNeurons.baseQualityOK & AllNeurons.cvQuartetOK), ...
    sum(AllNeurons.baseQualityOK));
for m = 1:numel(P.morphOrder)
    idx = strcmp(AllNeurons.morph, P.morphOrder{m});
    fprintf('  %-8s: %3d fish, %6d neurons\n', P.morphOrder{m}, ...
        numel(unique(AllNeurons.fishKey(idx))), sum(idx));
end

%% ======================== OUTPUT DIRECTORY ========================
if P.addTimestampToOutput
    outputDir = fullfile(P.outputParent, datestr(now, 'yyyymmdd_HHMMSS'));
else
    outputDir = P.outputParent;
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end
fprintf('Output directory:\n  %s\n', outputDir);

%% ======================== FISH-LEVEL SUMMARIES ========================
FishMetrics = buildFishSummaryTable(AllNeurons, Metrics, P, 'all');

% Exclude a complete fish if its fish-median held-out R2 is below the
% threshold for at least one of the three fitted models.
R2FishExclusionReport = table();

if P.fishR2Exclusion.enabled
    r2MetricNames = P.modelMetricNames;
    r2Values = FishMetrics{:,r2MetricNames};

    % NaN values do not trigger this particular exclusion. A fish is
    % excluded when at least one available model median is strictly < -1.
    failedModelMask = isfinite(r2Values) & ...
        r2Values < P.fishR2Exclusion.threshold;
    excludeFish = any(failedModelMask,2);

    reportVariables = [ ...
        {'morph','animalID','recordingID','fishKey'}, ...
        r2MetricNames];

    R2FishExclusionReport = ...
        FishMetrics(excludeFish,reportVariables);

    excludedRows = find(excludeFish);
    failedModels = cell(numel(excludedRows),1);

    for iExcluded = 1:numel(excludedRows)
        row = excludedRows(iExcluded);
        failedModels{iExcluded} = strjoin( ...
            P.modelNames(failedModelMask(row,:)),', ');
    end

    R2FishExclusionReport.failedModels = failedModels;
    R2FishExclusionReport.threshold = repmat( ...
        P.fishR2Exclusion.threshold, ...
        height(R2FishExclusionReport),1);

    excludedFishKeys = FishMetrics.fishKey(excludeFish);

    fprintf(['\nFish-level model-R2 exclusion: %d/%d fish excluded ' ...
        '(at least one median held-out R2 < %.3f).\n'], ...
        nnz(excludeFish),height(FishMetrics), ...
        P.fishR2Exclusion.threshold);

    if ~isempty(R2FishExclusionReport)
        disp(R2FishExclusionReport(:,{ ...
            'morph','animalID','recordingID', ...
            'median_cvR2Original','median_cvR2Behavior', ...
            'median_cvR2PhaseOnly','failedModels'}));
    end

    % Remove all neurons belonging to excluded animals. This propagates the
    % exclusion to ROC, class-stratified and any neuron-derived analyses.
    AllNeurons(ismember(AllNeurons.fishKey,excludedFishKeys),:) = [];

    % Remove the corresponding fish-level rows used by all primary plots
    % and statistical analyses.
    FishMetrics(excludeFish,:) = [];
end

assert(~isempty(FishMetrics), ...
    'All fish were excluded by the fish-level held-out R2 criterion.');

B1PositiveStats = runPositiveB1Tests(FishMetrics,P);
[RepresentativeExamples,AHVRegressorFish,AHVRegressorSamples] = ...
    buildRepresentativeExamplesAndAHVSummary(AllNeurons,P);
AHVRegressorStats = runAHVRegressorStats( ...
    AHVRegressorFish,AHVRegressorSamples,P);

FishRawCoefficientSummary = FishMetrics(:, {'morph','animalID', ...
    'recordingID','nNeuronsQuality','median_b0_raw','median_b1_raw', ...
    'median_b2_raw'});
FishRawCoefficientSummary.Properties.VariableNames{'nNeuronsQuality'} = ...
    'numberOfNeurons';
ModelDefinitions = buildModelDefinitionTable(P);

% Morph estimates are medians across fish with fish-level bootstrap CIs.
MorphEstimates = computeMorphEstimates(FishMetrics, Metrics, P);

% Inferential tests on fish-level values only and only for prespecified
% metrics that appear in figures.
[OmnibusStats, PairwiseStats] = runCrossMorphStats( ...
    FishMetrics, StatsMetrics, P, 'All selected phase-tuned neurons');
[OmnibusStats, PairwiseStats] = applyConfiguredCorrections( ...
    OmnibusStats, PairwiseStats, StatFamilies, P);
fprintf('\nPrimary planned two-group KW tests: raw coefficients, one fish median per animal\n');
rawPairRows = ismember(PairwiseStats.metric, primaryRawCoefficientMetricNames);
disp(PairwiseStats(rawPairRows, {'metricLabel','test','groupA','groupB','nA','nB','pRaw'}));

% Paired model tests use one aggregate CV score per fish. They do not treat
% neurons or folds as independent n. Model-difference metrics above are also
% compared across morphs and serve as interaction sensitivity analyses.
[WithinModelStats, ModelFriedmanStats] = runWithinMorphModelStats(FishMetrics, P);
WithinCoefficientStats = runWithinMorphPairedStats( ...
    FishMetrics,P.stats.withinCoefficientComparisons,P);
IncrementalR2Stats = runWithinMorphPairedStats(FishMetrics,{ ...
    'median_deltaR2Phase','median_deltaR2Behavior', ...
    'Delta R2 phase','Delta R2 behavior'},P);

%% ======================== PRIMARY FIGURES ========================
plotPairedModelCV(FishMetrics,WithinModelStats,ModelFriedmanStats,P, ...
    outputDir,'02_within_morph_paired_three_model_CV_R2');

plotRepresentativeExamples(RepresentativeExamples,P,outputDir, ...
    'representative_phase_tuned_neurons');
plotAHVRegressorDistributions(AHVRegressorFish,AHVRegressorSamples, ...
    AHVRegressorStats,P,outputDir, ...
    'cross_morph_abs_AHV_regressor_distribution');

plotFishLevelMetrics(FishMetrics, Metrics, primaryRawCoefficientMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    'Primary morph comparison: raw coefficients (fish medians)', ...
    '06_primary_raw_coefficients_fish_medians');

IncrementalR2PlotSummary = plotPairedIncrementalR2( ...
    FishMetrics,IncrementalR2Stats,P,outputDir, ...
    '10_paired_incremental_delta_R2_phase_behavior');

plotFishLevelMetrics(FishMetrics, Metrics, phaseDominanceMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    ['Cross-morph blocked-CV phase dominance: positive values indicate ' ...
        'stronger held-out phase than AHV-only prediction'], ...
    'cross_morph_phase_dominance_R2');

%% ======================== CLASS-STRATIFIED ANALYSIS ========================
FishClassMetrics = table();
ClassOmnibusStats = table();
ClassPairwiseStats = table();

if P.classAnalysis.enabled
    if ismember('class', AllNeurons.Properties.VariableNames)
        classMetrics = selectMetrics(Metrics, classMetricNames);

        FishClassMetrics = buildFishSummaryTable( ...
            AllNeurons, classMetrics, P, 'by_class');

        for c = 1:numel(P.classAnalysis.classes)
            cls = P.classAnalysis.classes{c};
            Fc = FishClassMetrics(strcmp(FishClassMetrics.class, cls), :);
            [Oc, Pc] = runCrossMorphStats(Fc, classMetrics, P, cls);
            classFamily = makeFamily(['class_' cls],classMetricNames);
            [Oc,Pc] = applyConfiguredCorrections(Oc,Pc,classFamily,P);
            ClassOmnibusStats = appendCompatibleTables(ClassOmnibusStats, Oc);
            ClassPairwiseStats = appendCompatibleTables(ClassPairwiseStats, Pc);
        end


        classCoefficientMetrics = selectMetrics(Metrics, ...
            classCoefficientMetricNames);
        plotClassStratifiedMetrics(FishClassMetrics, classCoefficientMetrics, ...
            ClassOmnibusStats, ClassPairwiseStats, P, outputDir, ...
            '10_class_stratified_coefficients');

        classPhaseDominanceMetrics = selectMetrics(Metrics, ...
            classPhaseDominanceMetricNames);
        plotClassStratifiedMetrics(FishClassMetrics, ...
            classPhaseDominanceMetrics, ...
            ClassOmnibusStats, ClassPairwiseStats, P, outputDir, ...
            '11_class_stratified_three_model_CV_R2');

    else
        warning('No class variable found. Skipping class-stratified analysis.');
    end
end

%% ======================== SAVE COMPLETE WORKSPACE ========================
IncrementalR2FigureMetadata = struct( ...
    'png',fullfile(outputDir,'10_paired_incremental_delta_R2_phase_behavior.png'), ...
    'svg',fullfile(outputDir,'10_paired_incremental_delta_R2_phase_behavior.svg'), ...
    'summaryCSV',fullfile(outputDir, ...
        '10_paired_incremental_delta_R2_phase_behavior_summary.csv'), ...
    'center','mean across fish-level neuron medians', ...
    'interval','fish-bootstrap 95% confidence interval', ...
    'tests','paired within-morph fish-level sign-rank tests');
save(fullfile(outputDir, ...
    'cross_morph_HD_AHV_phase_only_behavior_augmented_analysis.mat'), ...
    'P', 'Metrics', 'StatsMetrics', 'StatFamilies', ...
    'AllNeurons', 'FishMetrics', 'MorphEstimates', ...
    'OmnibusStats', 'PairwiseStats', 'FishClassMetrics', ...
    'ClassOmnibusStats', 'ClassPairwiseStats', 'WithinModelStats', ...
    'ModelFriedmanStats', 'WithinCoefficientStats', 'IncrementalR2Stats', ...
    'InputReport','ModelDefinitions','ExclusionReport', ...
    'R2FishExclusionReport','FishRawCoefficientSummary', ...
    'IncrementalR2PlotSummary','B1PositiveStats', ...
    'RepresentativeExamples','AHVRegressorFish','AHVRegressorSamples', ...
    'AHVRegressorStats', ...
    'IncrementalR2FigureMetadata','-v7.3');
writetable(FishRawCoefficientSummary, fullfile(outputDir, ...
    'fish_summary_raw_coefficients.csv'));
writetable(B1PositiveStats,fullfile(outputDir,"b1_positive_within_morph_tests.csv"));
writetable(AHVRegressorFish,fullfile(outputDir,"fish_AHV_summaries.csv"));
writetable(AHVRegressorStats,fullfile(outputDir,"cross_morph_AHV_summary_tests.csv"));
writetable(IncrementalR2PlotSummary,fullfile(outputDir, ...
    '10_paired_incremental_delta_R2_phase_behavior_summary.csv'));
if ~isempty(R2FishExclusionReport)
    writetable(R2FishExclusionReport,fullfile(outputDir, ...
        'fish_excluded_by_median_model_cvR2.csv'));
end

writeAnalysisNotes(outputDir, P);

fprintf('\nDone. Main inferential n is the number of fish, not neurons.\n');
fprintf('Results saved in:\n%s\n', outputDir);

results = struct('outputDir', outputDir, ...
    'allNeurons',AllNeurons,'fishMetrics', FishMetrics, 'morphEstimates', MorphEstimates, ...
    'omnibusStatistics', OmnibusStats, ...
    'pairwiseStatistics', PairwiseStats, ...
    'inputReport', InputReport, 'exclusionReport', ExclusionReport, ...
    'r2FishExclusionReport',R2FishExclusionReport, ...
    'b1PositiveStats',B1PositiveStats, ...
    'representativeExamples',RepresentativeExamples, ...
    'ahvRegressorFish',AHVRegressorFish,'ahvRegressorStats',AHVRegressorStats, ...
    'ahvRegressorSamples',AHVRegressorSamples, ...
    'ahvRegressorFigurePNG',fullfile(outputDir, ...
        'cross_morph_abs_AHV_regressor_distribution.png'), ...
    'ahvRegressorFigureSVG',fullfile(outputDir, ...
        'cross_morph_abs_AHV_regressor_distribution.svg'), ...
    'phaseDominanceStatistics',PairwiseStats(strcmp(PairwiseStats.metric, ...
        'median_phaseDominanceR2'),:), ...
    'phaseDominanceFigurePNG',fullfile(outputDir,'cross_morph_phase_dominance_R2.png'), ...
    'phaseDominanceFigureSVG',fullfile(outputDir,'cross_morph_phase_dominance_R2.svg'), ...
    'incrementalR2PlotSummary',IncrementalR2PlotSummary, ...
    'incrementalR2FigureMetadata',IncrementalR2FigureMetadata, ...
    'incrementalR2FigurePNG',fullfile(outputDir,'10_paired_incremental_delta_R2_phase_behavior.png'), ...
    'incrementalR2FigureSVG',fullfile(outputDir,'10_paired_incremental_delta_R2_phase_behavior.svg'), ...
    'incrementalR2SummaryCSV',fullfile(outputDir,'10_paired_incremental_delta_R2_phase_behavior_summary.csv'));
end

%% ======================== LOCAL FUNCTIONS ========================

function M = makeMetric(name, source, transform, summary, label, requested, referenceLine)
M = struct('name', name, 'source', source, 'transform', transform, ...
    'summary', summary, 'label', label, 'requested', requested, ...
    'referenceLine', referenceLine);
end

function T = buildModelDefinitionTable(P)
displayName = P.modelNames(:);
savedCVField = {'cvR2Original';'cvR2Behavior';'cvR2PhaseOnly'};
savedFitPrefix = {'original';'behavior';'phaseOnly'};
formula = { ...
    'cos(phi) + sin(phi) + |AHV| + AHV'; ...
    '|AHV| + AHV'; ...
    'cos(phi) + sin(phi)'};
sourceFittingScript = repmat({P.sourceFittingScript},numel(displayName),1);
T = table(displayName,savedFitPrefix,savedCVField,formula,sourceFittingScript);
end

function F = makeFamily(name, metricNames)
F = struct('name',name,'metrics',{metricNames});
end

function [AllNeurons, InputReport, ExclusionReport] = loadAllMorphResults(P)
AllNeurons = table();
ExclusionReport = emptyExclusionReport();

nInput = numel(P.inputs);
morph = cell(nInput,1);
resultFile = cell(nInput,1);
nFish = zeros(nInput,1);
nNeurons = zeros(nInput,1);
nFishBeforeExclusion = zeros(nInput,1);
nNeuronsBeforeExclusion = zeros(nInput,1);
nFishExcluded = zeros(nInput,1);
nNeuronsExcluded = zeros(nInput,1);
ahvTauSeconds = nan(nInput,1);
nBlockedCVFolds = nan(nInput,1);
minTestSamplesPerFold = nan(nInput,1);
minRotationsEachDirection = nan(nInput,1);
zscoreActivity = nan(nInput,1);
augmentedModelFormula = cell(nInput,1);
originalModelFormula = cell(nInput,1);
phaseOnlyModelFormula = cell(nInput,1);
behaviorModelFormula = cell(nInput,1);
sourceFittingScript = repmat({P.sourceFittingScript},nInput,1);
cvScheme = repmat({'contiguous blocked'},nInput,1);
animalIDSource = cell(nInput,1);

for i = 1:nInput
    morph{i} = P.inputs(i).morph;
    folder = P.inputs(i).folder;
    assert(exist(folder, 'dir') == 7, 'Input folder does not exist: %s', folder);
    if isfield(P.inputs, 'fileName') && ~isempty(P.inputs(i).fileName)
        preferredName = P.inputs(i).fileName;
    else
        preferredName = P.resultsFileName;
    end
    filePath = locateResultsFile(folder, preferredName);
    resultFile{i} = filePath;
    fprintf('  %-8s <- %s\n', morph{i}, filePath);

    S = load(filePath);
    if isfield(S, 'AllNeurons') && istable(S.AllNeurons)
        T = S.AllNeurons;
    elseif isfield(S, 'Results')
        T = concatenateResultFitTables(S.Results);
    else
        error('No AllNeurons table or Results cell found in %s', filePath);
    end
    validateThreeModelOutputTable(T,filePath,P);

    assert(ismember('session', T.Properties.VariableNames), ...
        'The table in %s has no session variable.', filePath);
    T.session = normalizeTextColumn(T.session);
    [T.animalID, animalIDSource{i}] = resolveAnimalIDs(T, P.animalIdentity);
    nFishBeforeExclusion(i) = numel(unique(T.session));
    nNeuronsBeforeExclusion(i) = height(T);

    [excludeMask, recNumber] = makeRecordingExclusionMask( ...
        T.session, P.exclusion);
    if any(excludeMask)
        excludedSessions = unique(T.session(excludeMask), 'stable');
        for e = 1:numel(excludedSessions)
            sessionMask = strcmp(T.session, excludedSessions{e});
            rec = unique(recNumber(sessionMask & excludeMask));
            rec = rec(isfinite(rec));
            if isempty(rec); rec = NaN; else; rec = rec(1); end
            newRow = table(morph(i), excludedSessions(e), ...
                strcat(morph(i), {'__'}, excludedSessions(e)), rec, ...
                sum(sessionMask), {P.exclusion.reason}, ...
                'VariableNames', ExclusionReport.Properties.VariableNames);
            ExclusionReport = [ExclusionReport; newRow]; %#ok<AGROW>
        end
        fprintf('  %-8s: excluded %d fish/session(s), %d neurons\n', ...
            morph{i}, numel(excludedSessions), sum(excludeMask));
        T(excludeMask,:) = [];
    end
    assert(~isempty(T), ...
        'All rows for %s were excluded. Check P.exclusion.recNumbers.', morph{i});

    T.morph = repmat(morph(i), height(T), 1);
    T.recordingID = T.session;
    T.fishKey = strcat(T.morph, {'__'}, T.animalID);

    % Normalize the class column when present.
    if ismember('class', T.Properties.VariableNames)
        T.class = normalizeTextColumn(T.class);
    end

    nFish(i) = numel(unique(T.fishKey));
    nNeurons(i) = height(T);
    nFishExcluded(i) = nFishBeforeExclusion(i) - nFish(i);
    nNeuronsExcluded(i) = nNeuronsBeforeExclusion(i) - nNeurons(i);

    if isfield(S, 'P')
        ahvTauSeconds(i) = getNestedNumeric(S.P, {'kernel','tauSeconds'});
        nBlockedCVFolds(i) = getNestedNumeric(S.P, {'cv','nBlockedFolds'});
        minTestSamplesPerFold(i) = getNestedNumeric(S.P, {'cv','minTestSamplesPerFold'});
        minRotationsEachDirection(i) = getNestedNumeric( ...
            S.P, {'lowTurnQC','minRotationsEachDirection'});
        zscoreActivity(i) = getNestedLogical(S.P, {'model','zscoreActivityWithinFitWindow'});
    end
    if isfield(S, 'Results')
        augmentedModelFormula{i} = getFirstFormulaFromResults( ...
            S.Results,{'augmentedModelFormula','modelFormula'});
        originalModelFormula{i} = getFormulaFromResults(S.Results,'originalModelFormula');
        phaseOnlyModelFormula{i} = getFormulaFromResults( ...
            S.Results,'phaseOnlyModelFormula');
        behaviorModelFormula{i} = getFormulaFromResults( ...
            S.Results,'behaviorModelFormula');
    else
        augmentedModelFormula{i} = '';
        originalModelFormula{i} = '';
        phaseOnlyModelFormula{i} = '';
        behaviorModelFormula{i} = '';
    end

    AllNeurons = appendCompatibleTables(AllNeurons, T);
end

InputReport = table(morph, resultFile, nFishBeforeExclusion, ...
    nNeuronsBeforeExclusion, nFishExcluded, nNeuronsExcluded, ...
    nFish, nNeurons, ahvTauSeconds, ...
    nBlockedCVFolds, minTestSamplesPerFold, minRotationsEachDirection, ...
    zscoreActivity, cvScheme, animalIDSource, sourceFittingScript, originalModelFormula, ...
    augmentedModelFormula, phaseOnlyModelFormula, behaviorModelFormula);

if P.exclusion.enabled && P.exclusion.warnIfRequestedRecNotFound
    found = unique(ExclusionReport.recNumber(isfinite(ExclusionReport.recNumber)));
    notFound = setdiff(unique(double(P.exclusion.recNumbers(:))), found);
    if ~isempty(notFound)
        warning('Requested recording exclusions not found in session labels: %s', ...
            mat2str(notFound(:).'));
    end
end

warnIfInconsistent(InputReport.ahvTauSeconds, 'AHV tau');
warnIfInconsistent(InputReport.nBlockedCVFolds, 'number of blocked CV folds');
warnIfInconsistent(InputReport.minTestSamplesPerFold, 'minimum test samples per fold');
warnIfInconsistent(InputReport.minRotationsEachDirection, ...
    'low-turn rotations-per-direction threshold');

warnIfInconsistent(InputReport.zscoreActivity, 'activity z-scoring setting');
finiteZ=InputReport.zscoreActivity(isfinite(InputReport.zscoreActivity));
if numel(finiteZ)~=nInput || numel(unique(finiteZ))~=1
    error(['Neural-response coefficients are not comparable: activity z-scoring ' ...
        'is missing or differs across morph files.']);
end
finiteTau = InputReport.ahvTauSeconds(isfinite(InputReport.ahvTauSeconds));
if numel(finiteTau) ~= nInput || numel(unique(finiteTau)) ~= 1
    error(['Raw AHV coefficients are not comparable: the saved AHV filter ' ...
        'tau is missing or differs across morphs.']);
end
if finiteZ(1)==1
    responseScale='per-neuron z-scored activity';
else
    responseScale='raw F/F0 activity (deltaFoF + 1)';
end
fprintf('Scaling check passed: %s; unstandardized AHV predictor with shared tau %.3g s.\n', ...
    responseScale,finiteTau(1));
nonemptyFormula = InputReport.augmentedModelFormula( ...
    ~cellfun(@isempty, InputReport.augmentedModelFormula));
if numel(unique(nonemptyFormula)) > 1
    warning(['Saved phase + AHV formulas are not ' ...
        'identical across morph files.']);
end
nonemptyFormula = InputReport.originalModelFormula( ...
    ~cellfun(@isempty,InputReport.originalModelFormula));
if numel(unique(nonemptyFormula)) > 1
    warning('Saved phase + AHV formulas are not identical across morph files.');
end
nonemptyFormula = InputReport.phaseOnlyModelFormula( ...
    ~cellfun(@isempty,InputReport.phaseOnlyModelFormula));
if numel(unique(nonemptyFormula)) > 1
    warning('Saved phase formulas are not identical across morph files.');
end
nonemptyFormula = InputReport.behaviorModelFormula( ...
    ~cellfun(@isempty,InputReport.behaviorModelFormula));
if numel(unique(nonemptyFormula)) > 1
    warning('Saved AHV formulas are not identical across morph files.');
end
if any(contains(lower(string(nonemptyFormula)), {'forward','vigor'}), 'all')
    error('Loaded outputs still contain forward or vigor in the AHV-only model. Refit the three models.');
end
finiteFolds = InputReport.nBlockedCVFolds(isfinite(InputReport.nBlockedCVFolds));
if P.cv.requireExpectedFoldCount && ...
        (numel(finiteFolds) ~= nInput || any(finiteFolds ~= P.cv.expectedBlockedFolds))
    error(['Expected every morph file to contain results from %d blocked CV folds. ' ...
        'Inspect the in-memory InputReport table or set ' ...
        'P.cv.requireExpectedFoldCount=false for exploratory loading.'], ...
        P.cv.expectedBlockedFolds);
end
end

function [animalID, source] = resolveAnimalIDs(T, A)
% Prefer a true animal identifier; use a supplied map to merge recordings.
candidates = {'animalID','animalId','fishID','fishId','animal','fish'};
source = '';
for k = 1:numel(candidates)
    if ismember(candidates{k}, T.Properties.VariableNames)
        animalID = normalizeTextColumn(T.(candidates{k}));
        source = ['table.' candidates{k}];
        return;
    end
end
if isfield(A,'sessionToAnimal') && istable(A.sessionToAnimal) && ~isempty(A.sessionToAnimal)
    M = A.sessionToAnimal;
    assert(all(ismember({'session','animalID'},M.Properties.VariableNames)), ...
        'Animal map needs session and animalID variables.');
    ms = normalizeTextColumn(M.session); ma = normalizeTextColumn(M.animalID);
    assert(numel(unique(ms)) == numel(ms),'Duplicate session in animal-ID map.');
    animalID = cell(height(T),1);
    for i = 1:height(T)
        hit = find(strcmp(ms,T.session{i}),1);
        assert(~isempty(hit),'Session %s is absent from the animal-ID map.',T.session{i});
        animalID{i} = ma{hit};
    end
    source = 'P.animalIdentity.sessionToAnimal';
    return;
end
animalID = T.session;
source = 'session/recording fallback';
warning(['No explicit animal ID is present; recording ID is used as animal ID. ' ...
    'If an animal has repeated recordings, fill P.animalIdentity.sessionToAnimal.']);
end

function T = emptyExclusionReport()
T = table(cell(0,1), cell(0,1), cell(0,1), zeros(0,1), ...
    zeros(0,1), cell(0,1), ...
    'VariableNames', {'morph','session','fishKey','recNumber', ...
    'nNeuronsExcluded','reason'});
end

function [excludeMask, recNumber] = makeRecordingExclusionMask(session, E)
% Extract the first integer following "rec" (case-insensitive). This avoids
% partial-string errors such as rec4 matching rec44.
session = normalizeTextColumn(session);
recNumber = nan(numel(session),1);
for k = 1:numel(session)
    token = regexp(session{k}, '(?i)rec[_\-\s]*0*(\d+)', ...
        'tokens', 'once');
    if ~isempty(token)
        recNumber(k) = str2double(token{1});
    end
end

excludeMask = false(numel(session),1);
if ~E.enabled || isempty(E.recNumbers)
    return;
end
requested = unique(double(E.recNumbers(:)));
requested = requested(isfinite(requested) & requested >= 0 & ...
    requested == floor(requested));
excludeMask = ismember(recNumber, requested);
end

function filePath = locateResultsFile(folder, preferredName)
candidate = fullfile(folder, preferredName);
if exist(candidate, 'file') == 2
    filePath = candidate;
    return;
end

files = dir(fullfile(folder, '*.mat'));
matches = {};
for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);
    vars = whos('-file', p);
    if any(strcmp({vars.name}, 'AllNeurons')) || any(strcmp({vars.name}, 'Results'))
        matches{end+1,1} = p; %#ok<AGROW>
    end
end
assert(~isempty(matches), 'No compatible result MAT file found in %s', folder);
assert(numel(matches) == 1, ...
    'Several compatible MAT files found in %s. Set P.resultsFileName explicitly.', folder);
filePath = matches{1};
end

function T = concatenateResultFitTables(Results)
T = table();
for i = 1:numel(Results)
    R = Results{i};
    if isstruct(R) && isfield(R, 'fitTable') && istable(R.fitTable) && ...
            height(R.fitTable) > 0
        T = appendCompatibleTables(T, R.fitTable);
    end
end
end

function out = normalizeTextColumn(in)
if iscell(in)
    out = cellfun(@localToChar, in, 'UniformOutput', false);
elseif isstring(in)
    out = cellstr(in);
elseif iscategorical(in)
    out = cellstr(in);
elseif ischar(in)
    if size(in,1) == 1
        out = {in};
    else
        out = cellstr(in);
    end
else
    out = arrayfun(@(x) localToChar(x), in, 'UniformOutput', false);
end
out = out(:);
end

function s = localToChar(x)
if ischar(x)
    s = x;
elseif isstring(x)
    s = char(x);
elseif iscategorical(x)
    s = char(x);
elseif isnumeric(x) || islogical(x)
    s = num2str(x);
else
    s = char(string(x));
end
end

function value = getNestedNumeric(S, fields)
value = NaN;
x = S;
for i = 1:numel(fields)
    if ~isstruct(x) || ~isfield(x, fields{i})
        return;
    end
    x = x.(fields{i});
end
if isnumeric(x) && isscalar(x)
    value = double(x);
end
end

function value = getNestedLogical(S, fields)
value = NaN;
x = S;
for i = 1:numel(fields)
    if ~isstruct(x) || ~isfield(x, fields{i})
        return;
    end
    x = x.(fields{i});
end
if (islogical(x) || isnumeric(x)) && isscalar(x)
    value = double(x);
end
end

function formula = getFormulaFromResults(Results, fieldName)
formula = '';
for i = 1:numel(Results)
    R = Results{i};
    if isstruct(R) && isfield(R, 'fit') && isfield(R.fit, fieldName)
        formula = R.fit.(fieldName);
        return;
    end
end
end

function formula = getFirstFormulaFromResults(Results,fieldNames)
formula = '';
for k = 1:numel(fieldNames)
    formula = getFormulaFromResults(Results,fieldNames{k});
    if ~isempty(formula)
        return;
    end
end
end

function validateThreeModelOutputTable(T,filePath,P)
% Require the three fitted models and their matched CV outputs.
% versions that computed AHV + forward + vigor only inside CV.
required = { ...
    'phaseOnlyR2','phaseOnlyAIC','phaseOnlyBIC', ...
    'originalR2','originalAIC','originalBIC', ...
    'r2','AIC','BIC', ...
    'behaviorR2','behaviorAIC','behaviorBIC', ...
    'behaviorB1','behaviorB2', ...
    'cvR2Original','cvR2PhaseOnly','cvR2Behavior', ...
    'cvSSEOriginal','cvSSEPhaseOnly','cvSSEBehavior'};
missing = setdiff(required,T.Properties.VariableNames);
if ~isempty(missing)
    error(['The file %s is not a complete output of %s. Missing fields: %s. ' ...
        'Re-run the three-model fitting script for this morph.'], ...
        filePath,P.sourceFittingScript,strjoin(missing,', '));
end
end

function warnIfInconsistent(x, label)
x = x(isfinite(x));
if numel(unique(x)) > 1
    warning('%s differs across morph result files.', label);
end
end

function T = appendCompatibleTables(T, U)
if isempty(T)
    T = U;
    return;
end

% All fitted morph outputs should have identical tables. This explicit check
% produces a useful error instead of an opaque vertical-concatenation error.
if ~isequal(T.Properties.VariableNames, U.Properties.VariableNames)
    missingInU = setdiff(T.Properties.VariableNames, U.Properties.VariableNames);
    extraInU = setdiff(U.Properties.VariableNames, T.Properties.VariableNames);
    error(['Cannot concatenate tables. Missing in new table: %s. ' ...
        'Extra in new table: %s.'], strjoin(missingInU, ', '), ...
        strjoin(extraInU, ', '));
end
T = [T; U]; %#ok<AGROW>
end

function ok = makeBaseQualityMask(T, P)
ok = true(height(T),1);

if ismember('nFitSamples', T.Properties.VariableNames)
    ok = ok & isfinite(T.nFitSamples) & T.nFitSamples >= P.quality.minFitSamples;
end
if P.quality.requireFullDesignRank && ...
        ismember('designRank', T.Properties.VariableNames)
    ok = ok & isfinite(T.designRank) & ...
        T.designRank >= P.quality.expectedDesignRank;
end
if ismember('conditionNumberStandardized', T.Properties.VariableNames) && ...
        isfinite(P.quality.maxConditionNumber)
    ok = ok & isfinite(T.conditionNumberStandardized) & ...
        T.conditionNumberStandardized <= P.quality.maxConditionNumber;
end
if P.quality.excludeByMaxVIF && ismember('maxVIF',T.Properties.VariableNames)
    ok = ok & isfinite(T.maxVIF) & T.maxVIF <= P.quality.maxVIF;
end
end

function T = addDerivedNeuronVariables(T)
% The current fitter has exactly three models (phase, AHV, phase + AHV).
% Retain legacy report columns explicitly: with no forward/vigor terms,
% all behavior is AHV and the effect of adding extra behavior is zero.
if ~ismember('cvUniqueAllBehaviorBeyondPhase', T.Properties.VariableNames)
    T.cvUniqueAllBehaviorBeyondPhase = T.cvUniqueAHVBeyondPhase;
end
if ~ismember('cvDeltaR2AllBehaviorBeyondPhase', T.Properties.VariableNames)
    T.cvDeltaR2AllBehaviorBeyondPhase = T.cvDeltaR2AHVBeyondPhase;
end
if ~ismember('cvUniqueAddedBehavior', T.Properties.VariableNames)
    T.cvUniqueAddedBehavior = zeros(height(T),1);
end
if ~ismember('cvDeltaR2AddedBehavior', T.Properties.VariableNames)
    T.cvDeltaR2AddedBehavior = zeros(height(T),1);
end
if ~ismember('partialR2AllBehaviorBeyondPhase', T.Properties.VariableNames)
    T.partialR2AllBehaviorBeyondPhase = T.partialR2AHVBeyondPhase;
end
if ~ismember('partialR2AddedBehavior', T.Properties.VariableNames)
    T.partialR2AddedBehavior = zeros(height(T),1);
end
zeroColumns = {'bF','bV','bF_std','bV_std', ...
    'deltaB0FromControls','deltaB1FromControls','deltaB2FromControls'};
for iColumn = 1:numel(zeroColumns)
    if ~ismember(zeroColumns{iColumn}, T.Properties.VariableNames)
        T.(zeroColumns{iColumn}) = zeros(height(T),1);
    end
end
nanColumns = {'vifForward','vifVigor'};
for iColumn = 1:numel(nanColumns)
    if ~ismember(nanColumns{iColumn}, T.Properties.VariableNames)
        T.(nanColumns{iColumn}) = nan(height(T),1);
    end
end

% Primary raw coefficients use the unstandardized phase + AHV fit.
T.a_raw = T.originalA;
T.d_raw = T.originalD;
T.b0_raw = hypot(T.a_raw,T.d_raw);
T.b1_raw = T.originalB1;
T.b2_raw = T.originalB2;
% These are deterministic re-expressions of already saved matched fits.
required = {'b1','b2','b1_std','b2_std','originalB1','originalB2','phaseOnlyB0', ...
    'originalB0','b0','phaseOnlyPrefRad','prefRad','originalPrefRad', ...
    'cvR2PhaseOnly','cvR2Original','cvR2Behavior','cvR2Full', ...
    'cvUniqueAHVBeyondPhase','cvDeltaR2AHVBeyondPhase', ...
    'cvUniqueAllBehaviorBeyondPhase','cvDeltaR2AllBehaviorBeyondPhase'};
missing = setdiff(required,T.Properties.VariableNames);
assert(isempty(missing),'Cannot derive comparison variables; missing: %s', ...
    strjoin(missing,', '));

T.gPositive = T.b1 + T.b2;
T.gNegative = T.b1 - T.b2;
T.directionalStrength = (abs(T.gPositive) + abs(T.gNegative)) ./ 2;
T.overallRawAHVStrength = hypot(T.b1,T.b2);
T.overallStandardizedAHVStrength = hypot(T.b1_std,T.b2_std);
T.deltaR2Phase = T.cvR2Original - T.cvR2Behavior;
T.deltaR2Behavior = T.cvR2Original - T.cvR2PhaseOnly;
T.phaseDominanceR2 = T.cvR2PhaseOnly - T.cvR2Behavior;
den = abs(T.gPositive) + abs(T.gNegative);
T.directionalAsymmetry = (T.gPositive - T.gNegative) ./ den;
T.directionalAsymmetry(~isfinite(T.directionalAsymmetry) | den <= eps) = NaN;

T.originalGPositive = T.originalB1 + T.originalB2;
T.originalGNegative = T.originalB1 - T.originalB2;
T.deltaAbsB1FromControls = abs(T.b1) - abs(T.originalB1);
T.deltaAbsB2FromControls = abs(T.b2) - abs(T.originalB2);
T.absPrefShiftFromOriginalDeg = abs(rad2deg(wrapToPiLocal( ...
    T.prefRad - T.originalPrefRad)));
T.deltaB0FromPhaseOnlyToOriginal = T.originalB0 - T.phaseOnlyB0;
T.deltaB0FromPhaseOnlyToFull = T.b0 - T.phaseOnlyB0;
T.absPrefShiftOriginalFromPhaseOnlyDeg = abs(rad2deg(wrapToPiLocal( ...
    T.originalPrefRad - T.phaseOnlyPrefRad)));
T.absPrefShiftFullFromPhaseOnlyDeg = abs(rad2deg(wrapToPiLocal( ...
    T.prefRad - T.phaseOnlyPrefRad)));

T.cvQuartetOK = isfinite(T.cvR2PhaseOnly) & isfinite(T.cvR2Original) & ...
    isfinite(T.cvR2Behavior) & isfinite(T.cvR2Full);
T.originalBetterPhaseOnly = nan(height(T),1);
T.fullBetterPhaseOnly = nan(height(T),1);
T.fullBetterOriginal = nan(height(T),1);
T.fullBetterBehavior = nan(height(T),1);
T.behaviorBetterOriginal = nan(height(T),1);
T.originalBetterPhaseOnly(T.cvQuartetOK) = double( ...
    T.cvR2Original(T.cvQuartetOK) > T.cvR2PhaseOnly(T.cvQuartetOK));
T.fullBetterPhaseOnly(T.cvQuartetOK) = double( ...
    T.cvR2Full(T.cvQuartetOK) > T.cvR2PhaseOnly(T.cvQuartetOK));
T.fullBetterOriginal(T.cvQuartetOK) = double( ...
    T.cvR2Full(T.cvQuartetOK) > T.cvR2Original(T.cvQuartetOK));
T.fullBetterBehavior(T.cvQuartetOK) = double( ...
    T.cvR2Full(T.cvQuartetOK) > T.cvR2Behavior(T.cvQuartetOK));
T.behaviorBetterOriginal(T.cvQuartetOK) = double( ...
    T.cvR2Behavior(T.cvQuartetOK) > T.cvR2Original(T.cvQuartetOK));
end

function Fish = buildFishSummaryTable(All, Metrics, P, mode)
fishKeys = unique(All.fishKey, 'stable');

if strcmp(mode, 'by_class')
    nRows = numel(fishKeys) * numel(P.classAnalysis.classes);
else
    nRows = numel(fishKeys);
end

morph = cell(nRows,1);
session = cell(nRows,1);
animalID = cell(nRows,1);
recordingID = cell(nRows,1);
fishKey = cell(nRows,1);
classLabel = cell(nRows,1);
nNeuronsAvailable = zeros(nRows,1);
nNeuronsQuality = zeros(nRows,1);

Fish = table(morph, session, animalID, recordingID, fishKey, classLabel, ...
    nNeuronsAvailable, nNeuronsQuality, ...
    'VariableNames', {'morph','session','animalID','recordingID','fishKey','class', ...
    'nNeuronsAvailable','nNeuronsQuality'});
for k = 1:numel(Metrics)
    Fish.(Metrics(k).name) = nan(nRows,1);
    Fish.(['n_' Metrics(k).name]) = zeros(nRows,1);
end

row = 0;
for f = 1:numel(fishKeys)
    idxFish = strcmp(All.fishKey, fishKeys{f});
    if strcmp(mode, 'by_class')
        levels = P.classAnalysis.classes;
    else
        levels = {''};
    end

    for c = 1:numel(levels)
        row = row + 1;
        idx = idxFish;
        if strcmp(mode, 'by_class')
            idx = idx & strcmp(All.class, levels{c});
        end

        first = find(idxFish, 1, 'first');
        Fish.morph{row} = All.morph{first};
        Fish.session{row} = All.session{first};
        Fish.animalID{row} = All.animalID{first};
        Fish.recordingID{row} = strjoin(unique(All.recordingID(idxFish),'stable'),';');
        Fish.fishKey{row} = All.fishKey{first};
        Fish.class{row} = levels{c};
        Fish.nNeuronsAvailable(row) = sum(idx);
        Fish.nNeuronsQuality(row) = sum(idx & All.baseQualityOK);

        if strcmp(mode, 'by_class')
            minimumN = P.classAnalysis.minNeuronsPerFishClass;
        else
            minimumN = P.quality.minNeuronsPerFish;
        end

        for k = 1:numel(Metrics)
            values = extractNeuronValues(All, Metrics(k), idx & All.baseQualityOK);
            Fish.(['n_' Metrics(k).name])(row) = numel(values);
            if numel(values) >= minimumN
                Fish.(Metrics(k).name)(row) = summarizeValues(values, Metrics(k).summary);
            end
        end
    end
end

% In all-neuron mode the class column is redundant; remove it from the CSV.
if strcmp(mode, 'all')
    Fish.class = [];
end
end

function values = extractNeuronValues(T, M, rowMask)
raw = double(T.(M.source));
valid = rowMask(:) & isfinite(raw);
cvSources = {'cvR2PhaseOnly','cvR2Original','cvR2Behavior','cvR2Full', ...
    'cvUniquePhase','cvDeltaR2Phase','cvUniqueAddedBehavior', ...
    'cvDeltaR2AddedBehavior','cvUniqueAHVBeyondPhase', ...
    'cvDeltaR2AHVBeyondPhase','cvUniqueAllBehaviorBeyondPhase','phaseDominanceR2', ...
    'cvDeltaR2AllBehaviorBeyondPhase','deltaR2Phase','deltaR2Behavior', ...
    'originalBetterPhaseOnly', ...
    'fullBetterPhaseOnly','fullBetterOriginal','fullBetterBehavior', ...
    'behaviorBetterOriginal'};
if ismember(M.source,cvSources) && ismember('cvQuartetOK',T.Properties.VariableNames)
    valid = valid & T.cvQuartetOK;
end
values = raw(valid);
switch M.transform
    case 'identity'
        % no change
    case 'abs'
        values = abs(values);
    otherwise
        error('Unknown metric transform: %s', M.transform);
end
values = values(:);
end

function value = summarizeValues(x, method)
switch method
    case 'median'
        value = median(x, 'omitnan');
    case 'fraction_positive'
        value = mean(x > 0);
    case 'mean'
        value = mean(x,'omitnan');
    otherwise
        error('Unknown fish summary method: %s', method);
end
end

function MorphEstimates = computeMorphEstimates(Fish, Metrics, P)
metric = {};
metricLabel = {};
morph = {};
nFish = [];
estimate = [];
ciLow = [];
ciHigh = [];

for k = 1:numel(Metrics)
    for m = 1:numel(P.morphOrder)
        x = Fish.(Metrics(k).name)(strcmp(Fish.morph, P.morphOrder{m}));
        x = x(isfinite(x));
        [est, ci] = bootstrapMedianCI(x, P.stats.nBootstrap, P.stats.bootstrapCI);
        metric{end+1,1} = Metrics(k).name; %#ok<AGROW>
        metricLabel{end+1,1} = Metrics(k).label; %#ok<AGROW>
        morph{end+1,1} = P.morphOrder{m}; %#ok<AGROW>
        nFish(end+1,1) = numel(x); %#ok<AGROW>
        estimate(end+1,1) = est; %#ok<AGROW>
        ciLow(end+1,1) = ci(1); %#ok<AGROW>
        ciHigh(end+1,1) = ci(2); %#ok<AGROW>
    end
end

MorphEstimates = table(metric, metricLabel, morph, nFish, estimate, ciLow, ciHigh);
end

function [Omnibus, Pairwise] = runCrossMorphStats(Fish, Metrics, P, scope)
metric = cell(numel(Metrics),1);
metricLabel = cell(numel(Metrics),1);
analysisScope = repmat({scope}, numel(Metrics),1);
nTotal = zeros(numel(Metrics),1);
nSurface = zeros(numel(Metrics),1);
nMolino = zeros(numel(Metrics),1);
nPachon = zeros(numel(Metrics),1);
H = nan(numel(Metrics),1);
df = nan(numel(Metrics),1);
pRaw = nan(numel(Metrics),1);

for k = 1:numel(Metrics)
    metric{k} = Metrics(k).name;
    metricLabel{k} = Metrics(k).label;
    x = Fish.(Metrics(k).name);
    g = Fish.morph;
    valid = isfinite(x) & ismember(g, P.morphOrder);
    x = x(valid);
    g = g(valid);
    counts = zeros(1,numel(P.morphOrder));
    for m = 1:numel(P.morphOrder)
        counts(m) = sum(strcmp(g, P.morphOrder{m}));
    end
    nTotal(k) = numel(x);
    nSurface(k) = counts(strcmp(P.morphOrder, 'Surface'));
    nMolino(k) = counts(strcmp(P.morphOrder, 'Molino'));
    nPachon(k) = counts(strcmp(P.morphOrder, 'Pachon'));

    % Deliberately no three-group omnibus test: it would include the
    % Molino-Pachon separation, which is outside the requested contrasts.
end

pAdjusted = pRaw;
correctionFamily = repmat({''},numel(Metrics),1);
Omnibus = table(metric, metricLabel, analysisScope, nTotal, nSurface, ...
    nMolino, nPachon, H, df, pRaw, pAdjusted, correctionFamily);

% Planned pairwise contrasts. Effect is always groupB - groupA.
nRows = numel(Metrics) * size(P.stats.pairwiseComparisons,1);
metric = cell(nRows,1);
metricLabel = cell(nRows,1);
analysisScope = repmat({scope}, nRows,1);
groupA = cell(nRows,1);
groupB = cell(nRows,1);
nA = zeros(nRows,1);
nB = zeros(nRows,1);
medianA = nan(nRows,1);
medianB = nan(nRows,1);
medianDifference_BminusA = nan(nRows,1);
ciLow = nan(nRows,1);
ciHigh = nan(nRows,1);
cliffsDelta_BvsA = nan(nRows,1);
pRaw = nan(nRows,1);
pAdjusted = nan(nRows,1);
test = repmat({'two-group Kruskal-Wallis'},nRows,1);

row = 0;
for k = 1:numel(Metrics)
    rowsThisMetric = zeros(size(P.stats.pairwiseComparisons,1),1);
    for j = 1:size(P.stats.pairwiseComparisons,1)
        row = row + 1;
        rowsThisMetric(j) = row;
        A = P.stats.pairwiseComparisons{j,1};
        B = P.stats.pairwiseComparisons{j,2};
        xa = Fish.(Metrics(k).name)(strcmp(Fish.morph, A));
        xb = Fish.(Metrics(k).name)(strcmp(Fish.morph, B));
        xa = xa(isfinite(xa));
        xb = xb(isfinite(xb));

        metric{row} = Metrics(k).name;
        metricLabel{row} = Metrics(k).label;
        groupA{row} = A;
        groupB{row} = B;
        nA(row) = numel(xa);
        nB(row) = numel(xb);
        if ~isempty(xa); medianA(row) = median(xa); end
        if ~isempty(xb); medianB(row) = median(xb); end

        if nA(row) >= P.stats.minFishPerMorph && ...
                nB(row) >= P.stats.minFishPerMorph
            medianDifference_BminusA(row) = medianB(row) - medianA(row);
            ci = bootstrapMedianDifferenceCI( ...
                xa, xb, P.stats.nBootstrap, P.stats.bootstrapCI);
            ciLow(row) = ci(1);
            ciHigh(row) = ci(2);
            cliffsDelta_BvsA(row) = cliffsDelta(xb, xa);
            kwGroup = [repmat({A},nA(row),1); repmat({B},nB(row),1)];
            pRaw(row) = kruskalwallis([xa; xb],kwGroup,'off');
        end
    end

    pAdjusted(rowsThisMetric) = pRaw(rowsThisMetric);
end

correctionFamily = repmat({''},nRows,1);
Pairwise = table(metric, metricLabel, analysisScope, test, groupA, groupB, ...
    nA, nB, medianA, medianB, medianDifference_BminusA, ciLow, ciHigh, ...
    cliffsDelta_BvsA, pRaw, pAdjusted, correctionFamily);
end

function [Omnibus,Pairwise] = applyConfiguredCorrections( ...
        Omnibus,Pairwise,Families,P)
% Each outcome contains the three planned pairwise morph contrasts.
% The default method is 'none', so they remain uncorrected. If 'bh' is
% requested as a sensitivity analysis, it is applied within one metric and
% analysis scope; biologically different outcomes are never pooled.
%
% The legacy Omnibus table is retained for compatibility but its p-values are
% NaN because three-group tests are intentionally disabled.
if ~isstruct(Families); error('Families must be a struct array.'); end
scopes = unique(Omnibus.analysisScope,'stable');
for s = 1:numel(scopes)
    for f = 1:numel(Families)
        idxO = strcmp(Omnibus.analysisScope,scopes{s}) & ...
            ismember(Omnibus.metric,Families(f).metrics);

        Omnibus.pAdjusted(idxO) = NaN;
        Omnibus.correctionFamily(idxO) = ...
            repmat({'disabled_three_group_test'},sum(idxO),1);

        % Apply the configured method to the three planned comparisons for
        % this metric. With the default 'none', pAdjusted equals pRaw.
        familyMetrics = Families(f).metrics;
        for k = 1:numel(familyMetrics)
            metricName = familyMetrics{k};
            idxP = strcmp(Pairwise.analysisScope,scopes{s}) & ...
                strcmp(Pairwise.metric,metricName);
            Pairwise.pAdjusted(idxP) = adjustPValues( ...
                Pairwise.pRaw(idxP),P.stats.multipleComparisonMethod);
            Pairwise.correctionFamily(idxP) = repmat( ...
                {[Families(f).name '_within_' metricName]},sum(idxP),1);
        end
    end
end
end

function pAdjusted = adjustPValues(pRaw,method)
switch lower(method)
    case 'bh'
        pAdjusted = bhAdjust(pRaw);
    case 'none'
        pAdjusted = pRaw;
    otherwise
        error('Unknown multiple-comparison method: %s',method);
end
end

function [PairStats,FriedmanStats] = runWithinMorphModelStats(Fish,P)
% Compare model scores within the same fish. Tests are performed separately
% in each morph and once across all fish as a descriptive overall check.
scopes = [P.morphOrder, {'AllMorphs'}];
PairStats=runWithinMorphPairedStats(Fish,P.stats.withinModelComparisons,P);

scope = scopes(:); nFish = zeros(numel(scopes),1); pFriedman = nan(numel(scopes),1);
chiSquare = nan(numel(scopes),1); kendallW = nan(numel(scopes),1);
for s=1:numel(scopes)
    if strcmp(scopes{s},'AllMorphs')
        idx=true(height(Fish),1);
    else
        idx=strcmp(Fish.morph,scopes{s});
    end
    X = nan(height(Fish),numel(P.modelMetricNames));
    for j = 1:numel(P.modelMetricNames)
        X(:,j) = Fish.(P.modelMetricNames{j});
    end
    X=X(idx,:); X=X(all(isfinite(X),2),:); nFish(s)=size(X,1);
    if size(X,1)>=P.stats.minFishPerMorph
        [pFriedman(s),tbl]=friedman(X,1,'off');
        chiSquare(s)=extractFriedmanChiSquare(tbl);
        if isfinite(chiSquare(s))
            kendallW(s)=chiSquare(s)/(size(X,1)*(size(X,2)-1));
        end
    end
end
pAdjusted = pFriedman;
FriedmanStats=table(scope,nFish,chiSquare,pFriedman,pAdjusted,kendallW);
end

function PairStats = runWithinMorphPairedStats(Fish,comparisons,P)
scopes = [P.morphOrder, {'AllMorphs'}];
nRows = numel(scopes)*size(comparisons,1);

scope = cell(nRows,1); metricA = cell(nRows,1); metricB = cell(nRows,1);
modelA = cell(nRows,1); modelB = cell(nRows,1); nFish = zeros(nRows,1);
medianA = nan(nRows,1); medianB = nan(nRows,1);
medianPairedDifference_BminusA = nan(nRows,1);
ciLow = nan(nRows,1); ciHigh = nan(nRows,1);
rankBiserial_BvsA = nan(nRows,1); pRaw = nan(nRows,1);

row = 0;
for s = 1:numel(scopes)
    if strcmp(scopes{s},'AllMorphs')
        idxScope = true(height(Fish),1);
    else
        idxScope = strcmp(Fish.morph,scopes{s});
    end
    for j = 1:size(comparisons,1)
        row = row+1;
        aName = comparisons{j,1}; bName = comparisons{j,2};
        a = Fish.(aName); b = Fish.(bName);
        valid = idxScope & isfinite(a) & isfinite(b);
        a = a(valid); b = b(valid); d = b-a;

        scope{row}=scopes{s}; metricA{row}=aName; metricB{row}=bName;
        modelA{row}=comparisons{j,3}; modelB{row}=comparisons{j,4};
        nFish(row)=numel(d);
        if ~isempty(d)
            medianA(row)=median(a); medianB(row)=median(b);
            medianPairedDifference_BminusA(row)=median(d);
        end
        if numel(d)>=P.stats.minFishPerMorph
            ci=bootstrapPairedMedianDifferenceCI(d,P.stats.nBootstrap,P.stats.bootstrapCI);
            ciLow(row)=ci(1); ciHigh(row)=ci(2);
            rankBiserial_BvsA(row)=pairedRankBiserial(d);
            pRaw(row)=safeSignrank(a,b);
        end
    end
end
% The comparisons form one prespecified family within each morph/scope.
pAdjusted = pRaw;
PairStats=table(scope,metricA,metricB,modelA,modelB,nFish,medianA,medianB, ...
    medianPairedDifference_BminusA,ciLow,ciHigh,rankBiserial_BvsA,pRaw,pAdjusted);
end

function chi2 = extractFriedmanChiSquare(tbl)
chi2=NaN;
try
    value=tbl{2,5};
    if iscell(value); value=value{1}; end
    if isnumeric(value) && isscalar(value); chi2=value; end
catch
end
end

function p = safeSignrank(a,b)
p=NaN;
d=b-a; d=d(isfinite(d));
if isempty(d); return; end
if all(abs(d)<=sqrt(eps)); p=1; return; end
try
    p=signrank(a,b);
catch
    p=NaN;
end
end

function r = pairedRankBiserial(d)
d=d(isfinite(d) & abs(d)>sqrt(eps));
if isempty(d); r=0; return; end
[~,ord]=sort(abs(d)); ranks=zeros(size(d));
sortedAbs=abs(d(ord)); i=1;
while i<=numel(d)
    j=i;
    while j<numel(d) && sortedAbs(j+1)==sortedAbs(i); j=j+1; end
    ranks(ord(i:j))=mean(i:j); i=j+1;
end
wPlus=sum(ranks(d>0)); wMinus=sum(ranks(d<0));
r=(wPlus-wMinus)/(wPlus+wMinus);
end

function ci = bootstrapPairedMedianDifferenceCI(d,nBoot,limits)
d=d(isfinite(d));
if isempty(d); ci=[NaN NaN]; return; end
if numel(d)==1; ci=[d d]; return; end
n=numel(d);
boot=median(d(randi(n,n,nBoot)),1);
ci=prctile(boot(:),limits);
end

function H = extractKWChiSquare(tbl)
H = NaN;
try
    value = tbl{2,5};
    if iscell(value); value = value{1}; end
    if isnumeric(value) && isscalar(value)
        H = value;
    end
catch
end
end

function [estimate, ci] = bootstrapMedianCI(x, nBoot, limits)
x = x(isfinite(x));
if isempty(x)
    estimate = NaN;
    ci = [NaN NaN];
    return;
end
estimate = median(x);
if numel(x) < 2
    ci = [estimate estimate];
    return;
end
n = numel(x);
boot = median(x(randi(n,n,nBoot)),1);
ci = prctile(boot(:), limits);
end

function ci = bootstrapMedianDifferenceCI(a, b, nBoot, limits)
% Difference is median(b) - median(a).
na = numel(a);
nb = numel(b);
bootA = median(a(randi(na,na,nBoot)),1);
bootB = median(b(randi(nb,nb,nBoot)),1);
ci = prctile(bootB(:)-bootA(:), limits);
end

function d = cliffsDelta(x, y)
% Probability(x > y) - probability(x < y). Here x is group B.
if isempty(x) || isempty(y)
    d = NaN;
    return;
end
score = 0;
for i = 1:numel(x)
    score = score + sum(x(i) > y) - sum(x(i) < y);
end
d = score ./ (numel(x) * numel(y));
end

function q = bhAdjust(p)
q = nan(size(p));
valid = find(isfinite(p));
if isempty(valid); return; end
[ps, order] = sort(p(valid));
m = numel(ps);
qs = ps .* m ./ (1:m)';
qs = flipud(cummin(flipud(qs)));
qs = min(qs,1);
q(valid(order)) = qs;
end

function Stats = runPositiveB1Tests(Fish,P)
morph=string(P.morphOrder(:)); nFish=zeros(numel(morph),1); medianB1=nan(numel(morph),1); ciLow=nan(numel(morph),1); ciHigh=nan(numel(morph),1); pGreaterThanZero=nan(numel(morph),1);
for m=1:numel(morph)
 x=Fish.median_b1_raw(strcmp(Fish.morph,morph(m))); x=x(isfinite(x)); nFish(m)=numel(x);
 if isempty(x), continue; end
 medianB1(m)=median(x); [~,ci]=bootstrapMedianCI(x,P.stats.nBootstrap,P.stats.bootstrapCI); ciLow(m)=ci(1); ciHigh(m)=ci(2);
 if numel(x)>=P.stats.minFishPerMorph, pGreaterThanZero(m)=signrank(x,0,"tail","right"); end
end
Stats=table(morph,nFish,medianB1,ciLow,ciHigh,pGreaterThanZero);
end

function [Examples,FishAHV,AHVSamples] = buildRepresentativeExamplesAndAHVSummary(All,P)
SData=struct("morph",{},"session",{},"time",{},"phase",{},"ahv",{},"candidatePath",{},"rawMapping",{});
for m=1:numel(P.inputs)
 filePath=locateResultsFile(P.inputs(m).folder,P.inputs(m).fileName); L=load(filePath,"Results"); if ~isfield(L,"Results"), continue; end
 R=L.Results; if isstruct(R), R=arrayfun(@(x){x},R); end
 for r=1:numel(R)
  Q=R{r}; if ~isstruct(Q)||~isfield(Q,"name")||~isfield(Q,"tCa")||~isfield(Q,"phaseCa")||~isfield(Q,"ahvCa"), continue; end
  mapping=struct(); if isfield(Q,"candidateDiagnostics")&&isfield(Q.candidateDiagnostics,"rawTraceMapping"), mapping=Q.candidateDiagnostics.rawTraceMapping; end; SData(end+1)=struct("morph",string(P.inputs(m).morph),"session",string(Q.name),"time",double(Q.tCa(:)),"phase",double(Q.phaseCa(:)),"ahv",double(Q.ahvCa(:)),"candidatePath",string(Q.paths.candidatePath),"rawMapping",mapping); %#ok<AGROW>
 end
end
Examples=struct("morph",{},"session",{},"animalID",{},"neuronID",{},"cvR2",{},"phaseDominance",{},"windowStartSeconds",{},"time",{},"observed",{},"prediction",{},"phaseTerm",{},"ahvTerm",{},"absAHVTerm",{});
for m=1:numel(P.morphOrder)
 rows=find(strcmp(All.morph,P.morphOrder{m})&All.baseQualityOK&isfinite(All.cvR2Original)&isfinite(All.phaseDominanceR2)); if isempty(rows), continue; end
 [~,o1]=sort(All.cvR2Original(rows)); rankR2=zeros(numel(rows),1); rankR2(o1)=(1:numel(rows))/numel(rows);
 [~,o2]=sort(All.phaseDominanceR2(rows)); rankDom=zeros(numel(rows),1); rankDom(o2)=(1:numel(rows))/numel(rows); score=(rankR2+rankDom)/2;
 targetScores=quantile(score,P.examples.scoreQuantiles); available=true(numel(rows),1); order=zeros(0,1);
 for q=1:numel(targetScores), candidate=find(available); if isempty(candidate), break; end; [~,nearest]=min(abs(score(candidate)-targetScores(q))); chosen=candidate(nearest); order(end+1,1)=chosen; available(chosen)=false; end
 [~,fallbackOrder]=sort(score,"descend"); order=[order;fallbackOrder(available(fallbackOrder))];
 nChosen=0;
 for j=transpose(order)
  row=rows(j); dataIndex=find([SData.morph]==string(P.morphOrder{m})&[SData.session]==string(All.session(row)),1); if isempty(dataIndex), continue; end
  observed=traceCellValue(All.cvObserved,row); prediction=traceCellValue(All.cvPredFull,row); n=min([numel(observed),numel(prediction),numel(SData(dataIndex).time),numel(SData(dataIndex).phase),numel(SData(dataIndex).ahv)]); if n<2, continue; end
  D=SData(dataIndex); responseOK=responseFiniteMask(D,All.neuronID(row)); time=D.time(responseOK); phase=D.phase(responseOK); ahv=D.ahv(responseOK); n=min([numel(observed),numel(prediction),numel(time),numel(phase),numel(ahv)]); assert(n==numel(observed)&&n==numel(time),"OOF trace and finite-response regressor lengths do not match for %s.",D.session); time=time(1:n); phase=phase(1:n); ahv=ahv(1:n); observed=observed(1:n); prediction=prediction(1:n);
  phaseTerm=All.originalA(row).*cos(All.originalPrefRad(row)-phase); ahvTerm=All.originalB2(row).*ahv; absAHVTerm=All.originalB1(row).*abs(ahv);
  window=selectRepresentativeWindow(time,observed,prediction,phaseTerm,ahvTerm,absAHVTerm,P.examples.windowSeconds,nChosen+1); if ~any(window), continue; end
  nChosen=nChosen+1; windowStart=time(find(window,1))-time(1); Examples(end+1)=struct("morph",string(P.morphOrder{m}),"session",string(All.session(row)),"animalID",string(All.animalID(row)),"neuronID",All.neuronID(row),"cvR2",All.cvR2Original(row),"phaseDominance",All.phaseDominanceR2(row),"windowStartSeconds",windowStart,"time",time(window),"observed",observed(window),"prediction",prediction(window),"phaseTerm",phaseTerm(window),"ahvTerm",ahvTerm(window),"absAHVTerm",absAHVTerm(window)); %#ok<AGROW>
  if nChosen>=P.examples.numberPerMorph, break; end
 end
end
fishKeys=unique(string(All.fishKey),"stable"); morph=strings(numel(fishKeys),1); animalID=strings(numel(fishKeys),1); meanAbsAHV=nan(numel(fishKeys),1); sdSignedAHV=nan(numel(fishKeys),1); nSamples=zeros(numel(fishKeys),1);
AHVSamples=struct("morph",cell(numel(fishKeys),1),"fishKey",cell(numel(fishKeys),1),"values",cell(numel(fishKeys),1));
for f=1:numel(fishKeys)
 q=string(All.fishKey)==fishKeys(f); first=find(q,1); morph(f)=string(All.morph(first)); animalID(f)=string(All.animalID(first)); sessions=unique(string(All.session(q))); x=[];
 for s=1:numel(sessions), k=find([SData.morph]==morph(f)&[SData.session]==sessions(s),1); if ~isempty(k), x=[x;SData(k).ahv(isfinite(SData(k).ahv))]; end; end %#ok<AGROW>
 nSamples(f)=numel(x); if ~isempty(x), meanAbsAHV(f)=mean(abs(x)); sdSignedAHV(f)=std(x,0); end
 AHVSamples(f).morph=morph(f); AHVSamples(f).fishKey=fishKeys(f); AHVSamples(f).values=x;
end
FishAHV=table(morph,animalID,fishKeys,meanAbsAHV,sdSignedAHV,nSamples);
end

function window = selectRepresentativeWindow(time,varargin)
windowSeconds=varargin{end-1}; windowRank=varargin{end}; signals=varargin(1:end-2);
window=false(size(time)); finiteTime=find(isfinite(time)); if isempty(finiteTime), return; end
starts=time(finiteTime(1)):windowSeconds:time(finiteTime(end)); scores=nan(size(starts)); masks=cell(size(starts));
for k=1:numel(starts)
 masks{k}=time>=starts(k)&time<starts(k)+windowSeconds;
 if nnz(masks{k})<2, continue; end
 z=0; for s=1:numel(signals), v=signals{s}(masks{k}); v=v(isfinite(v)); if ~isempty(v), z=z+std(v,0); end; end; scores(k)=z;
end
[~,order]=sort(scores,"descend","MissingPlacement","last"); order=order(isfinite(scores(order))); if isempty(order), return; end
window=masks{order(1+mod(windowRank-1,min(3,numel(order))))};
end

function ok = responseFiniteMask(D,neuronID)
C=load(D.candidatePath,"calcium_traces","candidate_cell_ids","time_s"); ids=double(C.candidate_cell_ids(:)); candidatePosition=find(ids==double(neuronID),1); assert(~isempty(candidatePosition),"Representative neuron ID is absent from its candidate file.");
if ~isempty(fieldnames(D.rawMapping))
    M=D.rawMapping; originalID=M.originalRoiIds(candidatePosition); R=load(M.rasterFile,"deltaFoF"); Y=double(R.deltaFoF);
    if size(Y,2)>=originalID, y=Y(M.rawFrameStart:M.rawFrameEnd,originalID); elseif size(Y,1)>=originalID, y=transpose(Y(originalID,M.rawFrameStart:M.rawFrameEnd)); else, error("RASTER dimensions do not contain mapped ROI."); end
else
    Y=double(C.calcium_traces); if size(Y,2)==numel(ids), y=Y(:,candidatePosition); elseif size(Y,1)==numel(ids), y=transpose(Y(candidatePosition,:)); else, error("Candidate trace dimensions do not match IDs."); end
end
[found,location]=ismember(D.time,double(C.time_s(:))); assert(all(found),"Saved fit times do not map exactly to candidate times."); ok=isfinite(y(location));
end

function x = traceCellValue(column,row)
if iscell(column), x=double(column{row}(:)); else, x=double(column(row,:)); x=x(:); end
end

function Stats = runAHVRegressorStats(T,Samples,P)
metrics = ["meanAbsAHV";"sdSignedAHV"];
metricLabels = ["Mean |AHV| per fish";"SD of signed AHV per fish"];
distributionComparisons = {"Surface","Molino";"Surface","Pachon"};
metricOut=strings(0,1); metricLabel=strings(0,1); test=strings(0,1);
groupA=strings(0,1); groupB=strings(0,1); nA=zeros(0,1); nB=zeros(0,1);
statistic=nan(0,1); omnibusP=nan(0,1); pRaw=nan(0,1); pAdjusted=nan(0,1);
for k=1:numel(metrics)
    values=[]; groups=strings(0,1); counts=zeros(numel(P.morphOrder),1);
    for m=1:numel(P.morphOrder)
        x=T.(metrics(k))(T.morph==string(P.morphOrder{m})); x=x(isfinite(x));
        counts(m)=numel(x); values=[values;x]; %#ok<AGROW>
        groups=[groups;repmat(string(P.morphOrder{m}),numel(x),1)]; %#ok<AGROW>
    end
    if all(counts>=P.stats.minFishPerMorph)
        [pKW,~,kwStats]=kruskalwallis(values,groups,"off");
        [posthoc,~,~,groupNames]=multcompare(kwStats, ...
            "CType","dunn-sidak","Display","off");
        for j=1:size(posthoc,1)
            ia=posthoc(j,1); ib=posthoc(j,2);
            metricOut(end+1,1)=metrics(k); metricLabel(end+1,1)=metricLabels(k); %#ok<AGROW>
            test(end+1,1)="Dunn-Sidak post hoc after one-way Kruskal-Wallis"; %#ok<AGROW>
            groupA(end+1,1)=string(groupNames{ia}); groupB(end+1,1)=string(groupNames{ib}); %#ok<AGROW>
            nA(end+1,1)=counts(strcmp(P.morphOrder,groupNames{ia})); %#ok<AGROW>
            nB(end+1,1)=counts(strcmp(P.morphOrder,groupNames{ib})); %#ok<AGROW>
            statistic(end+1,1)=posthoc(j,4); omnibusP(end+1,1)=pKW; %#ok<AGROW>
            pRaw(end+1,1)=NaN; pAdjusted(end+1,1)=posthoc(j,6); %#ok<AGROW>
        end
    end
end
[D,~,~]=normalizedFishAHVDensities(Samples,P);
for j=1:size(distributionComparisons,1)
    a=distributionComparisons{j,1}; b=distributionComparisons{j,2};
    ia=find([Samples.morph]==a); ib=find([Samples.morph]==b);
    ia=ia(all(isfinite(D(ia,:)),2)); ib=ib(all(isfinite(D(ib,:)),2));
    [distance,pv]=fishLabelPermutationCDFTest(D(ia,:),D(ib,:), ...
        P.ahvDistribution.nPermutations);
    metricOut(end+1,1)="signedAHVDistribution"; %#ok<AGROW>
    metricLabel(end+1,1)="Equal-fish-weight signed-AHV distribution"; %#ok<AGROW>
    test(end+1,1)="fish-label permutation of mean-CDF sup distance"; %#ok<AGROW>
    groupA(end+1,1)=a; groupB(end+1,1)=b; nA(end+1,1)=numel(ia); nB(end+1,1)=numel(ib); %#ok<AGROW>
    statistic(end+1,1)=distance; omnibusP(end+1,1)=NaN; %#ok<AGROW>
    pRaw(end+1,1)=pv; pAdjusted(end+1,1)=pv; %#ok<AGROW>
end
groupA=cellstr(groupA); groupB=cellstr(groupB);
Stats=table(metricOut,metricLabel,test,groupA,groupB,nA,nB,statistic,omnibusP,pRaw,pAdjusted);
end

function [D,centers,edges] = normalizedFishAHVDensities(Samples,P)
allQuantiles=nan(numel(Samples),2);
for f=1:numel(Samples)
    x=Samples(f).values; x=x(isfinite(x));
    if ~isempty(x), allQuantiles(f,:)=prctile(x,P.ahvDistribution.tailPercentiles); end
end
lo=min(allQuantiles(:,1),[],"omitnan"); hi=max(allQuantiles(:,2),[],"omitnan");
limit=max(abs([lo hi]));
if ~isfinite(limit)||limit<=0, limit=1; end
edges=linspace(-limit,limit,P.ahvDistribution.nBins+1);
centers=(edges(1:end-1)+edges(2:end))/2;
D=nan(numel(Samples),P.ahvDistribution.nBins);
for f=1:numel(Samples)
    x=Samples(f).values; x=x(isfinite(x)); counts=histcounts(x,edges);
    if sum(counts)>0, D(f,:)=counts./sum(counts)./diff(edges); end
end
end

function [observed,p] = fishLabelPermutationCDFTest(A,B,nPermutations)
if isempty(A)||isempty(B), observed=NaN; p=NaN; return; end
cdfA=cumsum(A./sum(A,2),2); cdfB=cumsum(B./sum(B,2),2);
observed=max(abs(mean(cdfA,1)-mean(cdfB,1)));
pooled=[cdfA;cdfB]; nA=size(A,1); n=size(pooled,1); exceed=0;
for r=1:nPermutations
    order=randperm(n); delta=max(abs(mean(pooled(order(1:nA),:),1)- ...
        mean(pooled(order(nA+1:end),:),1)));
    exceed=exceed+(delta>=observed);
end
p=(exceed+1)/(nPermutations+1);
end

function plotRepresentativeExamples(E,P,outputDir,baseName)
if isempty(E), warning("No compatible OOF traces were available for representative examples."); return; end
rowFields=["observed","phaseTerm","ahvTerm","absAHVTerm"]; rowLabels=["Activity / OOF prediction","Phase term","AHV term","|AHV| term"];
for exampleNumber=1:P.examples.numberPerMorph
 fig=figure("Visible",P.figure.visible,"Color","w","Position",[30 30 1450 620]);
 for m=1:3
  morphExamples=find([E.morph]==string(P.morphOrder{m})); e=[]; if numel(morphExamples)>=exampleNumber, e=morphExamples(exampleNumber); end
  for r=1:4
   ax=subplot(4,3,(r-1)*3+m,"Parent",fig); hold(ax,"on"); if isempty(e), title(ax,P.morphOrder{m}+" (no compatible example "+exampleNumber+")"); continue; end
   keep=1:max(1,ceil(numel(E(e).time)/1500)):numel(E(e).time); t=(E(e).time(keep)-E(e).time(1))/60;
   if r==1
    plot(ax,t,E(e).observed(keep),"Color",[.25 .25 .25],"LineWidth",.8,"DisplayName","Observed"); plot(ax,t,E(e).prediction(keep),"Color",P.morphColors(m,:),"LineWidth",1.4,"DisplayName","OOF full prediction");
    title(ax,sprintf("%s | %s | neuron %g | start %.1f min | CV R2 %.3f | dominance %.3f",P.morphOrder{m},E(e).session,E(e).neuronID,E(e).windowStartSeconds/60,E(e).cvR2,E(e).phaseDominance),"Interpreter","none");
   else
    values=E(e).(rowFields(r)); plot(ax,t,values(keep),"Color",P.morphColors(m,:),"LineWidth",1);
   end
   if m==1, ylabel(ax,rowLabels(r)); end; if r==4, xlabel(ax,"Time (min)"); else, set(ax,"XTickLabel",[]); end; grid(ax,"on"); box(ax,"off");
  end
 end
 assert(isgraphics(fig,"figure"),"Representative figure was unexpectedly closed before export."); drawnow; assert(isgraphics(fig,"figure"),"Representative figure closed during drawnow."); saveRepresentativeFigure(fig,outputDir,sprintf('%s_%d',baseName,exampleNumber),P);
end
end

function saveRepresentativeFigure(fig,outputDir,baseName,P)
if P.figure.savePNG
 print(fig,fullfile(outputDir,[baseName '.png']),'-dpng',sprintf('-r%d',P.figure.dpi));
 print(fig,fullfile(outputDir,[baseName '.svg']),'-dsvg');
end
if strcmpi(P.figure.visible,"off"), close(fig); end
end

function plotAHVRegressorDistributions(T,Samples,Stats,P,outputDir,baseName)
fig=figure("Visible",P.figure.visible,"Color","w", ...
    "Position",[60 60 1080 400],"Renderer","painters");
tl=tiledlayout(fig,1,3,"TileSpacing","compact","Padding","compact");
vars=["meanAbsAHV","sdSignedAHV"];
labels=["Mean |AHV| per fish","SD of signed AHV per fish"];
for k=1:numel(vars)
    ax=nexttile(tl); hold(ax,"on"); tickLabels=cell(size(P.morphOrder));
    for m=1:numel(P.morphOrder)
        y=T.(vars(k))(T.morph==string(P.morphOrder{m})); y=y(isfinite(y));
        tickLabels{m}=sprintf("%s (n=%d fish)",P.morphOrder{m},numel(y));
        if isempty(y), continue; end
        [med,ci]=bootstrapMedianCI(y,P.stats.nBootstrap,P.stats.bootstrapCI);
        bar(ax,m,med,P.figure.morphBarWidth,"FaceColor",P.morphColors(m,:), ...
            "FaceAlpha",.78,"EdgeColor","k","LineWidth",.8,"HandleVisibility","off");
        errorbar(ax,m,med,med-ci(1),ci(2)-med,"k","LineStyle","none", ...
            "LineWidth",1.7,"CapSize",8,"HandleVisibility","off");
        scatter(ax,m+deterministicJitter(numel(y),P.figure.jitterWidth),y, ...
            max(18,round(P.figure.pointSize*.55)), ...
            "MarkerFaceColor",mixWithWhite(P.morphColors(m,:),.38), ...
            "MarkerEdgeColor","k","LineWidth",.45,"MarkerFaceAlpha",.82, ...
            "HandleVisibility","off");
    end
    set(ax,"XTick",1:numel(P.morphOrder),"XTickLabel",tickLabels, ...
        "TickDir","out","Box","off");
    xlim(ax,P.figure.morphXLimits); ylabel(ax,labels(k)); grid(ax,"on");
    sub=Stats(Stats.metricOut==vars(k),:);
    title(ax,plannedKWTitle(sub),"FontWeight","normal");
    addMorphPairwiseBars(ax,sub,P.morphOrder);
end
ax=nexttile(tl); hold(ax,"on");
[D,centers,~]=normalizedFishAHVDensities(Samples,P);
for m=1:numel(P.morphOrder)
    rows=find([Samples.morph]==string(P.morphOrder{m}) & all(isfinite(D),2)');
    pale=mixWithWhite(P.morphColors(m,:),.72);
    for f=rows
        plot(ax,centers,D(f,:),"Color",pale,"LineWidth",.45, ...
            "HandleVisibility","off");
    end
    if ~isempty(rows)
        plot(ax,centers,mean(D(rows,:),1),"Color",P.morphColors(m,:), ...
            "LineWidth",2.2,"DisplayName",sprintf("%s (n=%d fish)", ...
            P.morphOrder{m},numel(rows)));
    end
end
xline(ax,0,":","Color",[.35 .35 .35],"HandleVisibility","off");
xlabel(ax,"Signed AHV"); ylabel(ax,"Probability density (equal fish weight)");
set(ax,"TickDir","out","Box","off"); grid(ax,"on"); legend(ax,"Location","best");
distStats=Stats(Stats.metricOut=="signedAHVDistribution",:);
title(ax,"Per-fish-normalized signed-AHV distributions", ...
    "FontWeight","normal");
addDistributionTestText(ax,distStats);
title(tl,"Fish-level AHV summaries", ...
    "Interpreter","none","FontWeight","bold");
saveFigureBoth(fig,outputDir,baseName,P);
end

function addDistributionTestText(ax,T)
pieces=strings(height(T),1);
for i=1:height(T)
    pieces(i)=string(shortMorph(T.groupA{i}))+"-"+ ...
        string(shortMorph(T.groupB{i}))+": "+ ...
        string(significanceLabel(T.pAdjusted(i)))+ ...
        " (p="+string(formatP(T.pAdjusted(i)))+")";
end
text(ax,.02,.98,strjoin(cellstr(pieces),newline),"Units","normalized", ...
    "HorizontalAlignment","left","VerticalAlignment","top","FontSize",8, ...
    "Interpreter","none","BackgroundColor","w","Margin",2);
end

function Summary = plotPairedIncrementalR2(Fish,Stats,P,outputDir,baseName)
labels = {'Blocked-CV Delta R^2 phase','Blocked-CV Delta R^2 AHV'};
vars = {'median_deltaR2Phase','median_deltaR2Behavior'};
offsets = [-0.12 0 0.12];
fig = figure('Visible',P.figure.visible,'Color','w','Position',[80 80 720 500]);
ax = axes(fig); hold(ax,'on'); yline(ax,0,':','Color',[0.6 0.6 0.6], ...
    'HandleVisibility','off');
Summary = table(); pText = cell(numel(P.morphOrder),1);
for m = 1:numel(P.morphOrder)
    idx = strcmp(Fish.morph,P.morphOrder{m});
    X = [Fish.(vars{1}),Fish.(vars{2})];
    X = X(idx,:); X = X(all(isfinite(X),2),:);
    means = mean(X,1); ci = nan(2,2);
    for j = 1:2
        ci(j,:) = bootstrapMeanCI(X(:,j),P.stats.nBootstrap,P.stats.bootstrapCI);
    end
    positions = [1 2] + offsets(m);
    plot(ax,positions,means,'-o','Color',P.morphColors(m,:), ...
        'MarkerFaceColor',P.morphColors(m,:),'MarkerEdgeColor','k', ...
        'LineWidth',1.8,'DisplayName',P.morphOrder{m});
    for j = 1:2
        errorbar(ax,positions(j),means(j),means(j)-ci(j,1),ci(j,2)-means(j), ...
            'Color',P.morphColors(m,:),'LineStyle','none','LineWidth',1.5, ...
            'CapSize',8,'HandleVisibility','off');
        morph = string(P.morphOrder{m}); metric = string(vars{j});
        nFish = size(X,1); meanFishMedians = means(j); ciLow = ci(j,1); ciHigh = ci(j,2);
        Summary = appendCompatibleTables(Summary, ...
            table(morph,metric,nFish,meanFishMedians,ciLow,ciHigh));
    end
    S = Stats(strcmp(Stats.scope,P.morphOrder{m}),:);
    if isempty(S); p = NaN; else; p = S.pAdjusted(1); end
    pText{m} = sprintf('%s paired p=%s',P.morphOrder{m},formatP(p));
end
set(ax,'XTick',1:2,'XTickLabel',labels,'XTickLabelRotation',0, ...
    'TickDir','out','Box','off'); xlim(ax,[0.6 2.4]); grid(ax,'on');
ylabel(ax,'Mean fish-median blocked-CV Delta R^2');
legend(ax,'Location','best');
text(ax,0.02,0.98,strjoin(pText,newline),'Units','normalized', ...
    'VerticalAlignment','top','FontSize',9,'BackgroundColor','w','Margin',2);
title(ax,['Cross-validated incremental predictive contributions by morph' newline ...
    'Centers: mean across fish medians; error bars: fish-bootstrap 95% CI'], ...
    'FontWeight','normal');
saveFigureBoth(fig,outputDir,baseName,P);
end

function ci = bootstrapMeanCI(x,nBoot,limits)
x = x(isfinite(x));
if isempty(x); ci = [NaN NaN]; return; end
if numel(x)==1; ci = [x x]; return; end
n = numel(x); boot = mean(x(randi(n,n,nBoot)),1);
ci = prctile(boot(:),limits);
end

function plotPairedModelCV(Fish,PairStats,FriedmanStats,P,outputDir,baseName)
figureModelOrder = [1 3 2]; % phase + AHV, phase, AHV
models = P.modelNames(figureModelOrder);
vars = P.modelMetricNames(figureModelOrder);
allModelColors=[0.55 0.55 0.55; 0.45 0.20 0.70; ...
    0.10 0.62 0.56; 0.25 0.55 0.85];
modelColors = allModelColors(figureModelOrder,:);
fig=figure('Visible',P.figure.visible,'Color','w', ...
    'Position',[50 80 470*numel(P.morphOrder) 430]);
tl=tiledlayout(fig,1,numel(P.morphOrder),'TileSpacing','compact','Padding','compact');
for m=1:numel(P.morphOrder)
    ax=nexttile(tl); hold(ax,'on'); yline(ax,0,':','Color',[0.65 0.65 0.65]);
    idx=strcmp(Fish.morph,P.morphOrder{m});
    X=[Fish.(vars{1}),Fish.(vars{2}),Fish.(vars{3})];
    X=X(idx,:); X=X(all(isfinite(X),2),:);
    for f=1:size(X,1)
        plot(ax,1:numel(models),X(f,:),'-o','Color',mixWithWhite(P.morphColors(m,:),0.58), ...
            'MarkerFaceColor',mixWithWhite(P.morphColors(m,:),0.35), ...
            'MarkerEdgeColor','none','LineWidth',0.9,'HandleVisibility','off');
    end
    for j=1:numel(models)
        [med,ci]=bootstrapMedianCI(X(:,j),P.stats.nBootstrap,P.stats.bootstrapCI);
        plot(ax,[j j],ci,'k-','LineWidth',2.2,'HandleVisibility','off');
        scatter(ax,j,med,80,'d','MarkerFaceColor',modelColors(j,:), ...
            'MarkerEdgeColor','k','LineWidth',1.1,'HandleVisibility','off');
    end
    set(ax,'XTick',1:numel(models),'XTickLabel',models,'XTickLabelRotation',18, ...
        'TickDir','out','Box','off'); xlim(ax,[0.65 numel(models)+0.35]); grid(ax,'on');
    ylabel(ax,'Median held-out R^2 per fish');
    fr=FriedmanStats(strcmp(FriedmanStats.scope,P.morphOrder{m}),:);
    if isempty(fr); pF=NaN; else; pF=fr.pAdjusted(1); end
    title(ax,sprintf('%s | Friedman p=%s',P.morphOrder{m},formatP(pF)), ...
        'FontWeight','normal');
    sub=PairStats(strcmp(PairStats.scope,P.morphOrder{m}),:);
    lines=cell(height(sub),1);
    for j=1:height(sub)
        lines{j}=sprintf('%s - %s: p=%s',sub.modelB{j},sub.modelA{j}, ...
            formatP(sub.pAdjusted(j)));
    end
    if ~isempty(lines)
        text(ax,0.02,0.98,strjoin(lines,'\n'),'Units','normalized', ...
            'VerticalAlignment','top','FontSize',8,'BackgroundColor','w','Margin',2);
    end
end
title(tl,['Paired fish-level model comparison | identical neurons, frames and ' ...
    sprintf('%d contiguous folds',P.cv.expectedBlockedFolds)], ...
    'Interpreter','none','FontWeight','bold');
saveFigureBoth(fig,outputDir,baseName,P);
end

function plotOriginalVsAugmentedCoefficients(Fish,P,outputDir,baseName)
origVars={'median_originalB0','median_abs_originalB1','median_abs_originalB2'};
fullVars={'median_b0_aug','median_abs_b1_aug','median_abs_b2_aug'};
labels={'Median b_0','Median |b_1|','Median |b_2|'};
nRows=numel(P.morphOrder); nCols=numel(labels);
fig=figure('Visible',P.figure.visible,'Color','w', ...
    'Position',[30 30 390*nCols 300*nRows]);
tl=tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
for m=1:nRows
    idxMorph=strcmp(Fish.morph,P.morphOrder{m});
    for k=1:nCols
        ax=nexttile(tl); hold(ax,'on');
        a=Fish.(origVars{k})(idxMorph); b=Fish.(fullVars{k})(idxMorph);
        valid=isfinite(a)&isfinite(b); a=a(valid); b=b(valid);
        for f=1:numel(a)
            plot(ax,[1 2],[a(f) b(f)],'-o', ...
                'Color',mixWithWhite(P.morphColors(m,:),0.52), ...
                'MarkerFaceColor',P.morphColors(m,:),'MarkerEdgeColor','none', ...
                'LineWidth',0.9,'HandleVisibility','off');
        end
        if ~isempty(a)
            plot(ax,1,median(a),'kd','MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.3);
            plot(ax,2,median(b),'kd','MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.3);
        end
        set(ax,'XTick',[1 2], ...
            'XTickLabel',{'phase + AHV','phase + AHV'}, ...
            'TickDir','out','Box','off'); xlim(ax,[0.65 2.35]); grid(ax,'on');
        if k==1; ylabel(ax,P.morphOrder{m},'FontWeight','bold'); end
        p=safeSignrank(a,b);
        title(ax,sprintf('%s | paired p=%s',labels{k},formatP(p)), ...
            'FontWeight','normal');
    end
end
title(tl,['Matched coefficient sensitivity: same neurons and frames ' ...
    '(effect of adding forward + vigor to phase + AHV)'], ...
    'Interpreter','none','FontWeight','bold');
saveFigureBoth(fig,outputDir,baseName,P);
end

function plotFishLevelMetrics(Fish, Metrics, metricNames, Omnibus, Pairwise, ...
        P, outputDir, figureTitle, baseName)
selected = selectMetrics(Metrics, metricNames);
n = numel(selected);
if n == 4
    nCols = 2;
else
    nCols = min(3,n);
end
nRows = ceil(n/nCols);
fig = figure('Visible', P.figure.visible, 'Color', 'w', ...
    'Position', [60 60 360*nCols 360*nRows]);
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact','Padding','compact');

for k = 1:n
    ax = nexttile(tl); hold(ax,'on');
    plotOneFishMetric(ax, Fish, selected(k), Omnibus, Pairwise, P);
    if strcmp(baseName, '06_primary_raw_coefficients_fish_medians')
        ax.XTickLabelRotation = 0;
    end
end
title(tl, figureTitle, 'Interpreter','none', 'FontWeight','bold');
saveFigureBoth(fig, outputDir, baseName, P);
end

function plotOneFishMetric(ax, Fish, M, Omnibus, Pairwise, P)
if isfinite(M.referenceLine)
    yline(ax, M.referenceLine, ':', 'Color',[0.65 0.65 0.65], 'LineWidth',1);
end

for m = 1:numel(P.morphOrder)
    morph = P.morphOrder{m};
    y = Fish.(M.name)(strcmp(Fish.morph, morph));
    y = y(isfinite(y));
    if isempty(y); continue; end

    [med, ci] = bootstrapMedianCI(y, P.stats.nBootstrap, P.stats.bootstrapCI);
    bar(ax,m,med,P.figure.morphBarWidth,'FaceColor',P.morphColors(m,:), ...
        'FaceAlpha',0.78,'EdgeColor','k','LineWidth',0.8, ...
        'HandleVisibility','off');
    errorbar(ax,m,med,med-ci(1),ci(2)-med,'k','LineStyle','none', ...
        'LineWidth',1.7,'CapSize',8,'HandleVisibility','off');

    jitter = deterministicJitter(numel(y), P.figure.jitterWidth);
    scatter(ax, m + jitter, y, max(18,round(P.figure.pointSize*0.55)), ...
        'MarkerFaceColor', mixWithWhite(P.morphColors(m,:),0.38), ...
        'MarkerEdgeColor', 'k', 'LineWidth',0.45, ...
        'MarkerFaceAlpha',0.82);
end

xlim(ax,P.figure.morphXLimits);
tickLabels = cell(size(P.morphOrder));
for m = 1:numel(P.morphOrder)
    nFishShown = sum(strcmp(Fish.morph,P.morphOrder{m}) & isfinite(Fish.(M.name)));
    tickLabels{m} = sprintf('%s (n=%d fish)',P.morphOrder{m},nFishShown);
end
set(ax, 'XTick',1:numel(P.morphOrder), 'XTickLabel',tickLabels, ...
    'TickDir','out', 'Box','off');
ylabel(ax, M.label, 'Interpreter','tex');
grid(ax,'on');
if strcmp(M.summary,'fraction_positive') || startsWith(M.name,'fraction_')
    ylim(ax,[0 1]);
end

idxP = strcmp(Pairwise.metric, M.name);
sub = Pairwise(idxP,:);
title(ax, plannedKWTitle(sub), 'FontWeight','normal');
addMorphPairwiseBars(ax,sub,P.morphOrder);
if strcmp(M.name,"median_b1_raw")
    pieces=strings(numel(P.morphOrder),1);
    for m=1:numel(P.morphOrder)
        x=Fish.(M.name)(strcmp(Fish.morph,P.morphOrder{m})); x=x(isfinite(x)); pv=NaN;
        if numel(x)>=P.stats.minFishPerMorph, pv=signrank(x,0,"tail","right"); end
        pieces(m)=string(P.morphOrder{m})+": b_1>0 p="+string(formatP(pv));
    end
    text(ax,.02,.98,strjoin(cellstr(pieces),newline),"Units","normalized", ...
        "VerticalAlignment","top","BackgroundColor","w","Margin",2,"FontSize",8);
end
end

function jitter = deterministicJitter(n, width)
if n <= 1
    jitter = 0;
else
    jitter = linspace(-width, width, n)';
end
end

function plotPairwiseEstimation(Pairwise, Metrics, metricNames, P, ...
        outputDir, figureTitle, baseName)
selected = selectMetrics(Metrics, metricNames);
n = numel(selected);
nCols = min(3,n);
nRows = ceil(n/nCols);
fig = figure('Visible', P.figure.visible, 'Color','w', ...
    'Position',[70 70 440*nCols 340*nRows]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');

for k = 1:n
    ax = nexttile(tl); hold(ax,'on');
    yline(ax,0,':','Color',[0.55 0.55 0.55]);
    sub = Pairwise(strcmp(Pairwise.metric,selected(k).name),:);
    for j = 1:height(sub)
        colorIdx = find(strcmp(P.morphOrder,sub.groupB{j}),1);
        if isempty(colorIdx); colorIdx = 1; end
        plot(ax,[j j],[sub.ciLow(j) sub.ciHigh(j)],'-', ...
            'Color',P.morphColors(colorIdx,:),'LineWidth',2.5);
        scatter(ax,j,sub.medianDifference_BminusA(j),65, ...
            'MarkerFaceColor',P.morphColors(colorIdx,:), ...
            'MarkerEdgeColor','k');
        text(ax,j,sub.ciHigh(j),sprintf('  p=%s',formatP(sub.pAdjusted(j))), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',8);
    end
    labels = cell(height(sub),1);
    for j = 1:height(sub)
        labels{j} = sprintf('%s-%s',shortMorph(sub.groupB{j}),shortMorph(sub.groupA{j}));
    end
    xlim(ax,[0.5 max(1.5,height(sub)+0.5)]);
    set(ax,'XTick',1:height(sub),'XTickLabel',labels,'TickDir','out','Box','off');
    ylabel(ax,['Difference in ' selected(k).label], 'Interpreter','tex');
    title(ax,selected(k).label,'FontWeight','normal','Interpreter','tex');
    grid(ax,'on');
end
title(tl,figureTitle,'Interpreter','none','FontWeight','bold');
saveFigureBoth(fig,outputDir,baseName,P);
end

function plotEqualFishECDFs(All, Metrics, metricNames, P, outputDir, ...
        figureTitle, baseName)
selected = selectMetrics(Metrics, metricNames);
n = numel(selected);
nCols = min(3,n);
nRows = ceil(n/nCols);
fig = figure('Visible',P.figure.visible,'Color','w', ...
    'Position',[60 60 455*nCols 350*nRows]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');

for k = 1:n
    ax = nexttile(tl); hold(ax,'on');
    M = selected(k);

    pooled = extractNeuronValues(All,M,All.baseQualityOK);
    if isempty(pooled); continue; end
    lim = prctile(pooled,P.ecdf.percentileLimits);
    if lim(1) == lim(2)
        pad = max(1,abs(lim(1))) * 0.05;
        lim = lim + [-pad pad];
    end
    gridX = linspace(lim(1),lim(2),P.ecdf.nGridPoints);

    for m = 1:numel(P.morphOrder)
        fish = unique(All.fishKey(strcmp(All.morph,P.morphOrder{m})),'stable');
        F = nan(numel(fish),numel(gridX));
        for f = 1:numel(fish)
            mask = strcmp(All.fishKey,fish{f}) & All.baseQualityOK;
            x = extractNeuronValues(All,M,mask);
            if numel(x) < P.quality.minNeuronsPerFish; continue; end
            for q = 1:numel(gridX)
                F(f,q) = mean(x <= gridX(q));
            end
            if P.ecdf.showIndividualFish
                plot(ax,gridX,F(f,:),'-','Color', ...
                    mixWithWhite(P.morphColors(m,:),0.72),'LineWidth',0.6);
            end
        end

        meanF = mean(F,1,'omitnan');
        nF = sum(isfinite(F),1);
        semF = std(F,0,1,'omitnan') ./ sqrt(max(nF,1));
        valid = isfinite(meanF);
        if any(valid)
            patch(ax,[gridX(valid) fliplr(gridX(valid))], ...
                [meanF(valid)-semF(valid) fliplr(meanF(valid)+semF(valid))], ...
                P.morphColors(m,:), 'FaceAlpha',0.13,'EdgeColor','none');
            plot(ax,gridX(valid),meanF(valid),'-','Color',P.morphColors(m,:), ...
                'LineWidth',2.5,'DisplayName',P.morphOrder{m});
        end
    end

    if isfinite(M.referenceLine) && ~strcmp(M.summary,'fraction_positive')
        xline(ax,M.referenceLine,':','Color',[0.55 0.55 0.55]);
    end
    xlim(ax,lim); ylim(ax,[0 1]);
    xlabel(ax,neuronAxisLabel(M),'Interpreter','tex');
    ylabel(ax,'Cumulative fraction of neurons');
    title(ax,M.label,'FontWeight','normal','Interpreter','tex');
    grid(ax,'on'); set(ax,'TickDir','out','Box','off');
    if k == 1
        legend(ax,'Location','best');
    end
end
title(tl,figureTitle,'Interpreter','none','FontWeight','bold');
saveFigureBoth(fig,outputDir,baseName,P);
end

function label = neuronAxisLabel(M)
switch M.transform
    case 'abs'
        label = ['Neuron-level |' strrep(M.source,'_','\_') '|'];
    otherwise
        label = ['Neuron-level ' strrep(M.source,'_','\_')];
end
end

function color = mixWithWhite(color, fractionWhite)
color = color .* (1-fractionWhite) + [1 1 1] .* fractionWhite;
end

function plotClassStratifiedMetrics(FishClass, Metrics, Omnibus, Pairwise, ...
        P, outputDir, baseName)
nRows = numel(P.classAnalysis.classes);
nCols = numel(Metrics);
fig = figure('Visible',P.figure.visible,'Color','w', ...
    'Position',[30 30 330*nCols 305*nRows]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');

for c = 1:nRows
    cls = P.classAnalysis.classes{c};
    F = FishClass(strcmp(FishClass.class,cls),:);
    for k = 1:nCols
        ax = nexttile(tl); hold(ax,'on');
        M = Metrics(k);
        if isfinite(M.referenceLine)
            yline(ax,M.referenceLine,':','Color',[0.65 0.65 0.65]);
        end
        for m = 1:numel(P.morphOrder)
            y = F.(M.name)(strcmp(F.morph,P.morphOrder{m}));
            y = y(isfinite(y));
            if isempty(y); continue; end
            [med,ci] = bootstrapMedianCI(y,P.stats.nBootstrap,P.stats.bootstrapCI);
            bar(ax,m,med,P.figure.morphBarWidth,'FaceColor',P.morphColors(m,:), ...
                'FaceAlpha',0.78,'EdgeColor','k','LineWidth',0.7, ...
                'HandleVisibility','off');
            errorbar(ax,m,med,med-ci(1),ci(2)-med,'k','LineStyle','none', ...
                'LineWidth',1.5,'CapSize',7,'HandleVisibility','off');
            scatter(ax,m+deterministicJitter(numel(y),P.figure.jitterWidth),y, ...
                max(16,round(P.figure.pointSize*0.48)), ...
                'MarkerFaceColor',mixWithWhite(P.morphColors(m,:),0.38), ...
                'MarkerEdgeColor','k','LineWidth',0.4,'MarkerFaceAlpha',0.82);
        end
        set(ax,'XTick',1:numel(P.morphOrder),'XTickLabel',P.morphOrder, ...
            'TickDir','out','Box','off');
        xlim(ax,P.figure.morphXLimits); grid(ax,'on');
        if c < nRows; set(ax,'XTickLabel',[]); end
        if k == 1; ylabel(ax,cls,'FontWeight','bold'); end
        % Display all three pre-specified pairwise two-group KW tests.
        idxP = strcmp(Pairwise.analysisScope,cls) & ...
            strcmp(Pairwise.metric,M.name);
        sub = Pairwise(idxP,:);
        title(ax,sprintf('%s | %s',M.label,plannedKWTitle(sub)), ...
            'FontWeight','normal','Interpreter','tex','FontSize',9);
        addMorphPairwiseBars(ax,sub,P.morphOrder);
    end
end
title(tl,['Class-stratified fish summaries (descriptive; class is model-derived)' newline ...
    sprintf('Minimum %d neurons per fish/class',P.classAnalysis.minNeuronsPerFishClass)], ...
    'Interpreter','none','FontWeight','bold');
saveFigureBoth(fig,outputDir,baseName,P);
end

function addMorphPairwiseBars(ax,sub,morphOrder)
% Add the three planned pairwise morph comparisons when
% their fish-level two-group Kruskal-Wallis p-values are available.
if isempty(sub) || ~ismember('pAdjusted',sub.Properties.VariableNames)
    return;
end
sub = sub(isfinite(sub.pAdjusted),:);
if isempty(sub)
    return;
end

limits = ylim(ax);
if strcmpi(ax.YScale,'log')
    if limits(2) <= 0; return; end
    factor = 1.32;
    dataTop = limits(2);
    for j = 1:height(sub)
        x1 = find(strcmpi(morphOrder,sub.groupA{j}),1);
        x2 = find(strcmpi(morphOrder,sub.groupB{j}),1);
        if isempty(x1) || isempty(x2); continue; end
        y = dataTop .* factor.^j;
        yTick = y ./ 1.045;
        plot(ax,[x1 x1 x2 x2],[yTick y y yTick],'k-', ...
            'LineWidth',1.15,'HandleVisibility','off', ...
            'Clipping','off');
        text(ax,mean([x1 x2]),y.*1.025, ...
            significanceLabel(sub.pAdjusted(j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontWeight','bold');
    end
    ylim(ax,[limits(1) dataTop.*factor.^(height(sub)+0.55)]);
else
    dataRange = limits(2)-limits(1);
    if ~isfinite(dataRange) || dataRange <= 0
        dataRange = max(1,abs(limits(2)));
    end
    base = limits(2) + 0.08.*dataRange;
    step = 0.13.*dataRange;
    tick = 0.035.*dataRange;
    for j = 1:height(sub)
        x1 = find(strcmpi(morphOrder,sub.groupA{j}),1);
        x2 = find(strcmpi(morphOrder,sub.groupB{j}),1);
        if isempty(x1) || isempty(x2); continue; end
        y = base + (j-1).*step;
        plot(ax,[x1 x1 x2 x2],[y-tick y y y-tick],'k-', ...
            'LineWidth',1.15,'HandleVisibility','off', ...
            'Clipping','off');
        text(ax,mean([x1 x2]),y+0.02.*dataRange, ...
            significanceLabel(sub.pAdjusted(j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontWeight','bold');
    end
    ylim(ax,[limits(1) base + max(height(sub)-1,0).*step + 0.16.*dataRange]);
end
addMorphPValueText(ax,sub);
end

function addMorphPValueText(ax,sub)
pieces = strings(height(sub),1);
for i = 1:height(sub)
    pieces(i) = string(sub.groupA{i}) + "-" + string(sub.groupB{i}) + ...
        ": p=" + string(sprintf('%.3g',sub.pAdjusted(i)));
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

function selected = selectMetrics(Metrics, names)
selected = Metrics([]);
for i = 1:numel(names)
    idx = find(strcmp({Metrics.name},names{i}),1);
    assert(~isempty(idx),'Unknown metric name: %s',names{i});
    selected(end+1) = Metrics(idx); %#ok<AGROW>
end
end

function s = shortMorph(morph)
switch morph
    case 'Surface'; s = 'Surf';
    case 'Molino';  s = 'Mol';
    case 'Pachon';  s = 'Pach';
    otherwise;      s = morph;
end
end

function s = plannedKWTitle(sub)
if isempty(sub)
    s = 'Two-group KW: NA';
    return;
end
parts = cell(height(sub),1);
for i = 1:height(sub)
    parts{i} = sprintf('%s-%s p=%s',shortMorph(sub.groupA{i}), ...
        shortMorph(sub.groupB{i}),formatP(sub.pAdjusted(i)));
end
s = ['Two-group KW: ' strjoin(parts,'; ')];
end

function s = formatP(p)
if ~isfinite(p)
    s = 'NA';
elseif p < 1e-4
    s = sprintf('%.1e',p);
else
    s = sprintf('%.4f',p);
end
end

function saveFigureBoth(fig, outputDir, baseName, P)
if P.figure.savePNG
    file = fullfile(outputDir,[baseName '.png']);
    try
        exportgraphics(fig,file,'Resolution',P.figure.dpi);
    catch
        print(fig,file,'-dpng',sprintf('-r%d',P.figure.dpi));
    end
    svgFile = fullfile(outputDir,[baseName '.svg']);
    try
        exportgraphics(fig,svgFile,'ContentType','vector');
    catch
        print(fig,svgFile,'-dsvg');
    end
end
if strcmpi(P.figure.visible,'off')
    close(fig);
end
end

function writeAnalysisNotes(outputDir, P)
file = fullfile(outputDir,'README_analysis_notes.txt');
fid = fopen(file,'w');
if fid < 0
    warning('Could not write %s',file);
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'Cross-morph matched three-model HD/AHV analysis\n\n');
fprintf(fid,'Source fitting script: %s\n',P.sourceFittingScript);
fprintf(fid,'Independent replicate: fish/session.\n');
fprintf(fid,'Neurons are nested within fish and are not treated as independent n.\n');
fprintf(fid,'CV folds are temporal partitions and are not treated as independent n.\n');
fprintf(fid,'Expected CV scheme: %d contiguous blocked folds.\n',P.cv.expectedBlockedFolds);
if P.exclusion.enabled && ~isempty(P.exclusion.recNumbers)
    fprintf(fid,'Recording exclusions enabled: rec numbers %s.\n', ...
        mat2str(P.exclusion.recNumbers));
    fprintf(fid,'Exclusion reason: %s.\n',P.exclusion.reason);
    fprintf(fid,['These exclusions are applied before all summaries, plots, ROC ' ...
        'analyses and tests; the ExclusionReport table is retained in the MAT workspace.\n']);
else
    fprintf(fid,'Recording exclusions disabled.\n');
end
fprintf(fid,'Model 1 - phase + AHV: phase cosine + phase sine + absolute AHV + signed AHV.\n');
fprintf(fid,'Model 2 - AHV: absolute AHV + signed AHV.\n');
fprintf(fid,'Model 3 - phase: phase cosine + phase sine.\n');
fprintf(fid,'No three-group omnibus test is run.\n');
fprintf(fid,'Planned tests: two-group Kruskal-Wallis for Surface-Molino, Surface-Pachon, and Molino-Pachon fish medians.\n');
fprintf(fid,'Within-morph model tests are paired sign-rank tests on fish summaries.\n');
fprintf(fid,'cvUniquePhase = 1 - SSE_(phase+AHV) / SSE_AHV.\n');
fprintf(fid,'cvUniqueAHVBeyondPhase = 1 - SSE_(phase+AHV) / SSE_phase.\n');
fprintf(fid,'All three models use identical neurons, frames and CV folds.\n');
fprintf(fid,'Figure 10 centers are means across fish-level neuron medians; error bars are fish-bootstrap 95%% confidence intervals.\n');
fprintf(fid,'Figure 10 paired tests retain paired fish-level values and never test displayed morph averages.\n');
fprintf(fid,['Absolute phase-model performance is conditional on selecting ' ...
    'phase-tuned neurons upstream; nested comparisons use the same selected set.\n']);
fprintf(fid,'VIF and condition number are diagnostics unless exclusion is explicitly enabled.\n');
fprintf(fid,'Class-stratified results are secondary/descriptive because class was model-derived.\n');
if P.fishR2Exclusion.enabled
    fprintf(fid,[ ...
        'Fish-level model-performance exclusion enabled: a fish was removed ' ...
        'from all figures and tests if at least one of the three model median ' ...
        'held-out R2 values was below %.3f.\n'], ...
        P.fishR2Exclusion.threshold);
else
    fprintf(fid,'Fish-level model-performance exclusion disabled.\n');
end
end

function x = wrapToPiLocal(x)
x = mod(x+pi,2*pi)-pi;
end
