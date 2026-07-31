function results = compare_HD_AHV_model_with_phase_only_behavior_params_across_morphs_v3(cfg)
%COMPARE_HD_AHV_MODEL_WITH_PHASE_ONLY_BEHAVIOR_PARAMS_ACROSS_MORPHS_V3
% Compare four explicitly named neural-activity models across morphs.
%
% This script reads the morph-level output produced by:
%   fit_neuron_HD_AHV_four_models_behavior_tuned_cells_only_ROC.m
%
% It expects one result file per morph:
%   HD_AHV_behavior_augmented_phase_tuned_only_with_phase_only_results.mat
%
% FOUR MODELS ARE COMPARED ON THE SAME NEURONS, FRAMES AND CV FOLDS
% ----------------------------------------------------------------
%   phase + AHV:                   cos(phi)+sin(phi)+|AHV|+AHV
%   phase + AHV + forward + vigor: cos(phi)+sin(phi)+|AHV|+AHV+F+V
%   phase:                         cos(phi)+sin(phi)
%   AHV + forward + vigor:         |AHV|+AHV+F+V
%
% The saved variable names cvR2Original and cvR2Full are retained by the
% fitting script for backward compatibility. In every user-facing label,
% "Original" is called "phase + AHV" and "Full" is called
% "phase + AHV + forward + vigor".
%
% The fitting script uses five contiguous blocked folds. This comparison
% script checks the saved fold count and reports aggregate held-out R^2.
% Folds are not treated as independent replicates.
%
% PRIMARY INFERENTIAL UNIT
% ------------------------
% Fish/session is the independent biological replicate. Neurons are nested
% within fish and are therefore NOT entered as independent observations in
% cross-morph statistical tests. Each requested outcome is first summarized
% within fish.
%
% PRIMARY QUESTIONS
% -----------------
%   1. Do phase + AHV + forward + vigor coefficients differ by morph?
%   2. Does adding forward + vigor alter the phase + AHV coefficients?
%   3. How do the four explicitly named models compare in held-out R^2?
%   4. Does AHV add held-out information beyond phase?
%   5. Is phase predictive beyond all measured behavior?
%   6. Do all measured behavior predictors add information beyond phase?
%   7. Are phase-amplitude estimates stable after adding AHV and controls?
%
% IMPORTANT INTERPRETATION
% ------------------------
% * Because neural activity was z-scored, b0 is in activity-SD units.
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
% * Omnibus three-morph Kruskal-Wallis test for every plotted fish-level metric.
% * Two pre-specified pairwise Wilcoxon rank-sum comparisons by default:
%       Molino - Surface, Pachon - Surface.
% * Pairwise effect size: difference between fish-level medians, with a
%   fish-level bootstrap 95% CI, plus Cliff's delta.
% * Each Kruskal-Wallis test is the single omnibus test for its outcome and
%   is reported without pooling it with biologically different outcomes.
% * The two planned pairwise contrasts are reported without multiplicity
%   correction by default, as explicitly configured below.
% * Within each morph, model scores are paired by fish and compared with
%   sign-rank tests. Cross-morph tests on matched unique-contribution
%   summaries provide a non-parametric model-by-morph sensitivity analysis.
%
% OUTPUTS
% -------
%   cross_morph_HD_AHV_phase_only_behavior_augmented_analysis.mat
%   README_analysis_notes.txt
%   figures as PNG and SVG (no CSV or PDF files)

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
    'fit_neuron_HD_AHV_four_models_behavior_tuned_cells_only_ROC.m';
P.modelNames = { ...
    'phase + AHV', ...
    'phase + AHV + forward + vigor', ...
    'phase', ...
    'AHV + forward + vigor'};
P.modelMetricNames = { ...
    'median_cvR2Original', ...
    'median_cvR2Full', ...
    'median_cvR2PhaseOnly', ...
    'median_cvR2Behavior'};

% Input folders supplied by the user. Folder order does not determine plot
% order; P.morphOrder below does.
P.inputs = struct( ...
    'morph', {'Surface','Molino','Pachon'}, ...
    'folder', {cfg.ModelDataDir, cfg.ModelDataDir, cfg.ModelDataDir});

P.resultsFileName = ...
    'HD_AHV_behavior_augmented_phase_tuned_only_with_phase_only_results.mat';
for iInput = 1:numel(P.inputs)
    P.inputs(iInput).folder = cfg.ModelDataDir;
    P.inputs(iInput).fileName = sprintf([ ...
        'HD_AHV_behavior_augmented_phase_tuned_only_' ...
        'with_phase_only_results_%s.mat'], ...
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
P.quality.expectedDesignRank = 7;
P.quality.maxConditionNumber = 1e8;
P.quality.minNeuronsPerFish = 5;
P.quality.excludeByMaxVIF = false; % diagnostic by default, not an exclusion
P.quality.maxVIF = 10;

% Require the saved fitting run to use five contiguous blocked folds. The
% comparison script reads already-computed CV outputs; it does not refit.
P.cv.expectedBlockedFolds = 5;
P.cv.requireExpectedFoldCount = true;

% Statistics.
P.stats.alpha = 0.05;
P.stats.minFishPerMorph = 3;
P.stats.nBootstrap = 10000;
P.stats.bootstrapCI = [2.5 97.5];
P.stats.randomSeed = 7;
P.stats.pairwiseComparisons = { ...
    'Surface','Molino'; ...
    'Surface','Pachon'};
% The two contrasts are explicit planned comparisons and are left
% uncorrected. Set this to 'bh' only if a corrected sensitivity analysis is
% desired; Holm correction is intentionally not implemented in this script.
P.stats.multipleComparisonMethod = 'none';

% A Kruskal-Wallis test is already the single omnibus comparison across all
% three morphs for one outcome. Leave these p-values unadjusted by default.
% Set this true only if you explicitly want to correct omnibus tests across
% the different metrics within each figure family.
P.stats.correctOmnibusAcrossMetrics = false;
P.stats.withinModelComparisons = { ...
    'median_cvR2Original','median_cvR2Full', ...
        'phase + AHV','phase + AHV + forward + vigor'; ...
    'median_cvR2Original','median_cvR2PhaseOnly', ...
        'phase + AHV','phase'; ...
    'median_cvR2Original','median_cvR2Behavior', ...
        'phase + AHV','AHV + forward + vigor'; ...
    'median_cvR2Full','median_cvR2PhaseOnly', ...
        'phase + AHV + forward + vigor','phase'; ...
    'median_cvR2Full','median_cvR2Behavior', ...
        'phase + AHV + forward + vigor','AHV + forward + vigor'; ...
    'median_cvR2PhaseOnly','median_cvR2Behavior', ...
        'phase','AHV + forward + vigor'};
P.stats.withinCoefficientComparisons = { ...
    'median_phaseOnlyB0_std','median_b0_std', ...
        'phase: standardized b0', ...
        'phase + AHV + forward + vigor: standardized b0'; ...
    'median_phaseOnlyB0','median_originalB0', ...
        'phase: b0','phase + AHV: b0'; ...
    'median_phaseOnlyB0','median_b0_aug', ...
        'phase: b0','phase + AHV + forward + vigor: b0'; ...
    'median_originalB0','median_b0_aug', ...
        'phase + AHV: b0','phase + AHV + forward + vigor: b0'; ...
    'median_originalB1','median_b1_aug', ...
        'phase + AHV: signed b1', ...
        'phase + AHV + forward + vigor: signed b1'; ...
    'median_abs_originalB1','median_abs_b1_aug', ...
        'phase + AHV: |b1|','phase + AHV + forward + vigor: |b1|'; ...
    'median_originalB2','median_b2_aug', ...
        'phase + AHV: signed b2', ...
        'phase + AHV + forward + vigor: signed b2'; ...
    'median_abs_originalB2','median_abs_b2_aug', ...
        'phase + AHV: |b2|','phase + AHV + forward + vigor: |b2|'};

% Secondary class-stratified description. This is useful for asking whether
% an overall morph effect is due to class composition, but interpret it
% cautiously because class membership was derived from these coefficients.
P.classAnalysis.enabled = true;
P.classAnalysis.classes = {'Symmetric','CW','CCW'};
P.classAnalysis.minNeuronsPerFishClass = 5;

% Optional ROC analysis. This is only scientifically defined when the
% fitting output contains held-out (out-of-fold) observed activity and the
% four matched held-out prediction traces for every neuron. The binary
% target is observed z-scored activity > P.roc.activityThresholdZ. This is a
% secondary activity-state discrimination analysis; it does not replace R^2.
% Recommended table variable names are cvObserved, cvPredPhaseOnly,
% cvPredOriginal, cvPredBehavior and cvPredFull; each row must contain a
% numeric trace.
P.roc.enabled = true;
P.roc.requirePredictionTraces = false;
P.roc.activityThresholdZ = 0;
P.roc.fprGrid = linspace(0,1,101);
P.roc.minPositiveFrames = 25;
P.roc.minNegativeFrames = 25;
P.roc.minNeuronsPerFish = 3;
P.roc.observedVariableCandidates = { ...
    'cvObserved','cvYObserved','yObservedCV','yTrueCV','yTrueOOF', ...
    'cvYTrue','cvObservedActivity'};
P.roc.phaseOnlyPredictionCandidates = { ...
    'cvPredPhaseOnly','cvPredictionPhaseOnly','yHatPhaseOnlyCV', ...
    'yPredPhaseOnlyOOF','oofPredPhaseOnly'};
P.roc.originalPredictionCandidates = { ...
    'cvPredOriginal','cvPredictionOriginal','yHatOriginalCV', ...
    'yPredOriginalOOF','oofPredOriginal'};
P.roc.behaviorPredictionCandidates = { ...
    'cvPredBehavior','cvPredictionBehavior','yHatBehaviorCV', ...
    'yPredBehaviorOOF','oofPredBehavior'};
P.roc.fullPredictionCandidates = { ...
    'cvPredFull','cvPredictionFull','yHatFullCV', ...
    'yPredFullOOF','oofPredFull'};

% Figure export.
P.figure.visible = cfg.FigureVisible;       % 'on' for interactive inspection
P.figure.savePNG = true;
P.figure.dpi = 300;
P.figure.pointSize = 46;
P.figure.jitterWidth = 0.15;

rng(P.stats.randomSeed, 'twister');

%% ======================== METRIC DEFINITIONS ========================
% transform: identity | abs
% summary:   median | fraction_positive | mean
Metrics = [ ...
    ... % phase model
    makeMetric('median_phaseOnlyB0',       'phaseOnlyB0', 'identity','median', ...
        'Median b_0: phase', true, NaN), ...
    makeMetric('median_phaseOnlyB0_std',   'phaseOnlyB0_std','identity','median', ...
        'Median standardized b_0: phase', true, NaN), ...
    ... % phase + AHV + forward + vigor coefficients
    makeMetric('median_b0_aug',            'b0',          'identity','median', ...
        'Median b_0: phase + AHV + forward + vigor', true, NaN), ...
    makeMetric('median_b1_aug',            'b1',          'identity','median', ...
        'Median signed b_1: phase + AHV + forward + vigor', true, 0), ...
    makeMetric('median_abs_b1_aug',        'b1',          'abs','median', ...
        'Median |b_1|: phase + AHV + forward + vigor', true, NaN), ...
    makeMetric('median_b2_aug',            'b2',          'identity','median', ...
        'Median signed b_2: phase + AHV + forward + vigor', true, 0), ...
    makeMetric('median_abs_b2_aug',        'b2',          'abs','median', ...
        'Median |b_2|: phase + AHV + forward + vigor', true, NaN), ...
    makeMetric('fraction_b2_positive_aug', 'b2',          'identity','fraction_positive', ...
        'Fraction b_2 > 0: phase + AHV + forward + vigor', true, 0.5), ...
    ... % Exactly matched four-model five-fold contiguous blocked CV
    makeMetric('median_cvR2PhaseOnly',     'cvR2PhaseOnly','identity','median', ...
        'Median blocked-CV R^2: phase', true, 0), ...
    makeMetric('median_cvR2Original',      'cvR2Original','identity','median', ...
        'Median blocked-CV R^2: phase + AHV', true, 0), ...
    makeMetric('median_cvR2Behavior',      'cvR2Behavior','identity','median', ...
        'Median blocked-CV R^2: AHV + forward + vigor', true, 0), ...
    makeMetric('median_cvR2Full',          'cvR2Full',    'identity','median', ...
        'Median blocked-CV R^2: phase + AHV + forward + vigor', true, 0), ...
    makeMetric('median_cvDeltaR2AddedBehavior','cvDeltaR2AddedBehavior','identity','median', ...
        'Median DeltaCV-R^2: adding forward + vigor to phase + AHV', true, 0), ...
    makeMetric('median_cvUniqueAddedBehavior','cvUniqueAddedBehavior','identity','median', ...
        'Median unique forward + vigor beyond phase + AHV', true, 0), ...
    makeMetric('median_cvDeltaR2Phase',    'cvDeltaR2Phase','identity','median', ...
        'Median DeltaCV-R^2: adding phase to AHV + forward + vigor', true, 0), ...
    makeMetric('median_cvUniquePhase',     'cvUniquePhase','identity','median', ...
        'Median unique phase beyond AHV + forward + vigor', true, 0), ...
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
        'Fraction neurons: phase + AHV + forward + vigor > phase', false, 0.5), ...
    makeMetric('fraction_full_better_original','fullBetterOriginal','identity','mean', ...
        'Fraction neurons: phase + AHV + forward + vigor > phase + AHV', ...
        false, 0.5), ...
    makeMetric('fraction_full_better_behavior','fullBetterBehavior','identity','mean', ...
        'Fraction neurons: phase + AHV + forward + vigor > AHV + forward + vigor', ...
        false, 0.5), ...
    makeMetric('fraction_cvR2PhaseOnly_positive','cvR2PhaseOnly','identity', ...
        'fraction_positive','Fraction phase CV-R^2 > 0', false, 0.5), ...
    makeMetric('fraction_cvR2Original_positive','cvR2Original','identity','fraction_positive', ...
        'Fraction phase + AHV CV-R^2 > 0', false, 0.5), ...
    makeMetric('fraction_cvR2Behavior_positive','cvR2Behavior','identity','fraction_positive', ...
        'Fraction AHV + forward + vigor CV-R^2 > 0', false, 0.5), ...
    makeMetric('fraction_cvR2Full_positive','cvR2Full','identity','fraction_positive', ...
        'Fraction phase + AHV + forward + vigor CV-R^2 > 0', false, 0.5), ...
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
        'Median Delta b_0: phase + AHV + forward + vigor minus phase', false, 0), ...
    makeMetric('median_abs_prefShift_AHVBeyondPhase', ...
        'absPrefShiftOriginalFromPhaseOnlyDeg','identity','median', ...
        'Median |phase shift|: phase + AHV vs phase (deg)', false, 0), ...
    makeMetric('median_abs_prefShift_AllBehaviorBeyondPhase', ...
        'absPrefShiftFullFromPhaseOnlyDeg','identity','median', ...
        'Median |phase shift|: phase + AHV + forward + vigor vs phase (deg)', ...
        false, 0), ...
    ... % Cross-morph standardized phase/AHV sensitivity controls
    makeMetric('median_b0_std',            'b0_std',      'identity','median', ...
        'Median standardized b_0: phase + AHV + forward + vigor', false, NaN), ...
    makeMetric('median_b1_std',            'b1_std',      'identity','median', ...
        'Median standardized b_1: phase + AHV + forward + vigor', false, 0), ...
    makeMetric('median_abs_b1_std',        'b1_std',      'abs','median', ...
        'Median |standardized b_1|', false, NaN), ...
    makeMetric('median_b2_std',            'b2_std',      'identity','median', ...
        'Median standardized b_2: phase + AHV + forward + vigor', false, 0), ...
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
        'Median gain for positive AHV', false, 0), ...
    makeMetric('median_gNegative',         'gNegative','identity','median', ...
        'Median gain for negative AHV magnitude', false, 0), ...
    makeMetric('median_directionalStrength','directionalStrength','identity','median', ...
        'Median directional AHV strength', false, NaN), ...
    makeMetric('median_directionalAsymmetry','directionalAsymmetry','identity','median', ...
        'Median normalized directional asymmetry', false, 0)];

modelR2MetricNames = P.modelMetricNames;
phaseAmplitudeMetricNames = {'median_originalB0','median_b0_aug', ...
    'median_phaseOnlyB0'};
behaviorCoefficientMetricNames = {'median_bF_std','median_abs_bF_std', ...
    'median_bV_std','median_abs_bV_std'};

% Only these outcomes are plotted and included in the revised inferential
% families. Keeping the full Metrics array above preserves all requested CSV
% summaries without generating inferential rows for unused diagnostics.
uniqueVarianceMetricNames = {'median_cvUniquePhase', ...
    'median_cvUniqueAHVBeyondPhase', ...
    'median_cvUniqueAllBehaviorBeyondPhase', ...
    'median_cvUniqueAddedBehavior'};
signedCoefficientShiftMetricNames = {'median_deltaB0_controls', ...
    'median_deltaB1_controls','median_deltaB2_controls'};
signedStandardizedMetricNames = {'median_phaseOnlyB0_std','median_b0_std', ...
    'median_b1_std','median_b2_std'};
reducedGainMetricNames = {'median_gPositive','median_gNegative', ...
    'median_directionalStrength'};
classCoefficientMetricNames = {'median_phaseOnlyB0_std','median_b0_std', ...
    'median_b1_std','median_b2_std'};
classModelPerformanceMetricNames = modelR2MetricNames;
classUniqueVarianceMetricNames = {'median_cvUniquePhase', ...
    'median_cvUniqueAHVBeyondPhase','median_cvUniqueAllBehaviorBeyondPhase'};
classMetricNames = unique([classCoefficientMetricNames, ...
    classModelPerformanceMetricNames,classUniqueVarianceMetricNames],'stable');

StatFamilies = [ ...
    makeFamily('model_cv', modelR2MetricNames), ...
    makeFamily('unique_variance', uniqueVarianceMetricNames), ...
    makeFamily('phase_amplitude', phaseAmplitudeMetricNames), ...
    makeFamily('forward_vigor', behaviorCoefficientMetricNames), ...
    makeFamily('signed_coefficient_change', signedCoefficientShiftMetricNames), ...
    makeFamily('standardized_phase_AHV', signedStandardizedMetricNames), ...
    makeFamily('directional_AHV_gain', reducedGainMetricNames)];
statsMetricNames = unique([StatFamilies.metrics], 'stable');
StatsMetrics = selectMetrics(Metrics, statsMetricNames);

%% ======================== LOAD ALL MORPHS ========================
fprintf('\nLoading matched four-model outputs from %s...\n', ...
    P.sourceFittingScript);
[AllNeurons, InputReport, ExclusionReport] = loadAllMorphResults(P);
assert(~isempty(AllNeurons), 'No neuron tables were loaded. Check P.inputs.');

AllNeurons = addDerivedNeuronVariables(AllNeurons);
requiredVariables = unique([{Metrics.source}, {'session','morph','fishKey'}]);
missing = setdiff(requiredVariables, AllNeurons.Properties.VariableNames);
assert(isempty(missing), 'AllNeurons is missing required variables: %s', ...
    strjoin(missing, ', '));

AllNeurons.baseQualityOK = makeBaseQualityMask(AllNeurons, P);

fprintf('Loaded %d neurons from %d fish.\n', height(AllNeurons), ...
    numel(unique(AllNeurons.fishKey)));
fprintf('Matched finite four-model CV outputs: %d/%d quality-passing neurons.\n', ...
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
ModelDefinitions = buildModelDefinitionTable(P);

% Morph estimates are medians across fish with fish-level bootstrap CIs.
MorphEstimates = computeMorphEstimates(FishMetrics, Metrics, P);

% Inferential tests on fish-level values only and only for prespecified
% metrics that appear in figures.
[OmnibusStats, PairwiseStats] = runCrossMorphStats( ...
    FishMetrics, StatsMetrics, P, 'All selected phase-tuned neurons');
[OmnibusStats, PairwiseStats] = applyConfiguredCorrections( ...
    OmnibusStats, PairwiseStats, StatFamilies, P);

% Paired model tests use one aggregate CV score per fish. They do not treat
% neurons or folds as independent n. Model-difference metrics above are also
% compared across morphs and serve as interaction sensitivity analyses.
[WithinModelStats, ModelFriedmanStats] = runWithinMorphModelStats(FishMetrics, P);
WithinCoefficientStats = runWithinMorphPairedStats( ...
    FishMetrics,P.stats.withinCoefficientComparisons,P);

%% ======================== PRIMARY FIGURES ========================
plotFishLevelMetrics(FishMetrics, Metrics, modelR2MetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    'Matched four-model blocked cross-validation: cross-morph comparison', ...
    '01_cross_morph_four_model_CV_R2');

plotPairedModelCV(FishMetrics,WithinModelStats,ModelFriedmanStats,P, ...
    outputDir,'02_within_morph_paired_four_model_CV_R2');

% ROC requires saved held-out traces; aggregate R^2 values cannot be
% converted into a valid ROC curve. The analysis is skipped with an explicit
% warning when those traces are absent.
FishROC = table();
ROCPairedStats = table();
ROCFriedmanStats = table();
if P.roc.enabled
    [FishROC, rocColumns] = computeFishROCFromOOFTraces(AllNeurons, P);
    if isempty(FishROC)
        msg = ['ROC not generated: no compatible out-of-fold observed and ' ...
            'prediction trace columns were found. Save these traces in the ' ...
            'fitting output; aggregate CV-R^2 is insufficient for ROC/AUC.'];
        if P.roc.requirePredictionTraces; error(msg); else; warning(msg); end
    else
        [ROCPairedStats, ROCFriedmanStats] = runROCModelStats(FishROC,P);
        plotModelROC(FishROC,ROCFriedmanStats,P,outputDir, ...
            '03_four_model_ROC_AUC');
        fprintf(['ROC columns: observed=%s, phase + AHV=%s, ' ...
            'phase + AHV + forward + vigor=%s, phase=%s, ' ...
            'AHV + forward + vigor=%s\n'],rocColumns.observed, ...
            rocColumns.original,rocColumns.full,rocColumns.phaseOnly, ...
            rocColumns.behavior);
    end
end

plotFishLevelMetrics(FishMetrics, Metrics, uniqueVarianceMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    'Matched cross-validated unique phase and behavior contributions', ...
    '04_cross_validated_unique_variance');

plotFishLevelMetrics(FishMetrics, Metrics, phaseAmplitudeMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    'Phase amplitude across nested models', ...
    '05_phase_amplitude_across_nested_models');

plotFishLevelMetrics(FishMetrics, Metrics, behaviorCoefficientMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    'Behavior-control coefficients (standardized for cross-morph comparison)', ...
    '06_standardized_forward_vigor_coefficients');

plotFishLevelMetrics(FishMetrics, Metrics, signedCoefficientShiftMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    'Signed coefficient changes after adding behavior controls', ...
    '07_signed_coefficient_changes_after_behavior_controls');

plotFishLevelMetrics(FishMetrics, Metrics, signedStandardizedMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    ['Standardized coefficients: phase and ' ...
        'phase + AHV + forward + vigor'], ...
    '08_signed_standardized_phase_AHV_coefficients');

plotFishLevelMetrics(FishMetrics, Metrics, reducedGainMetricNames, ...
    OmnibusStats, PairwiseStats, P, outputDir, ...
    'Positive- and negative-AHV gain re-expression', ...
    '09_direction_specific_AHV_gains');

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

        classModelPerformanceMetrics = selectMetrics(Metrics, ...
            classModelPerformanceMetricNames);
        plotClassStratifiedMetrics(FishClassMetrics, ...
            classModelPerformanceMetrics, ...
            ClassOmnibusStats, ClassPairwiseStats, P, outputDir, ...
            '11_class_stratified_four_model_CV_R2');

        classUniqueVarianceMetrics = selectMetrics(Metrics, ...
            classUniqueVarianceMetricNames);
        plotClassStratifiedMetrics(FishClassMetrics, ...
            classUniqueVarianceMetrics, ...
            ClassOmnibusStats, ClassPairwiseStats, P, outputDir, ...
            '12_class_stratified_unique_variance');
    else
        warning('No class variable found. Skipping class-stratified analysis.');
    end
end

%% ======================== SAVE COMPLETE WORKSPACE ========================
save(fullfile(outputDir, ...
    'cross_morph_HD_AHV_phase_only_behavior_augmented_analysis.mat'), ...
    'P', 'Metrics', 'StatsMetrics', 'StatFamilies', ...
    'AllNeurons', 'FishMetrics', 'MorphEstimates', ...
    'OmnibusStats', 'PairwiseStats', 'FishClassMetrics', ...
    'ClassOmnibusStats', 'ClassPairwiseStats', 'WithinModelStats', ...
    'ModelFriedmanStats', 'WithinCoefficientStats', 'FishROC', ...
    'ROCPairedStats','ROCFriedmanStats','InputReport','ModelDefinitions', ...
    'ExclusionReport','-v7.3');

writeAnalysisNotes(outputDir, P);

fprintf('\nDone. Main inferential n is the number of fish, not neurons.\n');
fprintf('Results saved in:\n%s\n', outputDir);

results = struct('outputDir', outputDir, ...
    'fishMetrics', FishMetrics, 'morphEstimates', MorphEstimates, ...
    'omnibusStatistics', OmnibusStats, ...
    'pairwiseStatistics', PairwiseStats, ...
    'inputReport', InputReport, 'exclusionReport', ExclusionReport);
end

%% ======================== LOCAL FUNCTIONS ========================

function M = makeMetric(name, source, transform, summary, label, requested, referenceLine)
M = struct('name', name, 'source', source, 'transform', transform, ...
    'summary', summary, 'label', label, 'requested', requested, ...
    'referenceLine', referenceLine);
end

function T = buildModelDefinitionTable(P)
displayName = P.modelNames(:);
savedCVField = {'cvR2Original';'cvR2Full'; ...
    'cvR2PhaseOnly';'cvR2Behavior'};
savedFitPrefix = {'original';'augmented';'phaseOnly';'behavior'};
formula = { ...
    'cos(phi) + sin(phi) + |AHV| + AHV'; ...
    'cos(phi) + sin(phi) + |AHV| + AHV + forward + vigor'; ...
    'cos(phi) + sin(phi)'; ...
    '|AHV| + AHV + forward + vigor'};
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
    validateFourModelOutputTable(T,filePath,P);

    assert(ismember('session', T.Properties.VariableNames), ...
        'The table in %s has no session variable.', filePath);
    T.session = normalizeTextColumn(T.session);
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
    T.fishKey = strcat(T.morph, {'__'}, T.session);

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
    zscoreActivity, cvScheme, sourceFittingScript, originalModelFormula, ...
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
nonemptyFormula = InputReport.augmentedModelFormula( ...
    ~cellfun(@isempty, InputReport.augmentedModelFormula));
if numel(unique(nonemptyFormula)) > 1
    warning(['Saved phase + AHV + forward + vigor formulas are not ' ...
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
    warning(['Saved AHV + forward + vigor formulas are not identical ' ...
        'across morph files.']);
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

function validateFourModelOutputTable(T,filePath,P)
% These fields distinguish the complete four-model output from earlier
% versions that computed AHV + forward + vigor only inside CV.
required = { ...
    'phaseOnlyR2','phaseOnlyAIC','phaseOnlyBIC', ...
    'originalR2','originalAIC','originalBIC', ...
    'r2','AIC','BIC', ...
    'behaviorR2','behaviorAIC','behaviorBIC', ...
    'behaviorB1','behaviorB2','behaviorBF','behaviorBV', ...
    'cvR2Original','cvR2Full','cvR2PhaseOnly','cvR2Behavior', ...
    'cvSSEOriginal','cvSSEFull','cvSSEPhaseOnly','cvSSEBehavior'};
missing = setdiff(required,T.Properties.VariableNames);
if ~isempty(missing)
    error(['The file %s is not a complete output of %s. Missing fields: %s. ' ...
        'Re-run the four-model fitting script for this morph.'], ...
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
% These are deterministic re-expressions of already saved matched fits.
required = {'b1','b2','originalB1','originalB2','phaseOnlyB0', ...
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
fishKey = cell(nRows,1);
classLabel = cell(nRows,1);
nNeuronsAvailable = zeros(nRows,1);
nNeuronsQuality = zeros(nRows,1);

Fish = table(morph, session, fishKey, classLabel, ...
    nNeuronsAvailable, nNeuronsQuality, ...
    'VariableNames', {'morph','session','fishKey','class', ...
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
    'cvDeltaR2AHVBeyondPhase','cvUniqueAllBehaviorBeyondPhase', ...
    'cvDeltaR2AllBehaviorBeyondPhase','originalBetterPhaseOnly', ...
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

    if numel(unique(g)) == numel(P.morphOrder) && ...
            all(counts >= P.stats.minFishPerMorph)
        [pRaw(k), tbl] = kruskalwallis(x, g, 'off');
        H(k) = extractKWChiSquare(tbl);
        df(k) = numel(P.morphOrder) - 1;
    end
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
            pRaw(row) = ranksum(xa, xb);
        end
    end

    pAdjusted(rowsThisMetric) = pRaw(rowsThisMetric);
end

correctionFamily = repmat({''},nRows,1);
Pairwise = table(metric, metricLabel, analysisScope, groupA, groupB, ...
    nA, nB, medianA, medianB, medianDifference_BminusA, ciLow, ciHigh, ...
    cliffsDelta_BvsA, pRaw, pAdjusted, correctionFamily);
end

function [Omnibus,Pairwise] = applyConfiguredCorrections( ...
        Omnibus,Pairwise,Families,P)
% Each outcome contains two planned cavefish-versus-Surface contrasts.
% The default method is 'none', so they remain uncorrected. If 'bh' is
% requested as a sensitivity analysis, it is applied within one metric and
% analysis scope; biologically different outcomes are never pooled.
%
% Omnibus Kruskal-Wallis p-values are left raw by default because there is one
% three-group omnibus test per outcome. Cross-metric omnibus correction can be
% enabled explicitly with P.stats.correctOmnibusAcrossMetrics.
if ~isstruct(Families); error('Families must be a struct array.'); end
scopes = unique(Omnibus.analysisScope,'stable');
for s = 1:numel(scopes)
    for f = 1:numel(Families)
        idxO = strcmp(Omnibus.analysisScope,scopes{s}) & ...
            ismember(Omnibus.metric,Families(f).metrics);

        if P.stats.correctOmnibusAcrossMetrics
            Omnibus.pAdjusted(idxO) = adjustPValues( ...
                Omnibus.pRaw(idxO),P.stats.multipleComparisonMethod);
            Omnibus.correctionFamily(idxO) = ...
                repmat({[Families(f).name '_omnibus']},sum(idxO),1);
        else
            Omnibus.pAdjusted(idxO) = Omnibus.pRaw(idxO);
            Omnibus.correctionFamily(idxO) = ...
                repmat({'single_outcome_omnibus'},sum(idxO),1);
        end

        % Apply the configured method to the two planned comparisons for
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

function [FishROC,columns] = computeFishROCFromOOFTraces(All,P)
% Compute neuron-level ROC curves from truly held-out prediction traces,
% summarize neurons within fish, then retain fish as the inferential unit.
FishROC = table();
columns = struct('observed','','phaseOnly','','original','','behavior','','full','');
columns.observed = firstExistingVariable(All,P.roc.observedVariableCandidates);
columns.phaseOnly = firstExistingVariable( ...
    All,P.roc.phaseOnlyPredictionCandidates);
columns.original = firstExistingVariable(All,P.roc.originalPredictionCandidates);
columns.behavior = firstExistingVariable(All,P.roc.behaviorPredictionCandidates);
columns.full = firstExistingVariable(All,P.roc.fullPredictionCandidates);
if any(cellfun(@isempty,struct2cell(columns)))
    return;
end

fish = unique(All.fishKey,'stable');
morph = cell(numel(fish),1); session = cell(numel(fish),1);
fishKey = fish(:); nNeuronsROC = zeros(numel(fish),1);
aucPhaseOnly = nan(numel(fish),1); aucOriginal = nan(numel(fish),1);
aucBehavior = nan(numel(fish),1); aucFull = nan(numel(fish),1);
fprGrid = cell(numel(fish),1); tprPhaseOnly = cell(numel(fish),1);
tprOriginal = cell(numel(fish),1); tprBehavior = cell(numel(fish),1);
tprFull = cell(numel(fish),1);

for f = 1:numel(fish)
    rows = find(strcmp(All.fishKey,fish{f}) & All.baseQualityOK);
    first = find(strcmp(All.fishKey,fish{f}),1,'first');
    morph{f} = All.morph{first}; session{f} = All.session{first};
    curvesP = []; curvesO = []; curvesB = []; curvesF = [];
    aucP = []; aucO = []; aucB = []; aucF = [];
    for r = rows(:)'
        y = getTraceEntry(All,columns.observed,r);
        pP = getTraceEntry(All,columns.phaseOnly,r);
        pO = getTraceEntry(All,columns.original,r);
        pB = getTraceEntry(All,columns.behavior,r);
        pF = getTraceEntry(All,columns.full,r);
        lengths = [numel(y),numel(pP),numel(pO),numel(pB),numel(pF)];
        if any(lengths < 2) || numel(unique(lengths)) ~= 1; continue; end
        valid = isfinite(y) & isfinite(pP) & isfinite(pO) & ...
            isfinite(pB) & isfinite(pF);
        y = y(valid); pP = pP(valid); pO = pO(valid);
        pB = pB(valid); pF = pF(valid);
        label = y > P.roc.activityThresholdZ;
        if sum(label) < P.roc.minPositiveFrames || ...
                sum(~label) < P.roc.minNegativeFrames
            continue;
        end
        [curveP,aP] = binaryROCAtGrid(label,pP,P.roc.fprGrid);
        [curveO,aO] = binaryROCAtGrid(label,pO,P.roc.fprGrid);
        [curveB,aB] = binaryROCAtGrid(label,pB,P.roc.fprGrid);
        [curveF,aF] = binaryROCAtGrid(label,pF,P.roc.fprGrid);
        curvesP(end+1,:) = curveP; %#ok<AGROW>
        curvesO(end+1,:) = curveO; %#ok<AGROW>
        curvesB(end+1,:) = curveB; %#ok<AGROW>
        curvesF(end+1,:) = curveF; %#ok<AGROW>
        aucP(end+1,1) = aP; %#ok<AGROW>
        aucO(end+1,1) = aO; %#ok<AGROW>
        aucB(end+1,1) = aB; %#ok<AGROW>
        aucF(end+1,1) = aF; %#ok<AGROW>
    end
    nNeuronsROC(f) = size(curvesP,1);
    if nNeuronsROC(f) >= P.roc.minNeuronsPerFish
        fprGrid{f} = P.roc.fprGrid(:)';
        tprPhaseOnly{f} = median(curvesP,1,'omitnan');
        tprOriginal{f} = median(curvesO,1,'omitnan');
        tprBehavior{f} = median(curvesB,1,'omitnan');
        tprFull{f} = median(curvesF,1,'omitnan');
        aucPhaseOnly(f) = median(aucP,'omitnan');
        aucOriginal(f) = median(aucO,'omitnan');
        aucBehavior(f) = median(aucB,'omitnan');
        aucFull(f) = median(aucF,'omitnan');
    end
end

FishROC = table(morph,session,fishKey,nNeuronsROC,aucPhaseOnly,aucOriginal, ...
    aucBehavior,aucFull,fprGrid,tprPhaseOnly,tprOriginal,tprBehavior,tprFull);
keep = isfinite(FishROC.aucPhaseOnly) & isfinite(FishROC.aucOriginal) & ...
    isfinite(FishROC.aucBehavior) & isfinite(FishROC.aucFull);
FishROC = FishROC(keep,:);
end

function name = firstExistingVariable(T,candidates)
name = '';
variables = T.Properties.VariableNames;
for k = 1:numel(candidates)
    idx = find(strcmpi(variables,candidates{k}),1);
    if ~isempty(idx)
        name = variables{idx};
        return;
    end
end
end

function x = getTraceEntry(T,varName,row)
v = T.(varName);
x = [];
if iscell(v)
    x = v{row};
elseif isnumeric(v) || islogical(v)
    if ismatrix(v) && size(v,1) == height(T) && size(v,2) > 1
        x = v(row,:);
    end
end
if ~(isnumeric(x) || islogical(x)); x = []; return; end
x = double(x(:));
end

function [tprGrid,auc] = binaryROCAtGrid(label,score,fprGrid)
% ROC points respect tied prediction scores; AUC uses average ranks.
label = logical(label(:)); score = double(score(:));
[sortedScore,order] = sort(score,'descend');
sortedLabel = label(order); nPos = sum(label); nNeg = sum(~label);
lastInTie = [find(diff(sortedScore) ~= 0); numel(sortedScore)];
tp = cumsum(sortedLabel); fp = cumsum(~sortedLabel);
fpr = [0; fp(lastInTie)./nNeg; 1];
tpr = [0; tp(lastInTie)./nPos; 1];
[fprUnique,~,group] = unique(fpr);
tprUnique = accumarray(group,tpr,[],@max);
tprUnique = cummax(tprUnique);
tprGrid = interp1(fprUnique,tprUnique,fprGrid,'previous','extrap');
tprGrid = min(max(tprGrid,0),1);

% Mann-Whitney interpretation of AUC with half credit for ties.
[ascending,ordAsc] = sort(score,'ascend');
ranks = zeros(size(score)); i = 1;
while i <= numel(score)
    j = i;
    while j < numel(score) && ascending(j+1) == ascending(i); j = j+1; end
    ranks(ordAsc(i:j)) = mean(i:j);
    i = j+1;
end
auc = (sum(ranks(label)) - nPos*(nPos+1)/2) ./ (nPos*nNeg);
end

function [PairStats,FriedmanStats] = runROCModelStats(FishROC,P)
n = P.modelNames;
comparisons = { ...
    'aucOriginal','aucFull',n{1},n{2}; ...
    'aucOriginal','aucPhaseOnly',n{1},n{3}; ...
    'aucOriginal','aucBehavior',n{1},n{4}; ...
    'aucFull','aucPhaseOnly',n{2},n{3}; ...
    'aucFull','aucBehavior',n{2},n{4}; ...
    'aucPhaseOnly','aucBehavior',n{3},n{4}};
PairStats = runWithinMorphPairedStats(FishROC,comparisons,P);
scopes = [P.morphOrder,{'AllMorphs'}];
scope = scopes(:); nFish = zeros(numel(scopes),1);
chiSquare = nan(numel(scopes),1); pRaw = nan(numel(scopes),1);
kendallW = nan(numel(scopes),1);
for s = 1:numel(scopes)
    if strcmp(scopes{s},'AllMorphs')
        idx = true(height(FishROC),1);
    else
        idx = strcmp(FishROC.morph,scopes{s});
    end
    X = [FishROC.aucOriginal,FishROC.aucFull, ...
        FishROC.aucPhaseOnly,FishROC.aucBehavior];
    X = X(idx,:); X = X(all(isfinite(X),2),:); nFish(s) = size(X,1);
    if size(X,1) >= P.stats.minFishPerMorph
        [pRaw(s),tbl] = friedman(X,1,'off');
        chiSquare(s) = extractFriedmanChiSquare(tbl);
        if isfinite(chiSquare(s))
            kendallW(s) = chiSquare(s)/(size(X,1)*(size(X,2)-1));
        end
    end
end
pAdjusted = pRaw;
FriedmanStats = table(scope,nFish,chiSquare,pRaw,pAdjusted,kendallW);
end

function plotModelROC(FishROC,FriedmanStats,P,outputDir,baseName)
scopes = [{'AllMorphs'},P.morphOrder];
scopeLabels = [{'All fish'},P.morphOrder];
modelVars = {'tprOriginal','tprFull','tprPhaseOnly','tprBehavior'};
aucVars = {'aucOriginal','aucFull','aucPhaseOnly','aucBehavior'};
modelNames = P.modelNames;
colors = [0.42 0.42 0.42; 0.48 0.22 0.68; ...
    0.10 0.62 0.56; 0.20 0.52 0.82];
fig = figure('Visible',P.figure.visible,'Color','w', ...
    'Position',[50 50 920 760]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for s = 1:numel(scopes)
    ax = nexttile(tl); hold(ax,'on');
    if strcmp(scopes{s},'AllMorphs')
        idx = true(height(FishROC),1);
    else
        idx = strcmp(FishROC.morph,scopes{s});
    end
    plot(ax,[0 1],[0 1],':','Color',[0.58 0.58 0.58], ...
        'LineWidth',1.1,'HandleVisibility','off');
    legendText = cell(1,numel(modelVars));
    for j = 1:numel(modelVars)
        traceCells = FishROC.(modelVars{j})(idx);
        C = vertcat(traceCells{:});
        if isempty(C); continue; end
        curve = median(C,1,'omitnan');
        lo = prctile(C,25,1); hi = prctile(C,75,1);
        x = P.roc.fprGrid(:)';
        patch(ax,[x fliplr(x)],[lo fliplr(hi)],colors(j,:), ...
            'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
        plot(ax,x,curve,'Color',colors(j,:),'LineWidth',2.3, ...
            'DisplayName',modelNames{j});
        auc = FishROC.(aucVars{j})(idx);
        [med,ci] = bootstrapMedianCI(auc,P.stats.nBootstrap,P.stats.bootstrapCI);
        legendText{j} = sprintf('%s: AUC %.3f [%.3f, %.3f]', ...
            modelNames{j},med,ci(1),ci(2));
    end
    axis(ax,'square'); xlim(ax,[0 1]); ylim(ax,[0 1]); grid(ax,'on');
    set(ax,'TickDir','out','Box','off');
    xlabel(ax,'False-positive rate'); ylabel(ax,'True-positive rate');
    row = FriedmanStats(strcmp(FriedmanStats.scope,scopes{s}),:);
    if isempty(row); p = NaN; else; p = row.pAdjusted(1); end
    title(ax,sprintf('%s | paired model p=%s',scopeLabels{s},formatP(p)), ...
        'FontWeight','normal');
    legend(ax,legendText,'Location','southeast','FontSize',8);
end
title(tl,sprintf(['Held-out activity-state ROC | positive: observed z-score > %.2f' ...
    '\nCurves and AUC give equal weight to each fish'],P.roc.activityThresholdZ), ...
    'FontWeight','bold');
saveFigureBoth(fig,outputDir,baseName,P);
end

function plotPairedModelCV(Fish,PairStats,FriedmanStats,P,outputDir,baseName)
models = P.modelNames;
vars = P.modelMetricNames;
modelColors=[0.55 0.55 0.55; 0.45 0.20 0.70; ...
    0.10 0.62 0.56; 0.25 0.55 0.85];
fig=figure('Visible',P.figure.visible,'Color','w', ...
    'Position',[50 80 470*numel(P.morphOrder) 430]);
tl=tiledlayout(fig,1,numel(P.morphOrder),'TileSpacing','compact','Padding','compact');
for m=1:numel(P.morphOrder)
    ax=nexttile(tl); hold(ax,'on'); yline(ax,0,':','Color',[0.65 0.65 0.65]);
    idx=strcmp(Fish.morph,P.morphOrder{m});
    X=[Fish.(vars{1}),Fish.(vars{2}),Fish.(vars{3}),Fish.(vars{4})];
    X=X(idx,:); X=X(all(isfinite(X),2),:);
    for f=1:size(X,1)
        plot(ax,1:4,X(f,:),'-o','Color',mixWithWhite(P.morphColors(m,:),0.58), ...
            'MarkerFaceColor',mixWithWhite(P.morphColors(m,:),0.35), ...
            'MarkerEdgeColor','none','LineWidth',0.9,'HandleVisibility','off');
    end
    for j=1:4
        [med,ci]=bootstrapMedianCI(X(:,j),P.stats.nBootstrap,P.stats.bootstrapCI);
        plot(ax,[j j],ci,'k-','LineWidth',2.2,'HandleVisibility','off');
        scatter(ax,j,med,80,'d','MarkerFaceColor',modelColors(j,:), ...
            'MarkerEdgeColor','k','LineWidth',1.1,'HandleVisibility','off');
    end
    set(ax,'XTick',1:4,'XTickLabel',models,'XTickLabelRotation',18, ...
        'TickDir','out','Box','off'); xlim(ax,[0.65 4.35]); grid(ax,'on');
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
            'XTickLabel',{'phase + AHV','phase + AHV + forward + vigor'}, ...
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
    'Position', [60 60 440*nCols 360*nRows]);
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact','Padding','compact');

for k = 1:n
    ax = nexttile(tl); hold(ax,'on');
    plotOneFishMetric(ax, Fish, selected(k), Omnibus, Pairwise, P);
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
    bar(ax,m,med,0.68,'FaceColor',P.morphColors(m,:), ...
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

xlim(ax, [0.5 numel(P.morphOrder)+0.5]);
set(ax, 'XTick',1:numel(P.morphOrder), 'XTickLabel',P.morphOrder, ...
    'TickDir','out', 'Box','off');
ylabel(ax, M.label, 'Interpreter','tex');
grid(ax,'on');
if strcmp(M.summary,'fraction_positive') || startsWith(M.name,'fraction_')
    ylim(ax,[0 1]);
end

idxO = strcmp(Omnibus.metric, M.name);
if any(idxO)
    pKW = Omnibus.pAdjusted(find(idxO,1));
else
    pKW = NaN;
end
title(ax, sprintf('KW p=%s', formatP(pKW)), ...
    'FontWeight','normal');

idxP = strcmp(Pairwise.metric, M.name);
sub = Pairwise(idxP,:);
addMorphPairwiseBars(ax,sub,P.morphOrder);
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
    'Position',[30 30 380*nCols 305*nRows]);
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
            bar(ax,m,med,0.68,'FaceColor',P.morphColors(m,:), ...
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
        xlim(ax,[0.5 numel(P.morphOrder)+0.5]); grid(ax,'on');
        if c < nRows; set(ax,'XTickLabel',[]); end
        if k == 1; ylabel(ax,cls,'FontWeight','bold'); end
        idx = strcmp(Omnibus.analysisScope,cls) & strcmp(Omnibus.metric,M.name);
        p = NaN;
        if any(idx); p = Omnibus.pAdjusted(find(idx,1)); end
        title(ax,sprintf('%s | KW p=%s',M.label,formatP(p)), ...
            'FontWeight','normal','Interpreter','tex','FontSize',9);

        % Display both pre-specified cavefish-versus-Surface comparisons.
        % With the default settings these are the raw, uncorrected planned
        % p-values; pAdjusted is retained as the generic output column name.
        idxP = strcmp(Pairwise.analysisScope,cls) & ...
            strcmp(Pairwise.metric,M.name);
        sub = Pairwise(idxP,:);
        addMorphPairwiseBars(ax,sub,P.morphOrder);
    end
end
title(tl,['Class-stratified fish summaries (descriptive; class is model-derived)' newline ...
    sprintf('Minimum %d neurons per fish/class',P.classAnalysis.minNeuronsPerFishClass)], ...
    'Interpreter','none','FontWeight','bold');
saveFigureBoth(fig,outputDir,baseName,P);
end

function addMorphPairwiseBars(ax,sub,morphOrder)
% Add the two planned Surface-Molino and Surface-Pachon comparisons when
% their fish-level rank-sum p-values are available.
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
fprintf(fid,'Cross-morph matched four-model HD/AHV analysis\n\n');
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
fprintf(fid,['Model 2 - phase + AHV + forward + vigor: phase cosine + phase ' ...
    'sine + absolute AHV + signed AHV + forward bouts + vigor.\n']);
fprintf(fid,'Model 3 - phase: phase cosine + phase sine.\n');
fprintf(fid,['Model 4 - AHV + forward + vigor: absolute AHV + signed AHV + ' ...
    'forward bouts + vigor.\n']);
fprintf(fid,'Primary omnibus test: three-morph Kruskal-Wallis on fish summaries.\n');
fprintf(fid,'Planned contrasts: Molino-Surface and Pachon-Surface rank-sum tests.\n');
if strcmpi(P.stats.multipleComparisonMethod,'none')
    fprintf(fid,['Multiplicity control: none for the two explicitly planned ' ...
        'Surface-Molino and Surface-Pachon contrasts.\n']);
else
    fprintf(fid,['Multiplicity control: the two planned contrasts use %s ' ...
        'separately within each outcome and analysis scope.\n'], ...
        P.stats.multipleComparisonMethod);
end
if P.stats.correctOmnibusAcrossMetrics
    fprintf(fid,['Omnibus p-values are additionally corrected across metrics ' ...
        'within each prespecified figure family.\n']);
else
    fprintf(fid,['Omnibus p-values are unadjusted: each is the single three-morph ' ...
        'test for its biological outcome.\n']);
end
fprintf(fid,'Effect: group B minus group A difference in fish-level medians.\n');
fprintf(fid,'Effect CI: percentile bootstrap resampling fish within morph (%d draws).\n', ...
    P.stats.nBootstrap);
fprintf(fid,'Within-morph model tests are paired sign-rank tests on fish summaries.\n');
fprintf(fid,['Nested phase-amplitude and phase + AHV versus phase + AHV + ' ...
    'forward + vigor coefficient tests are paired at the fish level.\n']);
fprintf(fid,['All upstream-selected phase-tuned neurons passing numerical fit QC ' ...
    'were included; no downstream coefficient-significance filter was applied.\n']);
fprintf(fid,'Raw b1/b2 units: activity SD per rad/s when the upstream fitting script was unchanged.\n');
fprintf(fid,'b1 signed = net speed effect; abs(b1) = complementary speed encoding strength.\n');
fprintf(fid,'Standardized bF/bV are primary across morphs because raw F/V scaling may differ.\n');
fprintf(fid,'Negative cross-validated R^2 and unique contributions were retained.\n');
fprintf(fid,['cvUniquePhase = 1 - SSE_(phase+AHV+forward+vigor) / ' ...
    'SSE_(AHV+forward+vigor).\n']);
fprintf(fid,['cvUniqueAddedBehavior = 1 - SSE_(phase+AHV+forward+vigor) / ' ...
    'SSE_(phase+AHV).\n']);
fprintf(fid,['cvUniqueAHVBeyondPhase = 1 - SSE_(phase+AHV) / ' ...
    'SSE_phase.\n']);
fprintf(fid,['cvUniqueAllBehaviorBeyondPhase = 1 - ' ...
    'SSE_(phase+AHV+forward+vigor) / SSE_phase.\n']);
fprintf(fid,['All four CV scores and all unique-contribution summaries use only ' ...
    'neurons with finite matched outputs for every model.\n']);
fprintf(fid,['Matched nested contrasts add AHV to phase, add all behavior ' ...
    'to phase, add forward + vigor to phase + AHV, and add phase to ' ...
    'AHV + forward + vigor.\n']);
fprintf(fid,['Absolute phase-model performance is conditional on selecting ' ...
    'phase-tuned neurons upstream; nested comparisons use the same selected set.\n']);
fprintf(fid,'VIF and condition number are diagnostics unless exclusion is explicitly enabled.\n');
fprintf(fid,'Class-stratified results are secondary/descriptive because class was model-derived.\n');
fprintf(fid,['ROC is computed only from out-of-fold observed and predicted traces. ' ...
    'Positive activity state: observed z-score > %.3f.\n'],P.roc.activityThresholdZ);
fprintf(fid,['ROC/AUC is a secondary thresholded activity-state analysis and does ' ...
    'not replace the continuous blocked-CV R^2 analysis.\n']);
end

function x = wrapToPiLocal(x)
x = mod(x+pi,2*pi)-pi;
end
