# HD/AHV analysis pipeline

MATLAB pipeline for fitting neural activity encoding models in Surface, Molino, and Pachon fish using network phase and angular head velocity (AHV). It selects phase-tuned candidates, fits three matched models, classifies CW/CCW/Symmetric neurons, and produces fish-level model, behavior, and functional/anatomical comparisons.

Raw data are read from an external tree /home/blanche/data/{surface,molino,pachon}/rec*. Generated data and tables go to `data_processed/`; figures and reports go to `outputs/`.

## Workflow

`scripts/run_pipeline.m` runs five stages:

1. `fit_neuron_HD_AHV_three_models_behavior_tuned_cells_only_ROC` selects phase-tuned candidates, fits models, applies pooled classification, and exports classified copies.
2. `compare_HD_AHV_model_with_phase_only_behavior_params_across_morphs_v3` compares coefficients and held-out performance.
3. `compare_core_behavior_across_morphs` compares original pass-2 turning behavior.
4. `reproduce_poster_fig3A_fig4E_H_candidate_cells_v2` reproduces candidate-cell poster analyses.
5. `compare_class_features_across_morphs` compares CW, CCW, and Symmetric abundance, anatomy, AHV responses, speed relationships, and preferred direction.

Configuration is in `scripts/pipeline_config.m`; shared fitter helpers are in `scripts/models/+hdahv/`.

## Requirements

- MATLAB R2020b or newer
- Statistics and Machine Learning Toolbox
- A data root containing `surface/`, `molino/`, and `pachon/`

## Inputs

The default data root is `/home/blanche/data`:

```text
<DataRoot>/
  surface/rec*/
  molino/rec*/
  pachon/rec*/
```

Every `rec*` directory is a recording/session and must contain one source `*_candidate_neurons_clean.mat` and one `*_swimResults_pass2*.mat`. Files containing `_classified` or `_direct_model` are not source candidates. Candidate MAT files must include `calcium_traces`, `time_s`, `network_phase_rad`, and `candidate_cell_ids`.

The class-feature stage also locates recording-level raster, all-cell ROI, and `fov_correction.mat` inputs to map candidates to original ROIs and aligned anatomy. Mapping failures are recorded in its audit.

## Run

From the repository root:

```matlab
addpath('scripts')
results = run_pipeline;
```

Common options:

```matlab
results = run_pipeline('DataRoot', '/path/to/data');
results = run_pipeline('FigureVisible', 'on');
results = run_pipeline('ForceModelRefit', true);
results = run_pipeline( ...
    'DataRoot', '/path/to/data', ...
    'PipelineRoot', '/path/to/output/project');
```

Accepted options are `DataRoot`, `PipelineRoot`, `FigureVisible`, `OverwriteOutputs`, and `ForceModelRefit`. `pipeline_config` creates required output directories.

## Model and classification

Phase-tuned candidates are selected separately within each fish using `ORI_V15.m` procedure: timewise z-scoring, framewise 2nd–98th percentile clipping and population centering, positive rectification, 36-bin circularly smoothed Skaggs information, 500 circular shifts over 10–90% of the valid sequence, and within-fish Benjamini–Hochberg FDR (`q < 0.05`). A candidate can contribute to the population phase used in its own test, so selection is partly circular unless upstream phase was constructed leave-one-out.

For each selected neuron, the combined fixed-preference model is:

```text
y(t) = betaTheta * cos(sourcePrefRad - phase(t))
     + b1 * abs(AHV(t)) + b2 * AHV(t) + intercept + error(t)
```

`sourcePrefRad` comes from the stored source tuning curve before behavior is loaded and is not re-estimated. Three models use identical valid samples and the same five contiguous blocked cross-validation folds:

- phase + AHV: fixed phase, `abs(AHV)`, and signed `AHV`
- phase only: fixed phase
- AHV only: `abs(AHV)` and signed `AHV`

AHV events are convolved with a causal calcium kernel (`tau = 5 s`, duration `20 s`). Turns at or below `0.239 rad` absolute angle are excluded. Positive AHV denotes CCW. Directional rotation counts are descriptive QC only and do not exclude recordings.

After all morphs finish, coefficient p-values are pooled across the cohort. Bonferroni adjustment is applied separately to the `betaTheta`, `b1`, and `b2` families. CW/CCW/Symmetric labels derive from signed `b2` and are saved in morph results and non-destructive classified copies.

## Cache behavior

Models are reused only when the current output contract is complete. `run_pipeline` checks model version `three_models_fixed_source_pref_turn_bias_v2`, fit/CV/classification fields, absence of failed sessions, current normalization and rotation-QC settings, exact coverage of current `rec*` folders, classified exports for all morphs, and a valid pooled-Bonferroni summary.

If any check fails, all morphs are refitted. `ForceModelRefit=true` forces this explicitly.

## Downstream analyses

### Model comparison

Fish/animal—not neuron or CV fold—is the inferential unit. The analysis summarizes coefficients and held-out R² within fish, compares paired models, evaluates phase and AHV contributions, and reports planned nonparametric morph contrasts with effect sizes and bootstrap intervals. Optional recording exclusions are disabled by default and retained only for sensitivity analysis.

### Core behavior

This stage reads the complete original pass-2 file for each successfully fitted session. It reports eligible-turn frequency, median absolute turn angle, median quiet interval, and 90th-percentile quiet interval. Quiet intervals lie between consecutive eligible turns. Plots show morph medians, fish-bootstrap 95% intervals, and individual fish.

### Functional classes and anatomy

The final stage compares CW, CCW, and Symmetric populations. Candidate traces are globally matched to original raster ROIs; anatomy comes from the corresponding all-cell file and is aligned with the FOV correction. The default coordinate system is normalized FOV space. Anatomical density is normalized within fish so cell-rich fish do not carry greater weight.

## Run individual stages

Add `scripts`, `scripts/models`, and `scripts/visualize` to the MATLAB path, then create `cfg = pipeline_config(...)`. The five stage functions listed under Workflow can then be called directly. Model, behavior, and class comparisons require fitter outputs. Poster reproduction reads source candidates directly and can run independently. Function headers document stage-specific options such as explicit file lists, phase alignment, coordinate mode, and output directories.

Example:

```matlab
addpath('scripts', 'scripts/models', 'scripts/visualize')
cfg = pipeline_config('DataRoot', '/path/to/data');
fitResults = fit_neuron_HD_AHV_three_models_behavior_tuned_cells_only_ROC(cfg);
modelResults = compare_HD_AHV_model_with_phase_only_behavior_params_across_morphs_v3(cfg);
behaviorResults = compare_core_behavior_across_morphs(cfg);
posterResults = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
    'DataRoot', cfg.DataRoot, 'OutputDir', cfg.PosterFigureDir);
```

## Outputs

```text
data_processed/
  model_fits/
    HD_AHV_three_models_phase_tuned_only_results_<morph>.mat
    HD_AHV_three_models_phase_tuned_only_neurons_<morph>.csv
    HD_AHV_three_models_phase_tuned_only_fish_summary_<morph>.csv
    HD_AHV_phase_tuning_all_candidates_<morph>.csv
    HD_AHV_phase_tuning_fish_summary_<morph>.csv
    three_model_fit_pipeline_summary.mat
    sessions/<morph>/<session>_three_models_phase_tuned_only.mat
  classified_candidates/<morph>/.../*_candidate_neurons_clean_classified.mat
  candidate_raster_mappings/
  pipeline_run_summary.mat

outputs/
  model_fits/<morph>/
  model_comparison/
  core_behavior_comparison/
  poster_figures/
  class_comparison_functional/specified_features/
```

`model_comparison/` contains result MAT data, analysis notes, fish-level tables, and PNG/SVG figures. `core_behavior_comparison/` contains fish and bout audits, estimates, tests, `core_behavior_results.mat`, and figures. `poster_figures/` contains metric tables, alignment audits, and figures. The class directory contains cell/fish summaries, ROI/anatomy audits, statistical tables, normalized density data, `class_feature_comparison.mat`, and figures.

`data_processed/pipeline_run_summary.mat` stores the returned results structure and effective configuration.

## Reproducibility note

Most scientific parameters are defined near the top of their analysis function, not in `pipeline_config.m`. Retain saved configurations and audits with any alternative analysis.
