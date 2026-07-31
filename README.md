# HD/AHV analysis pipeline

This folder is the cleaned, self-contained entry point for the five selected
analyses. Raw data remain read-only in the external data tree; generated
tables and MAT files go to `data_processed/`, and figures/reports go to
`outputs/`.

## Requirements

- MATLAB R2020b or newer.
- Statistics and Machine Learning Toolbox (used by the model and statistical
  analyses).
- A data root containing `surface/`, `molino/`, and `pachon/`.
- Each recording folder must contain one
  `*_candidate_neurons_clean.mat` file and one DLC pass-2
  `*_swimResults_pass2*.mat` file.
- The anatomy step additionally expects
  `<morph>/HD_AHV_neuron_model_outputs_candidate_clean/HD_AHV_model_results.mat`.

## Run

From MATLAB:

```matlab
cd('/path/to/HD_AHV_pipeline')
addpath(genpath('scripts'))
results = run_pipeline;
```

The default raw-data root is `/home/blanche/data`. To use another location:

```matlab
results = run_pipeline('DataRoot', '/path/to/data');
```

For a visible interactive run:

```matlab
results = run_pipeline('FigureVisible', 'on');
```

All paths are centralized in `scripts/pipeline_config.m`; no analysis script
needs machine-specific path edits.

## Run individual stages

```matlab
cfg = pipeline_config;
fitResults = fit_neuron_HD_AHV_four_models_behavior_tuned_cells_only_ROC(cfg);
modelStats = compare_HD_AHV_model_with_phase_only_behavior_params_across_morphs_v3(cfg);
anatomyStats = compare_anatomical_populations_across_morphs_adapted_stats(cfg);
poster = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
    'DataRoot', cfg.DataRoot, 'OutputDir', cfg.PosterFigureDir);
```

The class comparison consumes the centralized classified copies produced by
the fit stage. `run_pipeline` supplies those files automatically.

## Output layout

```text
data_processed/
  model_fits/
    HD_AHV_behavior_augmented_phase_tuned_only_with_phase_only_results_surface.mat
    HD_AHV_behavior_augmented_phase_tuned_only_with_phase_only_results_molino.mat
    HD_AHV_behavior_augmented_phase_tuned_only_with_phase_only_results_pachon.mat
    sessions/<morph>/
  classified_candidates/<morph>/
  pipeline_run_summary.mat

outputs/
  model_fits/<morph>/
  model_comparison/
  anatomy_comparison/
  poster_figures/
  class_comparison_functional/
```

Model figures contain the morph in both the filename and figure title.
Cross-morph figures live in one stable output directory per analysis.

## Reproducibility notes

- The fit uses identical valid samples and contiguous blocked-CV folds for
  the four models.
- Phase tuning is selected independently within each fish with the existing
  ORI_V15 shuffle test and within-fish BH-FDR.
- CW/CCW/Symmetric labels are produced by the augmented model and stored in
  the morph result files and centralized classified candidate copies.
- Cross-morph inference uses fish/session as the biological replicate.
- The previous hard-coded exclusions (`rec44`, `rec70`, `rec52`, `rec80`,
  `rec81`) are not applied by default because they were based on an outcome.
  They remain available in the comparison script as an explicitly labelled
  sensitivity analysis.
