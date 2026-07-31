function results = run_pipeline(varargin)
%RUN_PIPELINE Run the complete HD/AHV analysis in dependency order.
%
%   results = run_pipeline
%   results = run_pipeline('DataRoot', '/path/to/data')
%   results = run_pipeline('ForceModelRefit', true)
%
% Complete model and classification outputs already present under
% data_processed/ are reused by default. Use ForceModelRefit to rebuild them.
%
% Steps:
%   1) fit four models and classify phase-tuned neurons;
%   2) compare fitted parameters across morphs;
%   3) compare anatomical populations;
%   4) reproduce poster candidate-cell figures;
%   5) compare CW, CCW and Symmetric cells.

scriptsDir = fileparts(mfilename('fullpath'));
addpath(scriptsDir);
addpath(fullfile(scriptsDir, 'models'));
addpath(fullfile(scriptsDir, 'visualize'));

cfg = pipeline_config(varargin{:});
[reuseModels, cachedModelFits] = findReusableModelOutputs(cfg);
fitModels = cfg.ForceModelRefit || ~reuseModels;
validatePipelineInputs(cfg, fitModels);
results = struct();

if fitModels
    fprintf('%s', newline);
    if cfg.ForceModelRefit
        fprintf('[1/5] ForceModelRefit=true; fitting four matched models for all morphs...');
    else
        fprintf('[1/5] Processed model inputs are incomplete; fitting four matched models for all morphs...');
        fprintf('%s', newline);
        for iMissing = 1:numel(cachedModelFits.missing)
            fprintf('  %s', cachedModelFits.missing{iMissing});
            fprintf('%s', newline);
        end
    end
    fprintf('%s', newline);
    results.modelFits = ...
        fit_neuron_HD_AHV_four_models_behavior_tuned_cells_only_ROC(cfg);
else
    fprintf('%s', newline);
    fprintf('[1/5] Reusing complete model outputs in data_processed; model fitting skipped.');
    fprintf('%s', newline);
    results.modelFits = cachedModelFits;
    results.modelFits.skipped = true;
end

fprintf('\n[2/5] Comparing model parameters across morphs...\n');
results.modelComparison = ...
    compare_HD_AHV_model_with_phase_only_behavior_params_across_morphs_v3(cfg);

fprintf('\n[3/5] Comparing anatomical populations across morphs...\n');
results.anatomyComparison = ...
    compare_anatomical_populations_across_morphs_adapted_stats(cfg);

fprintf('\n[4/5] Reproducing poster candidate-cell figures...\n');
results.posterFigures = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
    'DataRoot', cfg.DataRoot, ...
    'OutputDir', cfg.PosterFigureDir, ...
    'Visible', cfg.FigureVisible);

fprintf('\n[5/5] Comparing CW, CCW and Symmetric classes...\n');
classifiedFiles = dir(fullfile(cfg.ClassifiedCandidateDir, '**', ...
    '*_candidate_neurons_clean_classified.mat'));
classifiedFiles = arrayfun(@(x) fullfile(x.folder, x.name), ...
    classifiedFiles, 'UniformOutput', false);
if isempty(classifiedFiles)
    error(['No centralized classified candidate files were produced. ' ...
        'Inspect the model-fitting output before running the class comparison.']);
end
results.classComparison = compare_CW_CCW_Symmetric_across_morphs_v4( ...
    'Files', classifiedFiles, ...
    'OutputDir', cfg.ClassComparisonDir, ...
    'Visible', cfg.FigureVisible);

save(fullfile(cfg.DataProcessedDir, 'pipeline_run_summary.mat'), ...
    'results', 'cfg', '-v7.3');
fprintf('\nPipeline complete. Outputs:\n  %s\n', cfg.OutputsDir);
end


function [isComplete, cache] = findReusableModelOutputs(cfg)
cache = struct();
cache.skipped = false;
cache.reason = '';
cache.modelFiles = cell(numel(cfg.Morphs),1);
cache.classifiedFiles = {};
cache.missing = {};

for m = 1:numel(cfg.Morphs)
    morph = cfg.Morphs(m);
    modelFile = fullfile(cfg.ModelDataDir, sprintf([ ...
        'HD_AHV_behavior_augmented_phase_tuned_only_' ...
        'with_phase_only_results_%s.mat'], morph.name));
    cache.modelFiles{m} = modelFile;
    if ~isfile(modelFile)
        cache.missing{end+1,1} = sprintf( ...
            'Missing %s model results: %s', morph.display, modelFile);
    else
        try
            variables = whos('-file', modelFile);
            names = {variables.name};
            if ~any(ismember({'AllNeurons','Results'}, names))
                cache.missing{end+1,1} = sprintf( ...
                    '%s has neither AllNeurons nor Results: %s', ...
                    morph.display, modelFile);
            end
        catch ME
            cache.missing{end+1,1} = sprintf( ...
                'Cannot inspect %s model results (%s): %s', ...
                morph.display, ME.message, modelFile);
        end
    end

    classified = dir(fullfile(cfg.ClassifiedCandidateDir, morph.name, ...
        '**', '*_candidate_neurons_clean_classified.mat'));
    if isempty(classified)
        cache.missing{end+1,1} = sprintf( ...
            'No %s classified candidate files under: %s', ...
            morph.display, fullfile(cfg.ClassifiedCandidateDir, morph.name));
        continue;
    end
    for f = 1:numel(classified)
        classifiedFile = fullfile(classified(f).folder, classified(f).name);
        cache.classifiedFiles{end+1,1} = classifiedFile;
        try
            variables = whos('-file', classifiedFile);
            names = {variables.name};
            if ~all(ismember({'candidate_cell_ids','neuron_class_label'}, names))
                cache.missing{end+1,1} = sprintf( ...
                    'Incomplete classified candidate file: %s', ...
                    classifiedFile);
            end
        catch ME
            cache.missing{end+1,1} = sprintf( ...
                'Cannot inspect classified candidate file (%s): %s', ...
                ME.message, classifiedFile);
        end
    end
end

isComplete = isempty(cache.missing);
if isComplete
    cache.reason = ['Complete morph-level model results and classified ' ...
        'candidate files already exist under data_processed.'];
end
end


function validatePipelineInputs(cfg, requireModelInputs)
requiredCandidateVariables = { ...
    'calcium_traces','time_s','network_phase_rad','candidate_cell_ids'};
for m = 1:numel(cfg.Morphs)
    morph = cfg.Morphs(m);
    if exist(morph.dataDir, 'dir') ~= 7
        error('Missing morph data directory: %s', morph.dataDir);
    end
    candidates = dir(fullfile(morph.dataDir, '**', ...
        '*_candidate_neurons_clean.mat'));
    candidateFolders = unique(string({candidates.folder}));
    if isempty(candidateFolders)
        error('%s has no clean candidate files under: %s', ...
            morph.display, morph.dataDir);
    end
    if requireModelInputs
        behaviors = dir(fullfile(morph.dataDir, '**', ...
            '*swimResults_pass2*.mat'));
        behaviorFolders = unique(string({behaviors.folder}));
        foldersToValidate = intersect(candidateFolders, behaviorFolders);
        if isempty(foldersToValidate)
            error(['%s has no recording folder containing both a clean ' ...
                'candidate file and DLC pass-2 behavior.'], morph.display);
        end
    else
        foldersToValidate = candidateFolders;
    end
    for f = 1:numel(foldersToValidate)
        candidate = dir(fullfile(char(foldersToValidate(f)), ...
            '*_candidate_neurons_clean.mat'));
        variables = whos('-file', fullfile(candidate(1).folder, candidate(1).name));
        missing = setdiff(requiredCandidateVariables, {variables.name});
        if ~isempty(missing)
            error('%s is missing variables: %s', ...
                fullfile(candidate(1).folder,candidate(1).name), ...
                strjoin(missing, ', '));
        end
    end
    if ~isfile(morph.anatomyModelFile)
        error('Missing anatomy model for %s: %s', ...
            morph.display, morph.anatomyModelFile);
    end
    if requireModelInputs
        fprintf('Preflight %s: %d model-ready paired recording folders.', ...
            morph.display, numel(foldersToValidate));
    else
        fprintf(['Preflight %s: %d candidate recording folders; ' ...
            'DLC/model-input checks skipped because cached models are reusable.'], ...
            morph.display, numel(foldersToValidate));
    end
    fprintf('%s', newline);
end
end
