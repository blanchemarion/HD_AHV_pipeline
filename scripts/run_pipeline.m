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
%   1) fit three models and classify phase-tuned neurons;
%   2) compare fitted parameters across morphs;
%   3) compare core behavior across morphs;
%   4) reproduce poster candidate-cell figures;
%   5) compare class abundance, anatomy, AHV response, speed, behavior and
%      preferred-direction features across morphs.

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
        fprintf('[1/5] ForceModelRefit=true; fitting three matched models for all morphs...');
    else
        fprintf('[1/5] Processed model inputs are incomplete; fitting three matched models for all morphs...');
        fprintf('%s', newline);
        for iMissing = 1:numel(cachedModelFits.missing)
            fprintf('  %s', cachedModelFits.missing{iMissing});
            fprintf('%s', newline);
        end
    end
    fprintf('%s', newline);
    results.modelFits = ...
        fit_neuron_HD_AHV_three_models_behavior_tuned_cells_only_ROC(cfg);
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


fprintf('\n[3/5] Comparing core behavior across morphs...\n');
results.coreBehaviorComparison = compare_core_behavior_across_morphs(cfg);

fprintf('\n[4/5] Reproducing poster candidate-cell figures...\n');
results.posterFigures = reproduce_poster_fig3A_fig4E_H_candidate_cells_v2( ...
    'DataRoot', cfg.DataRoot, ...
    'OutputDir', cfg.PosterFigureDir, ...
    'Visible', cfg.FigureVisible);


fprintf('\n[5/5] Comparing class features across morphs...\n');
classifiedFiles = dir(fullfile(cfg.ClassifiedCandidateDir, '**', ...
    '*_candidate_neurons_clean_classified.mat'));
classifiedFiles = arrayfun(@(x) fullfile(x.folder, x.name), ...
    classifiedFiles, 'UniformOutput', false);
if isempty(classifiedFiles)
    error(['No centralized classified candidate files were produced. ' ...
        'Inspect the model-fitting output before running the class comparison.']);
end

results.classFeatureComparison = compare_class_features_across_morphs( ...
    'Files', classifiedFiles, ...
    'DataRoot', cfg.DataRoot, ...
    'OutputDir', fullfile(cfg.ClassComparisonDir,'specified_features'), ...
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
    modelFile = fullfile(cfg.ModelDataDir, sprintf( ...
        'HD_AHV_three_models_phase_tuned_only_results_%s.mat', ...
        morph.name));
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
            elseif ~all(ismember({'P','sessions'}, names))
                cache.missing{end+1,1} = sprintf( ...
                    '%s cached results lack cohort/QC metadata: %s', ...
                    morph.display, modelFile);
            else
                cached = load(modelFile, 'P', 'sessions', 'AllNeurons', 'Results');
                sessionFailure = ~isfield(cached,'Results') || ~iscell(cached.Results) || ...
                    any(cellfun(@(r) isstruct(r) && isfield(r,'error'), cached.Results));
                if sessionFailure
                    cache.missing{end+1,1} = sprintf( ...
                        '%s cached results contain a failed or missing session: %s', ...
                        morph.display, modelFile);
                end
                requiredFitColumns = {'pA','pB1','pB2','nFitSamples', ...
                    'failureReason','pB2Adjusted','class', ...
                    'cvR2PhaseOnly','cvR2Behavior','cvR2Full'};
                validFitTable = isfield(cached,'AllNeurons') && ...
                    istable(cached.AllNeurons) && height(cached.AllNeurons)>0 && ...
                    all(ismember(requiredFitColumns, ...
                    cached.AllNeurons.Properties.VariableNames));
                if ~validFitTable
                    cache.missing{end+1,1} = sprintf( ...
                        '%s cached results lack completed pooled classification: %s', ...
                        morph.display, modelFile);
                elseif any(strcmp(string(cached.AllNeurons.class),'PendingPooledBonferroni'))
                    cache.missing{end+1,1} = sprintf( ...
                        '%s cached results still have pending pooled labels: %s', ...
                        morph.display, modelFile);
                end
                threeModelCache = isfield(cached.P, 'model') && ...
                    isfield(cached.P.model, 'version') && ...
                    strcmp(cached.P.model.version, ...
                    'three_models_fixed_source_pref_turn_bias_v2');
                responseScaleMatches = isfield(cached.P.model, 'zscoreActivityWithinFitWindow') && ...
                    logical(cached.P.model.zscoreActivityWithinFitWindow);
                if ~responseScaleMatches
                    cache.missing{end+1,1} = sprintf( ...
                        '%s cached response normalization does not match the current fitter: %s', ...
                        morph.display, modelFile);
                end
                if ~threeModelCache
                    cache.missing{end+1,1} = sprintf( ...
                        '%s cached results do not match the current three-model design: %s', ...
                        morph.display, modelFile);
                end
                rotationGateMatches = isfield(cached.P, 'lowTurnQC') && ...
                    isfield(cached.P.lowTurnQC, 'enabled') && ...
                    ~logical(cached.P.lowTurnQC.enabled) && ...
                    isfield(cached.P.lowTurnQC, 'minRotationsEachDirection') && ...
                    cached.P.lowTurnQC.minRotationsEachDirection==0;
                if ~rotationGateMatches
                    cache.missing{end+1,1} = sprintf( ...
                        '%s cached results do not match the current rotation-QC settings: %s', ...
                        morph.display, modelFile);
                end
                expectedNames = recordingFolderNames(morph.dataDir);
                cachedNames = string({cached.sessions.name});
                if ~isequal(sort(lower(cachedNames(:))), ...
                        sort(lower(expectedNames(:))))
                    cache.missing{end+1,1} = sprintf( ...
                        ['%s cached cohort does not cover every rec folder ' ...
                         '(cached %d, current %d): %s'], morph.display, ...
                        numel(cachedNames), numel(expectedNames), modelFile);
                end
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

summaryFile=fullfile(cfg.ModelDataDir,'three_model_fit_pipeline_summary.mat');
try
    summary=load(summaryFile,'pipelineResults');
    validSummary=isfield(summary,'pipelineResults') && ...
        isfield(summary.pipelineResults,'pooledBonferroni') && ...
        isfield(summary.pipelineResults.pooledBonferroni,'nTestsB2') && ...
        summary.pipelineResults.pooledBonferroni.nTestsB2>0;
catch
    validSummary=false;
end
if ~validSummary
    cache.missing{end+1,1} = sprintf('Missing or incomplete pooled Bonferroni summary: %s',summaryFile);
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
    recNames = recordingFolderNames(morph.dataDir);
    if isempty(recNames)
        error('%s has no rec* folders under: %s', morph.display, morph.dataDir);
    end
    foldersToValidate = fullfile(string(morph.dataDir), recNames);
    for f = 1:numel(foldersToValidate)
        candidate = dir(fullfile(char(foldersToValidate(f)), ...
            '*_candidate_neurons_clean.mat'));
        candidateNames = lower(string({candidate.name}));
        candidate = candidate(~contains(candidateNames, ...
            ["_classified","_direct_model"]));
        if isempty(candidate)
            error('%s is missing a clean candidate file.', ...
                char(foldersToValidate(f)));
        end
        if requireModelInputs
            behavior = dir(fullfile(char(foldersToValidate(f)), ...
                '*swimResults_pass2*.mat'));
            if isempty(behavior)
                error('%s is missing DLC pass-2 behavior.', ...
                    char(foldersToValidate(f)));
            end
        end
        variables = whos('-file', fullfile(candidate(1).folder, candidate(1).name));
        missing = setdiff(requiredCandidateVariables, {variables.name});
        if ~isempty(missing)
            error('%s is missing variables: %s', ...
                fullfile(candidate(1).folder,candidate(1).name), ...
                strjoin(missing, ', '));
        end
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

function names = recordingFolderNames(morphDir)
d = dir(fullfile(morphDir, 'rec*'));
d = d([d.isdir]);
names = string({d.name});
names = names(:);
end
