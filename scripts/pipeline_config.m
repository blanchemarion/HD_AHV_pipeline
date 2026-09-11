function cfg = pipeline_config(varargin)
%PIPELINE_CONFIG Central paths and run settings for the HD/AHV pipeline.
%
% Edit DataRoot below, or override it without editing this file:
%   cfg = pipeline_config('DataRoot', '/path/to/data');
%
% DataRoot must contain the morph folders surface/, molino/ and pachon/.

pipelineRoot = fileparts(fileparts(mfilename('fullpath')));

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'DataRoot', '/home/blanche/data', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'PipelineRoot', pipelineRoot, ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'FigureVisible', 'off', ...
    @(x) any(strcmpi(string(x), ["on","off"])));
addParameter(p, 'OverwriteOutputs', true, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'ForceModelRefit', false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});

cfg = p.Results;
cfg.DataRoot = char(cfg.DataRoot);
cfg.PipelineRoot = char(cfg.PipelineRoot);
cfg.FigureVisible = char(cfg.FigureVisible);
cfg.OverwriteOutputs = logical(cfg.OverwriteOutputs);
cfg.ForceModelRefit = logical(cfg.ForceModelRefit);

cfg.Morphs = struct( ...
    'name', {'surface','molino','pachon'}, ...
    'display', {'Surface','Molino','Pachon'});
for i = 1:numel(cfg.Morphs)
    cfg.Morphs(i).dataDir = fullfile(cfg.DataRoot, cfg.Morphs(i).name);
end

cfg.DataProcessedDir = fullfile(cfg.PipelineRoot, 'data_processed');
cfg.OutputsDir = fullfile(cfg.PipelineRoot, 'outputs');

cfg.ModelDataDir = fullfile(cfg.DataProcessedDir, 'model_fits');
cfg.SessionModelDataDir = fullfile(cfg.ModelDataDir, 'sessions');
cfg.CandidateRasterMappingDir = fullfile(cfg.DataProcessedDir, 'candidate_raster_mappings');
cfg.ClassifiedCandidateDir = ...
    fullfile(cfg.DataProcessedDir, 'classified_candidates');

cfg.ModelFigureDir = fullfile(cfg.OutputsDir, 'model_fits');
cfg.ModelComparisonDir = fullfile(cfg.OutputsDir, 'model_comparison');
cfg.PosterFigureDir = fullfile(cfg.OutputsDir, 'poster_figures');
cfg.CoreBehaviorComparisonDir = fullfile(cfg.OutputsDir, ...
    'core_behavior_comparison');
cfg.ClassComparisonDir = fullfile(cfg.OutputsDir, ...
    'class_comparison_functional');

requiredDirs = {cfg.DataProcessedDir, cfg.OutputsDir, cfg.ModelDataDir, ...
    cfg.SessionModelDataDir, cfg.CandidateRasterMappingDir, cfg.ClassifiedCandidateDir, ...
    cfg.ModelFigureDir, cfg.ModelComparisonDir, ...
    cfg.PosterFigureDir, cfg.CoreBehaviorComparisonDir, ...
    cfg.ClassComparisonDir};
for i = 1:numel(requiredDirs)
    if exist(requiredDirs{i}, 'dir') ~= 7
        mkdir(requiredDirs{i});
    end
end
end
