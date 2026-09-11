function results = compare_class_features_across_morphs(varargin)
%COMPARE_CLASS_FEATURES_ACROSS_MORPHS Fish-level CW/CCW/Symmetric analysis.
%
% Anatomical coordinates are read from the recording-specific ALL_CELLS
% file and aligned with that recording's fov_correction.mat. The expected
% directory convention is:
%
%   /home/blanche/data/surface/rec53
%   /home/blanche/data/molino/rec25
%   /home/blanche/data/pachon/rec137
%
% The saved FOV rotation is applied in raw image coordinates after
% subtracting the image centre. The displayed axes are positive rightward
% (mediolateral, ML) and positive rostrally (rostrocaudal, RC).
%
% RASTER and ALL_CELLS are treated as sharing the original ROI index. Every
% classified candidate trace is searched against every original RASTER
% column. A classified candidate ID is called an original ROI ID only when
% these global best matches reproduce the candidate IDs. The globally
% matched original ROI IDs are then used to retrieve anatomy from
% ALL_CELLS.
%
% Anatomical inference is performed on per-fish summaries. Fish-normalized
% density maps give every fish equal total weight, irrespective of how many
% neurons were classified in that fish.

baseCfg = pipeline_config();

p = inputParser;
addParameter(p,'Files',{},@(x)iscell(x)||isstring(x));
addParameter(p,'DataRoot','/home/blanche/data',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputDir',fullfile(baseCfg.ClassComparisonDir,'specified_features'),@(x)ischar(x)||isstring(x));
addParameter(p,'Visible',baseCfg.FigureVisible,@(x)ischar(x)||isstring(x));
addParameter(p,'NABins',21,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'Pseudocount',0.1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'AnatomyCoordinateMode','normalized',@(x)ischar(x)||isstring(x));
addParameter(p,'AnatomyGridSize',61,@(x)isnumeric(x)&&isscalar(x)&&x>=21);
addParameter(p,'AnatomySmoothBins',1.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MinCellsForDensity',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'RandomSeed',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(p,'MappingSampleFrames',500,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'MappingMinCorrelation',0.90,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'MappingMinMargin',0.01,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
addParameter(p,'RequireFovCorrection',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});

C = p.Results;
C.DataRoot = char(C.DataRoot);
C.OutputDir = char(C.OutputDir);
C.Visible = char(C.Visible);
C.AnatomyCoordinateMode = lower(string(C.AnatomyCoordinateMode));
assert(ismember(C.AnatomyCoordinateMode,["normalized","microns"]), ...
    'AnatomyCoordinateMode must be ''normalized'' or ''microns''.');

if ~isfolder(C.OutputDir), mkdir(C.OutputDir); end

C.MorphOrder = ["Surface","Molino","Pachon"];
C.ClassOrder = ["CW","CCW","Symmetric"];
C.SpeedOrder = ["positive","negative","independent"];
C.MorphColors = [.85 .20 .20; .20 .65 .30; .95 .72 .10];
C.ClassColors = [.25 .55 .85; .85 .35 .25; .55 .55 .55];
C.Alpha = .05;
C.MinFishPerMorphKW = 3;
C.PairwiseMorphComparisons = { ...
    'Surface','Molino'; ...
    'Surface','Pachon'; ...
    'Molino','Pachon'};
C.TurnThreshold = .239;
C.TauSeconds = 3;
C.KernelSeconds = 20;
C.AHVNativeToPaperSign = -1; % model convention: positive AHV is CCW
C.CandidateRasterMappingDir = baseCfg.CandidateRasterMappingDir;

if C.AnatomyCoordinateMode == "normalized"
    C.AnatomyXLabel = 'Aligned ML position (isotropic FOV units; right +)';
    C.AnatomyYLabel = 'Aligned RC position (isotropic FOV units; rostral +)';
else
    C.AnatomyXLabel = 'Aligned ML position (\mum; right +)';
    C.AnatomyYLabel = 'Aligned RC position (\mum; rostral +)';
end

files = cellstr(string(C.Files(:)));
if isempty(files)
    % Search only canonical morph directories; exclude archival roots.
    morphDirs = {'surface','molino','pachon'};
    for iMorph = 1:numel(morphDirs)
        d = dir(fullfile(C.DataRoot,morphDirs{iMorph},'**','*_candidate_neurons_clean_classified.mat'));
        files = [files; arrayfun(@(q)fullfile(q.folder,q.name),d,'UniformOutput',false)]; %#ok<AGROW>
    end
    if isempty(files)
        for iMorph = 1:numel(morphDirs)
            d = dir(fullfile(baseCfg.ClassifiedCandidateDir,morphDirs{iMorph},'**','*_candidate_neurons_clean_classified.mat'));
            files = [files; arrayfun(@(q)fullfile(q.folder,q.name),d,'UniformOutput',false)]; %#ok<AGROW>
        end
    end
end
assert(~isempty(files),'No classified candidate files found.');

F = repmat(emptyFish(),numel(files),1);
skippedFile=strings(0,1); skippedIdentifier=strings(0,1); skippedReason=strings(0,1);
for i = 1:numel(files)
    fprintf('Class-feature analysis %d/%d: %s\n',i,numel(files),files{i});
    try
        F(i) = loadFish(files{i},C);
    catch ME
        mappingFailure=any(contains(string(ME.message),[ ...
            "Candidate-to-original trace matching","Trace matching assigned", ...
            "Global one-to-one trace assignment","candidate-to-original mapping"]));
        if ~mappingFailure, rethrow(ME); end
        skippedFile(end+1,1)=string(files{i});
        skippedIdentifier(end+1,1)=string(ME.identifier);
        skippedReason(end+1,1)=string(ME.message);
        warning('Skipping recording from class-feature analysis because anatomical ROI mapping failed: %s\n%s', ...
            files{i},ME.message);
    end
end
skippedRecordings=table(skippedFile,skippedIdentifier,skippedReason, ...
    'VariableNames',{'File','ErrorIdentifier','Reason'});
F = F([F.nCells] > 0);
assert(~isempty(F),'No analyzable classified neurons.');

% Use one common signed-AHV grid. Extreme tails are clipped only for bin
% definition; every sample is assigned to the outermost applicable bin.
allA = vertcat(F.ahv);
lim = prctile(abs(allA(isfinite(allA))),99);
assert(isfinite(lim)&&lim>0,'No finite nonzero AHV samples.');
edges = linspace(-lim,lim,C.NABins+1);
centers = (edges(1:end-1)+edges(2:end))/2;
for i = 1:numel(F)
    F(i).ahvProfile = activityProfile(F(i),edges);
end

[fishSummary,cellTable,profileTable,speedTable,anatomyTable] = makeTables(F,C,centers);
candidateMappingTable=cellTable(:,{'Morph','Fish','CandidatePosition','CandidateID', ...
    'OriginalRoiID','AnatomyMatchCorrelation','AnatomyMatchMargin','Class'});
anatomyDensity = buildAnatomyDensity(F,C);
anatomyAudit = makeAnatomyAudit(F);

abundancePairedStats = abundancePairedTests(fishSummary,C);
[anatomySkewFish,anatomySkewStats] = computeAnatomySkew(cellTable,C);
rayleighFish = rayleighTable(F,C);
rayleighCombined = combineRayleigh(rayleighFish,C);
[ringRayleighFish,ringRayleighSummary,ringRayleighComparisons] = ...
    summarizeRingRayleigh(rayleighFish,C);
crossMorphKWStats = buildClassFeatureCrossMorphKW( ...
    fishSummary,anatomySkewFish,ringRayleighFish,C);

writetable(fishSummary,fullfile(C.OutputDir,'fish_summary.csv'));
writetable(skippedRecordings,fullfile(C.OutputDir,'skipped_recordings_anatomy_mapping.csv'));
writetable(cellTable,fullfile(C.OutputDir,'classified_cells.csv'));
writetable(candidateMappingTable,fullfile(C.OutputDir,'candidate_to_original_roi_mapping.csv'));
writetable(profileTable,fullfile(C.OutputDir,'fish_class_ahv_profiles.csv'));
writetable(speedTable,fullfile(C.OutputDir,'fish_direction_x_speed_proportions.csv'));
writetable(anatomyTable,fullfile(C.OutputDir,'fish_class_anatomy.csv'));
writetable(anatomyAudit,fullfile(C.OutputDir,'anatomy_alignment_audit.csv'));
writetable(abundancePairedStats,fullfile(C.OutputDir,'class_abundance_CW_vs_CCW_paired_tests.csv'));
writetable(anatomySkewFish,fullfile(C.OutputDir,'fish_CW_CCW_left_right_skewness.csv'));
writetable(anatomySkewStats,fullfile(C.OutputDir,'CW_CCW_left_right_skewness_tests.csv'));
writetable(rayleighFish,fullfile(C.OutputDir,'rayleigh_by_fish_and_class.csv'));
writetable(rayleighCombined,fullfile(C.OutputDir,'rayleigh_combined_pvalues.csv'));
writetable(ringRayleighFish,fullfile(C.OutputDir,'ring_rayleigh_fish_bonferroni.csv'));
writetable(ringRayleighSummary,fullfile(C.OutputDir,'ring_rayleigh_morph_summary.csv'));
writetable(ringRayleighComparisons,fullfile(C.OutputDir,'ring_rayleigh_cave_vs_surface_tests.csv'));
writetable(crossMorphKWStats,fullfile(C.OutputDir,'class_feature_cross_morph_pairwise_KW_tests.csv'));
save(fullfile(C.OutputDir,'fish_normalized_anatomy_density.mat'),'anatomyDensity','-v7.3');

figures = struct();
figures.abundance = plotAbundance(fishSummary,abundancePairedStats,crossMorphKWStats,C);
figures.anatomySkew = plotAnatomySkew(anatomySkewFish,anatomySkewStats,crossMorphKWStats,C);
figures.anatomyDensityDifference = plotAnatomyDensityDifference(anatomyDensity,C);
figures.anatomyQC = plotAnatomyQC(F,C);
figures.ahvProfiles = plotAHVProfiles(profileTable,C,centers);
figures.speed = plotSpeed(speedTable,C);
figures.preferredDirection = plotPreferred(cellTable,C);
figures.ringRayleigh = ...
    plotRingRayleigh(ringRayleighFish,rayleighCombined,crossMorphKWStats,C);

results = struct('config',C,'fish',fishSummary,'cells',cellTable, ...
    'candidateMapping',candidateMappingTable, ...
    'profiles',profileTable,'speed',speedTable, ...
    'anatomy',anatomyTable, ...
    'anatomyDensity',anatomyDensity,'abundancePairedStats',abundancePairedStats, ...
    'anatomySkewFish',anatomySkewFish,'anatomySkewStats',anatomySkewStats, ...
    'rayleighFish',rayleighFish,'rayleighCombined',rayleighCombined, ...
    'ringRayleighFish',ringRayleighFish,'ringRayleighSummary',ringRayleighSummary, ...
    'ringRayleighComparisons',ringRayleighComparisons, ...
    'crossMorphKWStats',crossMorphKWStats, ...
    'anatomyAudit',anatomyAudit,'skippedRecordings',skippedRecordings,'figures',figures);
save(fullfile(C.OutputDir,'class_feature_comparison.mat'),'results','-v7.3');
end

function E = emptyFish()
E = struct('morph',"",'fish',"",'file',"",'recordingFolder',"", ...
    'nCells',0,'ids',[],'candidatePositions',[],'originalIds',[], ...
    'mappingCorrelation',[],'mappingMargin',[],'mappingMethod',"", ...
    'labels',strings(0,1),'b1',[], ...
    'pB1Adjusted',[],'pref',[],'cvR2Full',[],'cvR2PhaseOnly',[],'cvR2AHVOnly',[],'x',[],'y',[],'xUm',[],'yUm',[], ...
    'xNorm',[],'yNorm',[],'traces',[],'ahv',[],'turnCW',NaN, ...
    'turnCCW',NaN,'ahvProfile',[],'anatomySource',"", ...
    'fovSource',"",'anatomy',emptyAnatomy());
end

function A = emptyAnatomy()
A = struct('available',false,'allCellsFile',"",'rasterFile',"",'fovFile',"", ...
    'candidateMappingMethod',"",'mappingCorrelation',[], ...
    'mappingSecondBestCorrelation',[],'mappingMargin',[], ...
    'rotationApplied',false,'fovAngleDeg',NaN,'R',eye(2), ...
    'linePos',[],'imageCenter',[NaN NaN],'imageSize',[NaN NaN], ...
    'pixelLengthX',NaN,'pixelLengthY',NaN,'background',[], ...
    'xyRawAll',[],'xUmAll',[],'yUmAll',[],'xNormAll',[], ...
    'yNormAll',[],'xAll',[],'yAll',[]);
end

function F = loadFish(file,C)
S = load(file);
[morph,fish] = pathIdentity(file,S);
F = emptyFish();
F.morph = morph;
F.fish = fish;
F.file = string(file);
F.recordingFolder = string(resolveRecordingFolder(file,S,morph,fish,C));

assert(isfield(S,'candidate_cell_ids')&&isfield(S,'neuron_class_label'), ...
    'Incomplete classified file: %s',file);

allIds = double(S.candidate_cell_ids(:));
allLabels = normalizeLabels(S.neuron_class_label);
assert(numel(allLabels)==numel(allIds), ...
    'candidate_cell_ids and neuron_class_label have different lengths: %s',file);

Y = double(S.calcium_traces);
if size(Y,2)~=numel(allIds) && size(Y,1)==numel(allIds), Y=Y'; end
assert(size(Y,2)==numel(allIds),'Trace/candidate mismatch: %s',file);

mapping = resolveCandidateToOriginalMapping( ...
    S,file,F.recordingFolder,allIds,Y,C);

keep = ismember(allLabels,C.ClassOrder);
F.ids = allIds(keep);
F.candidatePositions = find(keep);
F.originalIds = mapping.originalRoiIds(keep);
F.mappingCorrelation = mapping.correlation(keep);
F.mappingMargin = mapping.margin(keep);
F.mappingMethod = mapping.method;
F.labels = allLabels(keep);
F.nCells = nnz(keep);
if F.nCells == 0, return; end

F.traces = loadMappedDeltaFoF(mapping,keep);

F.b1 = nan(F.nCells,1);
F.pB1Adjusted = nan(F.nCells,1);
F.pref = nan(F.nCells,1);
F.cvR2Full = nan(F.nCells,1);
F.cvR2PhaseOnly = nan(F.nCells,1);
F.cvR2AHVOnly = nan(F.nCells,1);
if isfield(S,'three_model_fit_table') && istable(S.three_model_fit_table)
    T = S.three_model_fit_table;
    for j = 1:F.nCells
        k = find(double(T.neuronID)==F.ids(j),1);
        if isempty(k), continue; end
        F.b1(j) = double(T.b1(k));
        if ismember('pB1Adjusted',T.Properties.VariableNames)
            F.pB1Adjusted(j) = double(T.pB1Adjusted(k));
        end
        if ismember('cvR2Full',T.Properties.VariableNames), F.cvR2Full(j)=double(T.cvR2Full(k)); end
        if ismember('cvR2PhaseOnly',T.Properties.VariableNames), F.cvR2PhaseOnly(j)=double(T.cvR2PhaseOnly(k)); end
        if ismember('cvR2Behavior',T.Properties.VariableNames), F.cvR2AHVOnly(j)=double(T.cvR2Behavior(k)); end
        if ismember('sourcePrefRad',T.Properties.VariableNames)
            F.pref(j) = double(T.sourcePrefRad(k));
        elseif ismember('prefRad',T.Properties.VariableNames)
            F.pref(j) = double(T.prefRad(k));
        end
    end
end

[F.x,F.y,F.xUm,F.yUm,F.xNorm,F.yNorm,F.anatomy] = ...
    loadAnatomy(mapping,F.recordingFolder,F.originalIds,C);
F.anatomySource = F.anatomy.allCellsFile;
F.fovSource = F.anatomy.fovFile;

[F.ahv,F.turnCW,F.turnCCW] = loadAHV(S,file,F.recordingFolder,size(F.traces,1),C);
n = min(size(F.traces,1),numel(F.ahv));
F.traces = F.traces(1:n,:);
F.ahv = F.ahv(1:n);
end

function [morph,fish] = pathIdentity(file,S)
q = lower(string(file));
if isfield(S,'source_candidate_path')
    q = q + " " + lower(string(S.source_candidate_path));
end
if contains(q,'surface')
    morph = "Surface";
elseif contains(q,'molino')
    morph = "Molino";
elseif contains(q,'pachon')
    morph = "Pachon";
else
    error('Could not identify morph from classified file: %s',file);
end

t = regexp(char(q),'rec\d+','match','once','ignorecase');
if isempty(t), error('Could not identify rec number from: %s',file); end
fish = string(lower(t));
end

function recFolder = resolveRecordingFolder(classFile,S,morph,fish,C)
canonical = fullfile(C.DataRoot,lower(char(morph)),char(fish));
if isfolder(canonical)
    recFolder = canonical;
    return
end

possible = string.empty(0,1);
possible(end+1,1) = string(fileparts(classFile));
if isfield(S,'source_candidate_path')
    possible(end+1,1) = string(fileparts(char(string(S.source_candidate_path))));
end

for i = 1:numel(possible)
    folder = char(possible(i));
    while isfolder(folder)
        [parent,name] = fileparts(folder);
        if strcmpi(name,char(fish))
            recFolder = folder;
            warning('Canonical recording folder was absent; using %s.',recFolder);
            return
        end
        if isempty(parent) || strcmp(parent,folder), break; end
        folder = parent;
    end
end

error(['Recording folder not found. Expected: %s\n' ...
       'Set DataRoot or supply classified files inside their rec folders.'],canonical);
end

function L = normalizeLabels(x)
L = string(x(:));
L(ismember(lower(L),["ci","center","centre","symmetric"])) = "Symmetric";
L(strcmpi(L,'cw')) = "CW";
L(strcmpi(L,'ccw')) = "CCW";
end

function M = resolveCandidateToOriginalMapping(S,classFile,recFolder,candidateIds,candidateTraces,C)
sourceCandidatePath=classFile;
if isfield(S,'source_candidate_path')&&isfile(char(string(S.source_candidate_path)))
    sourceCandidatePath=char(string(S.source_candidate_path));
end
% Resolve candidate columns to the original RASTER/ALL_CELLS ROI index.
% candidate_cell_ids are never assumed to be original IDs. Every candidate
% trace is globally matched against every original RASTER column first.

allCellsFile=locateAllCellsFile(S,classFile,recFolder);
A=load(allCellsFile,'cells','cell_per');
if isfield(A,'cells')
    nAll=numel(A.cells);
elseif isfield(A,'cell_per')
    nAll=numel(A.cell_per);
else
    error('Cannot determine the ROI count in %s.',allCellsFile);
end

[M,cacheHit]=loadCandidateRasterMappingCache(sourceCandidatePath,recFolder,candidateIds,candidateTraces,allCellsFile,nAll,C);
if cacheHit, return; end

rasterFile=locateRasterFile(S,classFile,recFolder,allCellsFile,nAll);
W=whos('-file',rasterFile); rasterVariables=string({W.name});
if ismember("deltaFoF",rasterVariables)
    Q=load(rasterFile,'deltaFoF'); rasterTraces=double(Q.deltaFoF); traceField="deltaFoF";
elseif ismember("F0",rasterVariables)
    Q=load(rasterFile,'F0'); rasterTraces=double(Q.F0); traceField="F0";
else
    error('RASTER file has neither deltaFoF nor F0 for identity matching: %s',rasterFile);
end
if size(rasterTraces,2)~=nAll && size(rasterTraces,1)==nAll
    rasterTraces=rasterTraces';
end
assert(size(rasterTraces,2)==nAll, ...
    'RASTER and ALL_CELLS ROI counts differ (%d versus %d): %s', ...
    size(rasterTraces,2),nAll,rasterFile);
candidateFrameCount=size(candidateTraces,1);
rasterFrameCount=size(rasterTraces,1);
assert(candidateFrameCount<=rasterFrameCount, ...
    ['Candidate traces have more time frames than RASTER (%d versus %d), so ' ...
     'candidate-to-original mapping cannot be verified safely for %s.'], ...
    candidateFrameCount,rasterFrameCount,classFile);


nCandidates=numel(candidateIds);
assert(size(candidateTraces,2)==nCandidates, ...
    'Candidate trace count does not match candidate_cell_ids in %s.',classFile);
assert(nCandidates<=nAll,'More candidate traces than original ROIs in %s.',classFile);

% Candidate data are defined to overlap RASTER frames 1:nFrames. Match ROI
% columns only; never scan temporal offsets.
nFrames=candidateFrameCount;
firstWindow=rasterTraces(1:nFrames,:);
probeFrames=unique(round(linspace(1,nFrames,min(nFrames,round(C.MappingSampleFrames)))));
firstScore=corr(candidateTraces(probeFrames,:),firstWindow(probeFrames,:),'rows','pairwise');
firstScore(~isfinite(firstScore))=-Inf;
[firstSorted,firstIdx]=sort(firstScore,2,'descend');
if nAll>=2
    firstMargin=firstSorted(:,1)-firstSorted(:,2);
else
    firstMargin=Inf(nCandidates,1);
end
zeroOffsetProven=all(firstSorted(:,1)>=C.MappingMinCorrelation) && ...
    all(firstMargin>=C.MappingMinMargin) && ...
    numel(unique(firstIdx(:,1)))==nCandidates;
selectedOffset=0;
score=firstScore;
firstWindowScore=median(firstSorted(:,1),'omitnan');
rasterFrameStart=selectedOffset+1;
rasterFrameEnd=selectedOffset+nFrames;
rasterMatchTraces=rasterTraces(rasterFrameStart:rasterFrameEnd,:);
fprintf('  Fixed RASTER alignment: frames %d:%d, score=%.5f\n', ...
    rasterFrameStart,rasterFrameEnd,firstWindowScore);
[sortedScore,sortedIdx]=sort(score,2,'descend');
originalIds=double(sortedIdx(:,1));
if numel(unique(originalIds))<nCandidates && exist('matchpairs','file')==2
    pairs=matchpairs(-score,1e6,'min');
    assert(size(pairs,1)==nCandidates,'Global one-to-one trace assignment is incomplete.');
    [~,pairOrder]=sort(pairs(:,1));
    originalIds=double(pairs(pairOrder,2));
end
probeBest=score(sub2ind(size(score),(1:nCandidates)',originalIds));
fullBest=nan(nCandidates,1);
for candidateIndex=1:nCandidates
    q=corr(candidateTraces(:,candidateIndex),rasterMatchTraces(:,originalIds(candidateIndex)),'rows','pairwise');
    fullBest(candidateIndex)=q;
end
secondBest=nan(nCandidates,1);
for candidateIndex=1:nCandidates
    alternatives=score(candidateIndex,:);
    alternatives(originalIds(candidateIndex))=-Inf;
    secondBest(candidateIndex)=max(alternatives);
end
margin=probeBest-secondBest;

assert(numel(unique(originalIds))==nCandidates, ...
    ['Trace matching assigned the same original ROI to multiple candidates in %s. ' ...
     'No anatomical mapping was accepted.'],classFile);

bestCorr=fullBest;
badCorrelation=~isfinite(bestCorr)|bestCorr<C.MappingMinCorrelation;
badMargin=~isfinite(margin)|margin<C.MappingMinMargin;
if any(badCorrelation|badMargin)
    bad=find(badCorrelation|badMargin);
    error(['Candidate-to-original trace matching was ambiguous for %d/%d candidates ' ...
        'in %s. Worst full correlation = %.4f; smallest full-frame margin = %.4f. ' ...
        'No anatomical results were produced for this fish.'], ...
        numel(bad),nCandidates,classFile,min(bestCorr),min(margin));
end

candidateIdsCanIndexOriginal=all(isfinite(candidateIds)&candidateIds==round(candidateIds)& ...
    candidateIds>=1&candidateIds<=nAll) && numel(unique(candidateIds))==nCandidates;
if candidateIdsCanIndexOriginal && all(originalIds==candidateIds)
    method="candidate_ID_is_original_ROI_global_trace_verified";
else
    method="candidate_trace_matched_to_original_RASTER";
end

assert(all(originalIds>=1&originalIds<=nAll&originalIds==round(originalIds)), ...
    'Resolved original ROI IDs are invalid for %s.',allCellsFile);

M=struct();
M.candidatePositions=(1:nCandidates)';
M.candidateIds=candidateIds(:);
M.originalRoiIds=originalIds(:);
M.correlation=bestCorr(:);
M.secondBestCorrelation=secondBest(:);
M.margin=margin(:);
M.method=method;
M.traceField=traceField;
M.allCellsFile=string(allCellsFile);
M.rasterFile=string(rasterFile);
M.nOriginalRois=nAll;
M.candidateFrameCount=candidateFrameCount;
M.rasterFrameCount=rasterFrameCount;
M.rasterFrameStart=rasterFrameStart;
M.rasterFrameEnd=rasterFrameEnd;
M.selectedRasterOffset=selectedOffset;
M.selectedWindowScore=firstWindowScore;
M.zeroOffsetProven=zeroOffsetProven;
saveCandidateRasterMappingCache(M,sourceCandidatePath,candidateIds,candidateTraces,C);

fprintf(['  Candidate anatomy mapping: %s\n' ...
    '  RASTER field: %s | frames %d:%d of %d | median r = %.5f | minimum r = %.5f\n'], ...
    method,traceField,rasterFrameStart,rasterFrameEnd,rasterFrameCount, ...
    median(bestCorr,'omitnan'),min(bestCorr));
end

function traces=loadMappedDeltaFoF(mapping,keep)
% Retrieve the original, aligned RASTER traces in physical deltaF/F0 units.
W=whos("-file",char(mapping.rasterFile));
assert(ismember("deltaFoF",string({W.name})), ...
    "RASTER file lacks deltaFoF; cannot plot activity in deltaF/F0: %s", ...
    mapping.rasterFile);
Q=load(char(mapping.rasterFile),"deltaFoF"); raw=double(Q.deltaFoF);
if size(raw,2)~=mapping.nOriginalRois && size(raw,1)==mapping.nOriginalRois
    raw=transpose(raw);
end
rows=mapping.rasterFrameStart:mapping.rasterFrameEnd;
ids=mapping.originalRoiIds(keep);
traces=raw(rows,ids);
end

function [M,hit]=loadCandidateRasterMappingCache(candidatePath,recFolder,candidateIds,candidateTraces,allCellsFile,nAll,C)
M=struct(); hit=false; [cacheFile,candidateInfo]=candidateRasterCachePath(candidatePath,recFolder,C.CandidateRasterMappingDir);
if ~isfile(cacheFile), return; end
Q=load(cacheFile,"candidateRasterMapping"); if ~isfield(Q,"candidateRasterMapping"), return; end
D=Q.candidateRasterMapping; rasterInfo=dir(char(D.rasterFile));
valid=isfield(D,"version")&&D.version==2&&D.rawFrameStart==1&&strcmp(char(D.candidatePath),candidatePath)&& ...
    D.candidateBytes==candidateInfo.bytes&&D.candidateDatenum==candidateInfo.datenum&& ...
    isequal(D.candidateTraceSize,size(candidateTraces))&&isequal(D.candidateIds(:),candidateIds(:))&& ...
    ~isempty(rasterInfo)&&D.rasterBytes==rasterInfo.bytes&&D.rasterDatenum==rasterInfo.datenum&& ...
    numel(D.originalRoiIds)==numel(candidateIds)&&all(D.originalRoiIds>=1&D.originalRoiIds<=nAll);
if ~valid, warning("Ignoring stale candidate/RASTER mapping cache: %s",cacheFile); return; end
M=struct("candidatePositions",transpose(1:numel(candidateIds)),"candidateIds",candidateIds(:), ...
    "originalRoiIds",D.originalRoiIds(:),"correlation",D.correlation(:), ...
    "secondBestCorrelation",D.correlation(:)-D.margin(:),"margin",D.margin(:), ...
    "method","cached_candidate_trace_matched_to_original_RASTER","traceField","deltaFoF", ...
    "allCellsFile",string(allCellsFile),"rasterFile",string(D.rasterFile), ...
    "nOriginalRois",nAll,"candidateFrameCount",size(candidateTraces,1), ...
    "rasterFrameCount",NaN,"rasterFrameStart",D.rawFrameStart,"rasterFrameEnd",D.rawFrameEnd, ...
    "selectedRasterOffset",D.rawFrameStart-1,"selectedWindowScore",median(D.correlation,"omitnan"), ...
    "zeroOffsetProven",D.rawFrameStart==1);
hit=true; fprintf("  Reused candidate/RASTER mapping cache: %s\n",cacheFile);
end

function saveCandidateRasterMappingCache(M,candidatePath,candidateIds,candidateTraces,C)
[cacheFile,candidateInfo]=candidateRasterCachePath(candidatePath,fileparts(candidatePath),C.CandidateRasterMappingDir);
rasterInfo=dir(char(M.rasterFile));
candidateRasterMapping=struct("version",2,"candidatePath",string(candidatePath), ...
    "candidateBytes",candidateInfo.bytes,"candidateDatenum",candidateInfo.datenum, ...
    "candidateTraceSize",size(candidateTraces),"candidateIds",candidateIds(:), ...
    "rasterFile",string(M.rasterFile),"rasterBytes",rasterInfo.bytes,"rasterDatenum",rasterInfo.datenum, ...
    "rawFrameStart",M.rasterFrameStart,"rawFrameEnd",M.rasterFrameEnd, ...
    "originalRoiIds",M.originalRoiIds(:),"correlation",M.correlation(:),"margin",M.margin(:));
folder=fileparts(cacheFile); if ~isfolder(folder), mkdir(folder); end
save(cacheFile,"candidateRasterMapping","-v7"); fprintf("  Saved candidate/RASTER mapping cache: %s\n",cacheFile);
end

function [file,info]=candidateRasterCachePath(candidatePath,recFolder,cacheRoot)
[parent,recording]=fileparts(recFolder); [~,morph]=fileparts(parent); [~,candidateName]=fileparts(candidatePath);
file=fullfile(cacheRoot,morph,recording,candidateName + "_raster_mapping.mat");
info=dir(candidatePath); assert(~isempty(info),"Candidate file not found: %s",candidatePath); info=info(1);
end

function X = normalizeColumnsForCorrelation(X)
mu=mean(X,1,'omitnan'); X=X-mu; X(~isfinite(X))=0;
ss=sqrt(sum(X.^2,1));
bad=~isfinite(ss)|ss<=eps;
ss(bad)=1; X=X./ss; X(:,bad)=0;
end

function r = pairedColumnCorrelations(A,B)
assert(isequal(size(A),size(B)),'Paired correlation matrices must have equal size.');
A=normalizeColumnsForCorrelation(A); B=normalizeColumnsForCorrelation(B);
r=sum(A.*B,1)';
end

function [x,y,xUm,yUm,xNorm,yNorm,Aout] = loadAnatomy(mapping,recFolder,ids,C)
Aout = emptyAnatomy();
allCellsFile = char(mapping.allCellsFile);
fovFile = locateFovFile(recFolder,C.RequireFovCorrection);

A = load(allCellsFile,'cells','cell_per','pixelLengthX','pixelLengthY','avg','bkg');
if isfield(A,'avg')
    background = squeeze(double(A.avg));
elseif isfield(A,'bkg')
    background = squeeze(double(A.bkg));
else
    error('ALL_CELLS file has no avg or bkg image: %s',allCellsFile);
end
if ndims(background)>2, background=prctile(background,90,3); end
imgSize = size(background);
imgSize = imgSize(1:2);

if isfield(A,'cells') && iscell(A.cells)
    xyRaw = centroidsFromLinearPixels(A.cells,imgSize);
elseif isfield(A,'cell_per') && iscell(A.cell_per)
    xyRaw = centroidsFromContours(A.cell_per);
else
    error('ALL_CELLS file contains neither usable cells nor cell_per: %s',allCellsFile);
end

nAll = size(xyRaw,1);
assert(all(ids==round(ids) & ids>=1 & ids<=nAll), ...
    'Resolved original ROI IDs do not index the ALL_CELLS file: %s',allCellsFile);

sx = 1; sy = 1;
if isfield(A,'pixelLengthX') && isfinite(double(A.pixelLengthX))
    sx = double(A.pixelLengthX);
end
if isfield(A,'pixelLengthY') && isfinite(double(A.pixelLengthY))
    sy = double(A.pixelLengthY);
end

R = eye(2); fovAngleDeg = NaN; linePos = []; rotationApplied = false;
if isfile(fovFile)
    Q = load(fovFile);
    if isfield(Q,'fov_correction_saved'), Q=Q.fov_correction_saved; end
    assert(isfield(Q,'R_fov') && isequal(size(Q.R_fov),[2 2]), ...
        'fov_correction file has no 2x2 R_fov: %s',fovFile);
    R = double(Q.R_fov);
    assert(all(isfinite(R),'all') && abs(det(R)-1)<1e-3 && ...
        norm(R*R'-eye(2),'fro')<1e-3, ...
        'R_fov is not a valid rotation matrix: %s',fovFile);
    if isfield(Q,'fov_angle_deg'), fovAngleDeg=double(Q.fov_angle_deg); end
    if isfield(Q,'line_pos') && isequal(size(Q.line_pos),[2 2])
        linePos = double(Q.line_pos);
    end
    rotationApplied = true;
end

imageCenter = [(imgSize(2)+1)/2, (imgSize(1)+1)/2];
xyCenteredPixels = xyRaw-imageCenter;
xyAlignedPixels = (R*xyCenteredPixels')';

% With isotropic pixels, use the saved pixel-space rotation directly. If
% pixels are anisotropic and the manual axis is available, reconstruct the
% rotation in physical coordinates from that caudal-to-rostral line.
relativePixelMismatch = abs(sx-sy)/mean([sx sy]);
if relativePixelMismatch <= .01
    xyAlignedUm = xyAlignedPixels*mean([sx sy]);
elseif ~isempty(linePos)
    xyPhysical = xyCenteredPixels.*[sx sy];
    v = (linePos(2,:)-linePos(1,:)).*[sx sy];
    alpha = -pi/2-atan2(v(2),v(1));
    Rphysical = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
    xyAlignedUm = (Rphysical*xyPhysical')';
    warning('Anisotropic pixels in %s: rotation reconstructed from line_pos.',allCellsFile);
else
    xyPhysical = xyCenteredPixels.*[sx sy];
    xyAlignedUm = (R*xyPhysical')';
    warning('Anisotropic pixels but no line_pos in %s; using R_fov approximately.',fovFile);
end

% Convert image y-down to an anatomical RC axis with rostral positive.
xUmAll = xyAlignedUm(:,1);
yUmAll = -xyAlignedUm(:,2);
isotropicHalfFovUm = max([imgSize(2)*sx,imgSize(1)*sy])/2;
xNormAll = xUmAll/isotropicHalfFovUm;
yNormAll = yUmAll/isotropicHalfFovUm;

if C.AnatomyCoordinateMode == "normalized"
    xAll = xNormAll; yAll = yNormAll;
else
    xAll = xUmAll; yAll = yUmAll;
end

id = ids(:);
x = xAll(id); y = yAll(id);
xUm = xUmAll(id); yUm = yUmAll(id);
xNorm = xNormAll(id); yNorm = yNormAll(id);

Aout.available = true;
Aout.allCellsFile = string(allCellsFile);
Aout.rasterFile = mapping.rasterFile;
Aout.fovFile = string(fovFile);
Aout.candidateMappingMethod = mapping.method;
Aout.mappingCorrelation = mapping.correlation;
Aout.mappingSecondBestCorrelation = mapping.secondBestCorrelation;
Aout.mappingMargin = mapping.margin;
Aout.rotationApplied = rotationApplied;
Aout.fovAngleDeg = fovAngleDeg;
Aout.R = R;
Aout.linePos = linePos;
Aout.imageCenter = imageCenter;
Aout.imageSize = imgSize;
Aout.pixelLengthX = sx;
Aout.pixelLengthY = sy;
Aout.background = background;
Aout.xyRawAll = xyRaw;
Aout.xUmAll = xUmAll;
Aout.yUmAll = yUmAll;
Aout.xNormAll = xNormAll;
Aout.yNormAll = yNormAll;
Aout.xAll = xAll;
Aout.yAll = yAll;
end

function file = locateAllCellsFile(S,classFile,recFolder)
candidateFields = {'source_all_cells_path','all_cells_path','source_all_cells_file'};
for i = 1:numel(candidateFields)
    if isfield(S,candidateFields{i})
        q = char(string(S.(candidateFields{i})));
        if isfile(q), file=promoteToOriginalAllCells(q); return; end
        q2 = fullfile(recFolder,q);
        if isfile(q2), file=promoteToOriginalAllCells(q2); return; end
    end
end

d = dir(fullfile(recFolder,'*ALL_CELLS*.mat'));
d = d(~[d.isdir]);
% Candidate/subset ALL_CELLS files use their own local index space. The
% anatomical reference must be the original file paired with the original
% RASTER, so exclude subset products from automatic discovery.
namesLower=lower(string({d.name}));
isSubset=contains(namesLower,'selected_rois')|contains(namesLower,'non_selected_rois');
d=d(~isSubset);
if isempty(d)
    error('No original (non-subset) ALL_CELLS.mat file found in %s.',recFolder);
end
if numel(d)==1
    file = fullfile(d(1).folder,d(1).name);
    return
end

[~,className] = fileparts(classFile);
prefix = regexprep(className,'_candidate_neurons.*$','','ignorecase');
match = startsWith(string({d.name}),string(prefix),'IgnoreCase',true);
if nnz(match)==1
    k=find(match,1); file=fullfile(d(k).folder,d(k).name); return
end

names = strjoin(string({d.name}),', ');
error(['Several ALL_CELLS files were found in %s and could not be matched ' ...
       'unambiguously to %s: %s'],recFolder,className,names);
end

function file = promoteToOriginalAllCells(file)
nameLower=lower(string(file));
if ~(contains(nameLower,'selected_rois')||contains(nameLower,'non_selected_rois'))
    return
end
Q=load(file,'originalAllCellsFile');
if isfield(Q,'originalAllCellsFile') && isfile(char(string(Q.originalAllCellsFile)))
    file=char(string(Q.originalAllCellsFile));
else
    error(['An ALL_CELLS subset was referenced but its originalAllCellsFile ' ...
        'could not be resolved: %s'],file);
end
end

function file = locateRasterFile(S,classFile,recFolder,allCellsFile,nAll)
candidateFields={'source_raster_path','raster_path','source_raster_file'};
for i=1:numel(candidateFields)
    if isfield(S,candidateFields{i})
        q=char(string(S.(candidateFields{i})));
        if isfile(q) && rasterHasRoiCount(q,nAll), file=q; return; end
        q2=fullfile(recFolder,q);
        if isfile(q2) && rasterHasRoiCount(q2,nAll), file=q2; return; end
    end
end

[folder,allName]=fileparts(allCellsFile);
pairedName=regexprep(allName,'_ALL_CELLS$','_RASTER','ignorecase');
pairedFile=fullfile(folder,[pairedName '.mat']);
if isfile(pairedFile) && rasterHasRoiCount(pairedFile,nAll)
    file=pairedFile;
    return
end

d=dir(fullfile(recFolder,'*RASTER*.mat')); d=d(~[d.isdir]);
namesLower=lower(string({d.name}));
isSubset=contains(namesLower,'selected_rois')|contains(namesLower,'non_selected_rois');
d=d(~isSubset);

% Retain only RASTER files whose ROI count equals the original ALL_CELLS.
keep=false(numel(d),1);
for i=1:numel(d)
    q=fullfile(d(i).folder,d(i).name);
    keep(i)=rasterHasRoiCount(q,nAll);
end
d=d(keep);

if numel(d)==1
    file=fullfile(d(1).folder,d(1).name);
    return
end
if isempty(d)
    error('No original RASTER with %d ROIs was found in %s.',nAll,recFolder);
end
names=strjoin(string({d.name}),', ');
error(['Several original RASTER files with %d ROIs were found in %s and ' ...
    'could not be matched to %s: %s'],nAll,recFolder,classFile,names);
end

function tf = rasterHasRoiCount(file,nAll)
tf=false;
W=whos('-file',file); variableNames=string({W.name});
if ismember("deltaFoF",variableNames)
    info=W(variableNames=="deltaFoF"); dims=info.size;
elseif ismember("raster",variableNames)
    info=W(variableNames=="raster"); dims=info.size;
else
    return
end
tf=any(dims==nAll);
end

function file = locateFovFile(recFolder,required)
file = fullfile(recFolder,'fov_correction.mat');
if isfile(file), return; end
d = dir(fullfile(recFolder,'fov_correction*.mat'));
d = d(~[d.isdir]);
if numel(d)==1
    file=fullfile(d(1).folder,d(1).name);
elseif required
    error('Exactly one fov_correction.mat is required in %s.',recFolder);
else
    file='';
    warning('No FOV correction found in %s; using identity rotation.',recFolder);
end
end

function xy = centroidsFromLinearPixels(cells,imgSize)
n = numel(cells); xy=nan(n,2);
for i = 1:n
    idx=double(cells{i}(:));
    idx=idx(isfinite(idx) & idx>=1 & idx<=prod(imgSize));
    if isempty(idx), continue; end
    [yy,xx]=ind2sub(imgSize,round(idx));
    xy(i,:)=[mean(xx,'omitnan'),mean(yy,'omitnan')];
end
end

function xy = centroidsFromContours(contours)
n=numel(contours); xy=nan(n,2);
for i=1:n
    q=double(contours{i});
    if isempty(q), continue; end
    if size(q,2)<2 && size(q,1)==2, q=q'; end
    if size(q,2)>=2
        xy(i,:)=[mean(q(:,1),'omitnan'),mean(q(:,2),'omitnan')];
    end
end
end

function [ahv,nCW,nCCW] = loadAHV(S,classFile,recFolder,nTime,C)
ahv=nan(nTime,1); nCW=NaN; nCCW=NaN; bf='';
if isfield(S,"source_behavior_path") && isfile(char(string(S.source_behavior_path))), bf=char(string(S.source_behavior_path)); end
if isempty(bf)
 d=dir(fullfile(char(recFolder),"*swimResults_pass2*.mat"));
 if isempty(d), d=dir(fullfile(fileparts(classFile),"*swimResults_pass2*.mat")); end
 if ~isempty(d), bf=fullfile(d(1).folder,d(1).name); end
end
if isempty(bf) || ~isfile(bf), return; end
Q=load(bf,"pass2FileResult"); if ~isfield(Q,"pass2FileResult"), return; end
B=Q.pass2FileResult;
if isfield(B,"dtheta"), bias=double(B.dtheta(:)); elseif isfield(B,"LI"), bias=double(B.LI(:)); else, return; end
if max(abs(bias),[],"omitnan")>2*pi+.5, bias=deg2rad(bias); end
bias=mod(bias+pi,2*pi)-pi; sf=double(B.swimFrames);
if size(sf,1)==2, start=transpose(sf(1,:)); elseif size(sf,2)==2, start=sf(:,1); else, error("Invalid swimFrames shape."); end
n=min(numel(start),numel(bias)); start=round(start(1:n)); bias=bias(1:n);
qual=isfinite(start)&isfinite(bias)&abs(bias)>C.TurnThreshold;
nCW=nnz(qual&bias>0); nCCW=nnz(qual&bias<0); fps=double(B.fps(1));
if isfield(S,"time_s")&&numel(S.time_s)>=nTime, t=double(S.time_s(1:nTime)); t=t(:); else, t=transpose(0:nTime-1)/3.41; end
nBeh=inferBehaviorLengthForAHV(B); assert(nBeh>=2,"Could not infer behavior length.");
tb=chooseBehaviorTimeForAHV(B,nBeh,fps,t); validFrame=qual&start>=1&start<=nBeh;
boutTime=nan(n,1); boutTime(validFrame)=tb(start(validFrame));
idx=round(interp1(t,transpose(1:nTime),boutTime,"linear",NaN)); ok=qual&isfinite(idx)&idx>=1&idx<=nTime;
r=accumarray(idx(ok),C.AHVNativeToPaperSign*bias(ok),[nTime 1],@sum,0);
dt=median(diff(t)); kernelTime=transpose(0:floor(C.KernelSeconds/dt))*dt;
k=exp(-kernelTime/C.TauSeconds); k=k/sum(k); ahv=conv(r,k); ahv=ahv(1:nTime);
% Exact fitted predictor, scaled from rad/sample to deg/s.
ahv=rad2deg(ahv/dt); finiteTb=tb(isfinite(tb)); ahv(t<min(finiteTb)|t>max(finiteTb))=NaN;
end

function n=inferBehaviorLengthForAHV(B)
n=NaN; names={"behaviorTimeSec","heading_est","tail_angle"};
for i=1:numel(names), if isfield(B,names{i})&&~isempty(B.(names{i})), n=numel(B.(names{i})); return; end, end
if isfield(B,"nFrames")&&isfinite(double(B.nFrames(1))), n=round(double(B.nFrames(1))); end
end

function t=chooseBehaviorTimeForAHV(B,n,fps,tCa)
local=transpose(0:n-1)/fps; startFrame=1;
if isfield(B,"startFrame")&&isfinite(double(B.startFrame(1))), startFrame=round(double(B.startFrame(1))); end
fromStart=(transpose(startFrame:startFrame+n-1)-1)/fps; candidates={local,fromStart};
if isfield(B,"behaviorTimeSec")&&numel(B.behaviorTimeSec)==n, candidates{end+1}=double(B.behaviorTimeSec(:)); end
best=1; bestCoverage=-Inf; bestStart=Inf;
for i=1:numel(candidates)
 q=candidates{i}; fq=q(isfinite(q)); if numel(fq)<2, continue; end
 coverage=nnz(isfinite(tCa)&tCa>=min(fq)&tCa<=max(fq)); startDistance=abs(min(tCa(isfinite(tCa)))-min(fq));
 if coverage>bestCoverage||(coverage==bestCoverage&&startDistance<bestStart), best=i; bestCoverage=coverage; bestStart=startDistance; end
end
t=candidates{best};
end

function P = activityProfile(F,edges)
nb=numel(edges)-1; P=nan(3,nb);
bin=discretize(F.ahv,[-inf edges(2:end-1) inf]);
classOrder=["CW","CCW","Symmetric"];
for c=1:3
    q=F.labels==classOrder(c);
    for b=1:nb
        P(c,b)=mean(F.traces(bin==b,q),'all','omitnan');
    end
end
end

function [S,Cel,Prof,Speed,Ana] = makeTables(F,C,centers)
S=table(); Cel=table(); Prof=table(); Speed=table(); Ana=table();
for i=1:numel(F)
    n=F(i).nCells;
    nc=[nnz(F(i).labels=="CW"),nnz(F(i).labels=="CCW"),nnz(F(i).labels=="Symmetric")];
    pr=nc/n;
    lr=log((pr(2)+C.Pseudocount/n)/(pr(1)+C.Pseudocount/n));
    tr=log((F(i).turnCCW+C.Pseudocount)/(F(i).turnCW+C.Pseudocount));
    S=[S;table(F(i).morph,F(i).fish,n,nc(1),nc(2),nc(3), ...
        pr(1),pr(2),pr(3),lr,F(i).turnCW,F(i).turnCCW,tr, ...
        F(i).recordingFolder,F(i).anatomySource,F(i).anatomy.rasterFile, ...
        F(i).fovSource,F(i).mappingMethod, ...
        'VariableNames',{'Morph','Fish','NAll','NCW','NCCW','NSymmetric', ...
        'PropCW','PropCCW','PropSymmetric','CellLogCCWOverCW','NTurnCW', ...
        'NTurnCCW','TurnLogCCWOverCW','RecordingFolder','AnatomySource', ...
        'RasterSource','FovSource','CandidateMappingMethod'})];

    speed=repmat("independent",n,1);
    sig=F(i).pB1Adjusted<C.Alpha;
    speed(sig&F(i).b1>0)="positive";
    speed(sig&F(i).b1<0)="negative";

    Cel=[Cel;table(repmat(F(i).morph,n,1),repmat(F(i).fish,n,1), ...
        F(i).candidatePositions,F(i).ids,F(i).ids,F(i).originalIds, ...
        F(i).mappingCorrelation,F(i).mappingMargin,F(i).labels, ...
        F(i).b1,F(i).pB1Adjusted,speed,F(i).pref,F(i).cvR2Full,F(i).cvR2PhaseOnly, ...
        F(i).cvR2AHVOnly,F(i).cvR2Full-F(i).cvR2AHVOnly, ...
        F(i).cvR2Full-F(i).cvR2PhaseOnly, ...
        F(i).x,F(i).y,F(i).xUm,F(i).yUm,F(i).xNorm,F(i).yNorm, ...
        'VariableNames',{'Morph','Fish','CandidatePosition','CandidateID', ...
        'NeuronID','OriginalRoiID','AnatomyMatchCorrelation','AnatomyMatchMargin', ...
        'Class','Beta1','Beta1PAdjusted','SpeedClass','PreferredRad', ...
        'CvR2Full','CvR2PhaseOnly','CvR2AHVOnly','DeltaR2Phase','DeltaR2AHV','X','Y', ...
        'MLum','RCum','MLnorm','RCnorm'})];

    for c=1:3
        q=F(i).labels==C.ClassOrder(c);
        anatomical=q&isfinite(F(i).x)&isfinite(F(i).y);
        nRight=nnz(anatomical&F(i).x>0);
        nLeft=nnz(anatomical&F(i).x<0);
        if nRight+nLeft>0
            laterality=(nRight-nLeft)/(nRight+nLeft);
        else
            laterality=NaN;
        end
        nTopLeft=nnz(anatomical&F(i).x<0&F(i).y>0);
        nTopRight=nnz(anatomical&F(i).x>0&F(i).y>0);
        nBottomLeft=nnz(anatomical&F(i).x<0&F(i).y<0);
        nBottomRight=nnz(anatomical&F(i).x>0&F(i).y<0);

        Ana=[Ana;table(F(i).morph,F(i).fish,C.ClassOrder(c),nnz(q), ...
            nnz(anatomical),median(F(i).x(q),'omitnan'),median(F(i).y(q),'omitnan'), ...
            median(F(i).xUm(q),'omitnan'),median(F(i).yUm(q),'omitnan'), ...
            median(F(i).xNorm(q),'omitnan'),median(F(i).yNorm(q),'omitnan'), ...
            nRight,nLeft,laterality,nTopLeft,nTopRight,nBottomLeft,nBottomRight, ...
            F(i).anatomySource,F(i).fovSource, ...
            'VariableNames',{'Morph','Fish','Class','N','NWithAnatomy', ...
            'MedianX','MedianY','MedianMLum','MedianRCum','MedianMLnorm', ...
            'MedianRCnorm','NRight','NLeft','LateralityIndex', ...
            'NTopLeft','NTopRight','NBottomLeft','NBottomRight', ...
            'AnatomySource','FovSource'})];

        for b=1:numel(centers)
            Prof=[Prof;table(F(i).morph,F(i).fish,C.ClassOrder(c), ...
                centers(b),F(i).ahvProfile(c,b),'VariableNames', ...
                {'Morph','Fish','Class','AHVBin','DeltaFoF'})];
        end
        for s=1:3
            Speed=[Speed;table(F(i).morph,F(i).fish,C.ClassOrder(c), ...
                C.SpeedOrder(s),nnz(q&speed==C.SpeedOrder(s))/max(nnz(q),1), ...
                'VariableNames',{'Morph','Fish','DirectionClass','SpeedClass', ...
                'ProportionAllCells'})];
        end
    end
end
end

function D = buildAnatomyDensity(F,C)
allX=vertcat(F.x); allY=vertcat(F.y);
allX=allX(isfinite(allX)); allY=allY(isfinite(allY));
assert(~isempty(allX)&&~isempty(allY),'No finite aligned anatomical coordinates.');
xMax=max(abs(allX)); yMax=max(abs(allY));
if xMax==0, xMax=1; end
if yMax==0, yMax=1; end
xMax=1.04*xMax; yMax=1.04*yMax;

n=C.AnatomyGridSize;
xEdges=linspace(-xMax,xMax,n+1); yEdges=linspace(-yMax,yMax,n+1);
xCenters=(xEdges(1:end-1)+xEdges(2:end))/2;
yCenters=(yEdges(1:end-1)+yEdges(2:end))/2;

meanDensity=nan(n,n,3,3); nFish=zeros(3,3);
fishDensity=cell(3,3); differenceMean=nan(n,n,3); differenceNFish=zeros(3,1);

for m=1:3
    fishIdx=find([F.morph]==C.MorphOrder(m));
    for c=1:3
        M=[];
        for ii=fishIdx
            q=F(ii).labels==C.ClassOrder(c)&isfinite(F(ii).x)&isfinite(F(ii).y);
            if nnz(q)<C.MinCellsForDensity, continue; end
            H=oneFishDensity(F(ii).x(q),F(ii).y(q),xEdges,yEdges,C.AnatomySmoothBins);
            M=cat(3,M,H);
        end
        fishDensity{m,c}=M;
        nFish(m,c)=size(M,3);
        if ~isempty(M), meanDensity(:,:,m,c)=mean(M,3,'omitnan'); end
    end

    MDiff=[];
    for ii=fishIdx
        qCW=F(ii).labels=="CW"&isfinite(F(ii).x)&isfinite(F(ii).y);
        qCCW=F(ii).labels=="CCW"&isfinite(F(ii).x)&isfinite(F(ii).y);
        if nnz(qCW)<C.MinCellsForDensity || nnz(qCCW)<C.MinCellsForDensity, continue; end
        hCW=oneFishDensity(F(ii).x(qCW),F(ii).y(qCW),xEdges,yEdges,C.AnatomySmoothBins);
        hCCW=oneFishDensity(F(ii).x(qCCW),F(ii).y(qCCW),xEdges,yEdges,C.AnatomySmoothBins);
        MDiff=cat(3,MDiff,hCW-hCCW);
    end
    differenceNFish(m)=size(MDiff,3);
    if ~isempty(MDiff), differenceMean(:,:,m)=mean(MDiff,3,'omitnan'); end
end

D=struct('xEdges',xEdges,'yEdges',yEdges,'xCenters',xCenters, ...
    'yCenters',yCenters,'xLim',[-xMax xMax],'yLim',[-yMax yMax], ...
    'meanDensity',meanDensity,'fishDensity',{fishDensity},'nFish',nFish, ...
    'differenceMean',differenceMean,'differenceNFish',differenceNFish, ...
    'coordinateMode',C.AnatomyCoordinateMode,'minimumCells',C.MinCellsForDensity);
end

function H = oneFishDensity(x,y,xEdges,yEdges,sigmaBins)
H=histcounts2(y,x,yEdges,xEdges);
if sigmaBins>0
    radius=max(1,ceil(3*sigmaBins));
    [gx,gy]=meshgrid(-radius:radius,-radius:radius);
    K=exp(-(gx.^2+gy.^2)/(2*sigmaBins^2)); K=K/sum(K,'all');
    H=conv2(H,K,'same');
end
total=sum(H,'all');
if total>0, H=H/total; else, H(:)=NaN; end
end

function T = makeAnatomyAudit(F)
T=table();
for i=1:numel(F)
    A=F(i).anatomy;
    T=[T;table(F(i).morph,F(i).fish,F(i).recordingFolder, ...
        A.allCellsFile,A.rasterFile,A.fovFile,F(i).mappingMethod, ...
        median(F(i).mappingCorrelation,'omitnan'),min(F(i).mappingCorrelation), ...
        median(F(i).mappingMargin,'omitnan'),numel(unique(F(i).originalIds)), ...
        all(F(i).ids==F(i).originalIds),A.rotationApplied,A.fovAngleDeg, ...
        A.imageSize(2),A.imageSize(1),A.pixelLengthX,A.pixelLengthY, ...
        F(i).nCells,nnz(isfinite(F(i).x)&isfinite(F(i).y)), ...
        'VariableNames',{'Morph','Fish','RecordingFolder','AllCellsFile', ...
        'RasterFile','FovFile','CandidateMappingMethod','MedianMatchCorrelation', ...
        'MinimumMatchCorrelation','MedianMatchMargin','NUniqueOriginalRois', ...
        'CandidateIDsEqualOriginalIDs','RotationApplied','FovAngleDeg', ...
        'ImageWidth','ImageHeight','PixelLengthX','PixelLengthY','NClassified', ...
        'NWithAnatomy'})];
end
end

function T = abundancePairedTests(S,C)
morph=strings(0,1); nFish=zeros(0,1); meanCW=nan(0,1); cwCILow=nan(0,1); cwCIHigh=nan(0,1); meanCCW=nan(0,1); ccwCILow=nan(0,1); ccwCIHigh=nan(0,1); tStatistic=nan(0,1); degreesFreedom=nan(0,1); pValue=nan(0,1); cohensDz=nan(0,1);
for m=1:numel(C.MorphOrder)
 q=S.Morph==C.MorphOrder(m)&isfinite(S.PropCW)&isfinite(S.PropCCW); cw=S.PropCW(q); ccw=S.PropCCW(q); d=cw-ccw; n=numel(d);
 ciCW=bootstrapMeanCI(cw,2000); ciCCW=bootstrapMeanCI(ccw,2000); t=NaN; df=NaN; p=NaN; dz=NaN;
 if n>=2&&std(d)>0, [~,p,~,st]=ttest(cw,ccw); t=st.tstat; df=st.df; dz=mean(d)/std(d); elseif n>=2&&all(d==0), t=0; df=n-1; p=1; dz=0; end
 morph(end+1,1)=C.MorphOrder(m); nFish(end+1,1)=n; meanCW(end+1,1)=mean(cw,"omitnan"); cwCILow(end+1,1)=ciCW(1); cwCIHigh(end+1,1)=ciCW(2); meanCCW(end+1,1)=mean(ccw,"omitnan"); ccwCILow(end+1,1)=ciCCW(1); ccwCIHigh(end+1,1)=ciCCW(2); tStatistic(end+1,1)=t; degreesFreedom(end+1,1)=df; pValue(end+1,1)=p; cohensDz(end+1,1)=dz; %#ok<AGROW>
end
T=table(morph,nFish,meanCW,cwCILow,cwCIHigh,meanCCW,ccwCILow,ccwCIHigh,tStatistic,degreesFreedom,pValue,cohensDz);
end

function [Fish,Stats] = computeAnatomySkew(Cel,C)

keys = unique(Cel(:,{'Morph','Fish'}),"rows","stable");
Fish = table();

for i = 1:height(keys)

    q = Cel.Morph == keys.Morph(i) & ...
        Cel.Fish == keys.Fish(i) & ...
        isfinite(Cel.X) & Cel.X ~= 0;

    n = nnz(q);
    cwLeft = NaN;
    ccwRight = NaN;

    if n > 0
        right = Cel.X(q) > 0;
        left  = Cel.X(q) < 0;

        labels = Cel.Class(q);
        cw  = labels == "CW";
        ccw = labels == "CCW";

        cwLeft = mean(cw & left) - mean(cw)*mean(left);
        ccwRight = mean(ccw & right) - mean(ccw)*mean(right);
    end

    Fish = [Fish; table(...
        keys.Morph(i),keys.Fish(i),n,cwLeft,ccwRight,...
        'VariableNames',{'Morph','Fish','NValidAnatomy',...
        'CWLeftSkew','CCWRightSkew'})]; %#ok<AGROW>
end

Stats = table();
metrics = ["CWLeftSkew","CCWRightSkew"];

for m = 1:numel(C.MorphOrder)

    q = Fish.Morph == C.MorphOrder(m);

    % Test each covariance against zero
    for k = 1:numel(metrics)
        x = Fish.(metrics(k))(q);
        x = x(isfinite(x));

        [t,df,p,dz] = oneSampleTStats(x);

        Stats = [Stats; table(...
            C.MorphOrder(m),"one-sample vs zero",metrics(k),...
            numel(x),mean(x,"omitnan"),t,df,p,dz,...
            'VariableNames',{'Morph','Test','Metric','NFish',...
            'Estimate','TStatistic','DegreesFreedom','PValue',...
            'CohensDz'})]; %#ok<AGROW>
    end

    % Paired comparison within fish
    a = Fish.CWLeftSkew(q);
    b = Fish.CCWRightSkew(q);
    ok = isfinite(a) & isfinite(b);
    d = a(ok) - b(ok);

    [t,df,p,dz] = oneSampleTStats(d);

    Stats = [Stats; table(...
        C.MorphOrder(m),...
        "paired CW-Left vs CCW-Right",...
        "CWLeftSkew-CCWRightSkew",...
        nnz(ok),mean(d,"omitnan"),t,df,p,dz,...
        'VariableNames',{'Morph','Test','Metric','NFish',...
        'Estimate','TStatistic','DegreesFreedom','PValue',...
        'CohensDz'})]; %#ok<AGROW>
end
end

function [t,df,p,dz]=oneSampleTStats(x)
x=x(isfinite(x)); t=NaN; df=NaN; p=NaN; dz=NaN; n=numel(x);
if n>=2&&std(x)>0, [~,p,~,st]=ttest(x,0); t=st.tstat; df=st.df; dz=mean(x)/std(x); elseif n>=2&&all(x==0), t=0; df=n-1; p=1; dz=0; end
end

function T = rayleighTable(F,C)
T=table();
for i=1:numel(F)
    for c=1:3
        a=F(i).pref(F(i).labels==C.ClassOrder(c)); a=a(isfinite(a));
        [p,z]=rayleigh(a);
        T=[T;table(F(i).morph,F(i).fish,C.ClassOrder(c),numel(a),z,p, ...
            'VariableNames',{'Morph','Fish','Class','N','RayleighZ','P'})];
    end
end
end

function [p,z] = rayleigh(a)
n=numel(a); p=NaN; z=NaN;
if n<2, return; end
R=abs(sum(exp(1i*a))); z=R^2/n;
p=exp(-z)*(1+(2*z-z^2)/(4*n) ...
    -(24*z-132*z^2+76*z^3-9*z^4)/(288*n^2));
p=min(max(p,0),1);
end

function T = combineRayleigh(R,C)

T = table();

for m = C.MorphOrder
    for cls = ["CW","CCW"]

        q = R.Morph == m & R.Class == cls & isfinite(R.P);
        pFish = max(R.P(q),realmin);
        nFish = numel(pFish);

        fisherStatistic = NaN;
        combinedP = NaN;

        if nFish > 0
            fisherStatistic = -2*sum(log(pFish));
            combinedP = chi2cdf(...
                fisherStatistic,2*nFish,"upper");
        end

        T = [T; table(...
            m,cls,nFish,fisherStatistic,combinedP,...
            'VariableNames',{'Morph','Class','NFish',...
            'FisherStatistic','CombinedP'})]; %#ok<AGROW>
    end
end

T.HolmP = holmAdjust(T.CombinedP);
T.Significant = T.HolmP < C.Alpha;
end

function adjustedP = holmAdjust(p)

adjustedP = nan(size(p));
valid = find(isfinite(p));

[sortedP,order] = sort(p(valid));
n = numel(sortedP);

scaledP = (n-(1:n)'+1).*sortedP;
scaledP = cummax(scaledP);
scaledP = min(scaledP,1);

adjustedP(valid(order)) = scaledP;
end

function [Fish,Summary,Comparisons] = summarizeRingRayleigh(R,C)
Fish=R(ismember(R.Class,["CW","CCW"]),:); Fish.BonferroniAlpha=nan(height(Fish),1); Fish.Significant=false(height(Fish),1);
for m=1:numel(C.MorphOrder)
 nMorph=numel(unique(R.Fish(R.Morph==C.MorphOrder(m)))); threshold=C.Alpha/max(nMorph,1); q=Fish.Morph==C.MorphOrder(m); Fish.BonferroniAlpha(q)=threshold; Fish.Significant(q)=isfinite(Fish.P(q))&Fish.P(q)<threshold;
end
Summary=table();
for m=1:numel(C.MorphOrder)
 for cls=["CW","CCW"]
  q=Fish.Morph==C.MorphOrder(m)&Fish.Class==cls&isfinite(Fish.P); n=nnz(q); nSig=nnz(Fish.Significant(q)); frac=nSig/max(n,1); threshold=unique(Fish.BonferroniAlpha(q)); if isempty(threshold), threshold=NaN; else, threshold=threshold(1); end
  Summary=[Summary;table(C.MorphOrder(m),cls,n,nSig,frac,threshold,'VariableNames',{'Morph','Class','NFishTested','NSignificant','FractionSignificant','BonferroniAlpha'})]; %#ok<AGROW>
 end
end
Comparisons=table();
for cls=["CW","CCW"]
 surface=Summary(Summary.Morph=="Surface"&Summary.Class==cls,:);
 for cave=["Molino","Pachon"]
  other=Summary(Summary.Morph==cave&Summary.Class==cls,:); p=NaN; oddsRatio=NaN;
  if ~isempty(surface)&&~isempty(other)
   counts=[surface.NSignificant surface.NFishTested-surface.NSignificant; other.NSignificant other.NFishTested-other.NSignificant];
   if all(sum(counts,2)>0), [~,p,st]=fishertest(counts); if isfield(st,"OddsRatio"), oddsRatio=st.OddsRatio; end; end
  end
  Comparisons=[Comparisons;table(cls,"Surface",cave,p,oddsRatio,'VariableNames',{'Class','GroupA','GroupB','FisherP','OddsRatioBvsA'})]; %#ok<AGROW>
 end
end
end

function Stats=buildClassFeatureCrossMorphKW(S,AnatomyFish,RayleighFish,C)
Stats=table();
abundanceMetrics=["PropCW","PropCCW","PropSymmetric","CellLogCCWOverCW"];
for k=1:numel(abundanceMetrics)
 Stats=[Stats;pairwiseMorphKW(S.Morph,S.(abundanceMetrics(k)), ...
  "Abundance",abundanceMetrics(k),abundanceMetrics(k),C)]; %#ok<AGROW>
end
anatomyMetrics=["CWLeftSkew","CCWRightSkew"];
for k=1:numel(anatomyMetrics)
 Stats=[Stats;pairwiseMorphKW(AnatomyFish.Morph,AnatomyFish.(anatomyMetrics(k)), ...
  "AnatomySkew",anatomyMetrics(k),anatomyMetrics(k),C)]; %#ok<AGROW>
end
for cls=["CW","CCW"]
 q=RayleighFish.Class==cls; values=-log10(max(RayleighFish.P(q),realmin));
 Stats=[Stats;pairwiseMorphKW(RayleighFish.Morph(q),values, ...
  "RingRayleigh",cls,"negativeLog10RayleighP",C)]; %#ok<AGROW>
end
end

function T=pairwiseMorphKW(morph,values,figureName,panel,metric,C)
nRows=size(C.PairwiseMorphComparisons,1);
Figure=repmat(string(figureName),nRows,1); Panel=repmat(string(panel),nRows,1);
Metric=repmat(string(metric),nRows,1); Test=repmat("two-group Kruskal-Wallis",nRows,1);
groupA=strings(nRows,1); groupB=strings(nRows,1); nA=zeros(nRows,1); nB=zeros(nRows,1);
pRaw=nan(nRows,1); pAdjusted=nan(nRows,1); correctionMethod=repmat("none (planned comparisons)",nRows,1);
for j=1:nRows
 groupA(j)=string(C.PairwiseMorphComparisons{j,1}); groupB(j)=string(C.PairwiseMorphComparisons{j,2});
 a=values(string(morph)==groupA(j)); b=values(string(morph)==groupB(j)); a=a(isfinite(a)); b=b(isfinite(b));
 nA(j)=numel(a); nB(j)=numel(b);
 if nA(j)>=C.MinFishPerMorphKW&&nB(j)>=C.MinFishPerMorphKW
  groups=[repmat(groupA(j),nA(j),1);repmat(groupB(j),nB(j),1)];
  pRaw(j)=kruskalwallis([a;b],groups,"off");
 end
 pAdjusted(j)=pRaw(j);
end
T=table(Figure,Panel,Metric,Test,groupA,groupB,nA,nB,pRaw,pAdjusted,correctionMethod);
end

function path = plotAnatomySkew(Fish,Stats,KW,C)

fig = figure("Color","w","Visible",C.Visible,...
    "Position",[60 80 800 460]);
tl = tiledlayout(fig,1,2,...
    "TileSpacing","compact","Padding","compact");

metrics = ["CWLeftSkew","CCWRightSkew"];
labels = [...
    "p(CW \cap Left) - p(CW)p(Left)",...
    "p(CCW \cap Right) - p(CCW)p(Right)"];

for k = 1:2
    ax = nexttile(tl);
    fishBoxplots(ax,Fish,Fish.(metrics(k)),C);
    pbaspect(ax,[1.15 1 1]);
    yline(ax,0,"k--","LineWidth",1.1);
    ylabel(ax,labels(k));
    addClassFeatureKWBrackets(ax,KW(KW.Figure=="AnatomySkew"&KW.Panel==metrics(k),:),C.MorphOrder);

    lines = strings(3,1);
    for m = 1:3
        s = Stats(...
            Stats.Morph == C.MorphOrder(m) & ...
            Stats.Test == "one-sample vs zero" & ...
            Stats.Metric == metrics(k),:);

        if ~isempty(s)
            lines(m) = sprintf(...
                "%s: t(%g)=%.2f, p=%.3g, dz=%.2f",...
                C.MorphOrder(m),s.DegreesFreedom,...
                s.TStatistic,s.PValue,s.CohensDz);
        end
    end

    text(ax,.02,.02,strjoin(cellstr(lines),newline),...
        "Units","normalized","VerticalAlignment","bottom",...
        "BackgroundColor","w","Margin",3,"FontSize",8);
end

paired = Stats(Stats.Test == ...
    "paired CW-Left vs CCW-Right",:);

lines = strings(height(paired),1);
for i = 1:height(paired)
    lines(i) = sprintf(...
        "%s paired: t(%g)=%.2f, p=%.3g, dz=%.2f",...
        paired.Morph(i),paired.DegreesFreedom(i),...
        paired.TStatistic(i),paired.PValue(i),...
        paired.CohensDz(i));
end

title(tl,{...
    "Fish-level directional-class / side covariance",...
    "Within-morph paired tests: " + strjoin(lines,"; ")});

path = saveBoth(fig,C.OutputDir,...
    "08_CW_left_CCW_right_skewness");
end

function path = plotRingRayleigh(Fish,Combined,KW,C)

fig = figure("Color","w","Visible",C.Visible,...
    "Position",[70 70 1100 500]);

tl = tiledlayout(fig,1,2,...
    "TileSpacing","compact","Padding","compact");

classes = ["CW","CCW"];

for c = 1:2

    ax = nexttile(tl);
    hold(ax,"on");

    globalLines = strings(numel(C.MorphOrder),1);

    for m = 1:numel(C.MorphOrder)

        morph = C.MorphOrder(m);

        % Per-fish Rayleigh p-values
        q = Fish.Class == classes(c) & ...
            Fish.Morph == morph & ...
            isfinite(Fish.P);

        pFish = max(Fish.P(q),realmin);
        y = -log10(pFish);

        if ~isempty(y)

            boxchart(ax,repmat(m,numel(y),1),y,...
                "BoxFaceColor",C.MorphColors(m,:),...
                "BoxFaceAlpha",0.25,...
                "WhiskerLineColor",C.MorphColors(m,:),...
                "MarkerStyle","none");

            jitter = linspace(-0.10,0.10,numel(y))';
            if numel(y) == 1
                jitter = 0;
            end

            scatter(ax,m+jitter,y,32,...
                C.MorphColors(m,:),...
                "filled",...
                "MarkerEdgeColor","k",...
                "LineWidth",0.5);
        end

        % Morph-level Fisher combination
        g = Combined(...
            Combined.Morph == morph & ...
            Combined.Class == classes(c),:);

        if ~isempty(g)

            if g.Significant
                status = "significant";
            else
                status = "not significant";
            end

            globalLines(m) = sprintf(...
                "%s: %s; Fisher p=%s, Holm p=%s, n=%d",...
                string(morph),status,...
                char(pValueText(g.CombinedP)),...
                char(pValueText(g.HolmP)),...
                g.NFish);
        end
    end

    yline(ax,-log10(0.05),"k--","p = 0.05",...
        "LabelHorizontalAlignment","left");

    xlim(ax,[0.5 numel(C.MorphOrder)+0.5]);
    xticks(ax,1:numel(C.MorphOrder));
    xticklabels(ax,C.MorphOrder);

    ylabel(ax,"-\log_{10}(per-fish Rayleigh p)");
    grid(ax,"on");
    addClassFeatureKWBrackets(ax,KW(KW.Figure=="RingRayleigh"&KW.Panel==classes(c),:),C.MorphOrder);

    title(ax,{...
        classes(c) + " preferred-HD non-uniformity",...
        strjoin(globalLines,newline)},...
        "FontSize",10);
end

title(tl,...
    "Morph-level non-uniformity from Fisher-combined fish-level Rayleigh tests");

path = saveBoth(fig,C.OutputDir,...
    "07_ring_rayleigh_morph_global");

end

function path = plotAbundance(S,Stats,KW,C)
previousRng=rng; rngCleanup=onCleanup(@()rng(previousRng)); %#ok<NASGU>
rng(C.RandomSeed,"twister");
fig=figure("Color","w","Visible",C.Visible,"Position",[50 80 1450 560]);
t=tiledlayout(fig,1,4,"TileSpacing","compact","Padding","compact");
vars={"PropCW","PropCCW","PropSymmetric"};
labs={"CW proportion","CCW proportion","Symmetric proportion"};
for j=1:3
    ax=nexttile(t);
    plotAbundanceProportion(ax,S,S.(vars{j}),C);
    yline(ax,1/3,"k--","LineWidth",1.2,"Label","1/3 reference", ...
        "LabelHorizontalAlignment","left");
    ylabel(ax,labs{j}); ylim(ax,[0 1]);
    addClassFeatureKWBrackets(ax,KW(KW.Figure=="Abundance"&KW.Panel==string(vars{j}),:),C.MorphOrder);
end
ax=nexttile(t); fishBoxplots(ax,S,S.CellLogCCWOverCW,C);
yline(ax,0,"k:"); ylabel(ax,"log(CCW/CW cells)");
addClassFeatureKWBrackets(ax,KW(KW.Figure=="Abundance"&KW.Panel=="CellLogCCWOverCW",:),C.MorphOrder);
statsLines=strings(height(Stats),1);
for i=1:height(Stats)
 statsLines(i)=sprintf("%s: CW ring (%.1f%%, 95%% CI [%.1f%%, %.1f%%]) and CCW ring (%.1f%%, 95%% CI [%.1f%%, %.1f%%]); paired t test: t(%g)=%.2f, p=%.3f, dz=%.2f",Stats.morph(i),100*Stats.meanCW(i),100*Stats.cwCILow(i),100*Stats.cwCIHigh(i),100*Stats.meanCCW(i),100*Stats.ccwCILow(i),100*Stats.ccwCIHigh(i),Stats.degreesFreedom(i),Stats.tStatistic(i),Stats.pValue(i),Stats.cohensDz(i));
end
title(t,[{"Class abundance: fish-level proportions and within-morph paired CW vs CCW tests"};cellstr(statsLines)],"FontSize",10);
path=saveBoth(fig,C.OutputDir,"01_class_abundance");
end

function plotAbundanceProportion(ax,T,v,C)
hold(ax,"on");
for m=1:numel(C.MorphOrder)
    values=v(T.Morph==C.MorphOrder(m)); values=values(isfinite(values));
    if isempty(values), continue; end
    jitter=linspace(-.10,.10,numel(values))';
    if numel(values)==1, jitter=0; end
    scatter(ax,m+jitter,values,32,C.MorphColors(m,:),"filled", ...
        "MarkerEdgeColor","k","MarkerFaceAlpha",.65);
    estimate=mean(values,"omitnan");
    ci=bootstrapMeanCI(values,2000);
    errorbar(ax,m,estimate,estimate-ci(1),ci(2)-estimate,"k", ...
        "LineStyle","none","LineWidth",1.7,"CapSize",10);
    plot(ax,m,estimate,"kd","MarkerFaceColor",C.MorphColors(m,:), ...
        "MarkerSize",8,"LineWidth",1.1);
end
set(ax,"XTick",1:numel(C.MorphOrder),"XTickLabel",C.MorphOrder, ...
    "TickDir","out","Box","off");
xlim(ax,[.5 numel(C.MorphOrder)+.5]); grid(ax,"on");
end

function ci=bootstrapMeanCI(values,nBootstrap)
values=values(isfinite(values));
if isempty(values), ci=[NaN NaN]; return; end
if numel(values)==1, ci=[values values]; return; end
indices=randi(numel(values),numel(values),nBootstrap);
bootstrapMeans=mean(values(indices),1,"omitnan");
ci=prctile(bootstrapMeans,[2.5 97.5]);
end

function path = plotAnatomyDensityDifference(D,C)
fig=figure('Color','w','Visible',C.Visible,'Position',[60 80 1250 390]);
t=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
mx=max(abs(D.differenceMean),[],'all','omitnan');
if ~isfinite(mx)||mx<=0, mx=1; end
for m=1:3
    ax=nexttile(t);
    imagesc(ax,D.xCenters,D.yCenters,D.differenceMean(:,:,m));
    set(ax,'YDir','normal'); axis(ax,'equal'); xlim(ax,D.xLim); ylim(ax,D.yLim);
    caxis(ax,[-mx mx]); xline(ax,0,'k:'); yline(ax,0,'k:');
    xlabel(ax,C.AnatomyXLabel); ylabel(ax,C.AnatomyYLabel);
    title(ax,sprintf('%s | n=%d fish',C.MorphOrder(m),D.differenceNFish(m)));
end
colormap(fig,divergingMap(256)); colorbar(nexttile(t,3));
title(t,'Within-fish normalized CW density minus CCW density');
path=saveBoth(fig,C.OutputDir,'02d_CW_minus_CCW_anatomical_density');
end

function map = divergingMap(n)
if nargin<1, n=256; end
x=linspace(0,1,n)';
blue=[.15 .35 .85]; white=[1 1 1]; red=[.85 .20 .15];
map=zeros(n,3); lo=x<=.5; hi=~lo;
u=x(lo)/.5; map(lo,:)=blue.*(1-u)+white.*u;
u=(x(hi)-.5)/.5; map(hi,:)=white.*(1-u)+red.*u;
end

function paths = plotAnatomyQC(F,C)
out=fullfile(C.OutputDir,'anatomy_qc'); if ~isfolder(out), mkdir(out); end
paths=cell(numel(F),1);
for i=1:numel(F)
    A=F(i).anatomy;
    fig=figure('Color','w','Visible',C.Visible,'Position',[50 60 1300 560]);
    t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(t); imagesc(ax,A.background); colormap(ax,gray); axis(ax,'image');
    set(ax,'YDir','reverse'); hold(ax,'on');
    scatter(ax,A.xyRawAll(:,1),A.xyRawAll(:,2),5,[.8 .8 .8],'filled', ...
        'HandleVisibility','off');
    if ~isempty(A.linePos)
        plot(ax,A.linePos(:,1),A.linePos(:,2),'-','Color',[.1 .9 .2], ...
            'LineWidth',2,'DisplayName','Caudal to rostral');
        scatter(ax,A.linePos(1,1),A.linePos(1,2),45,[.1 .9 .2],'o','filled', ...
            'HandleVisibility','off');
        scatter(ax,A.linePos(2,1),A.linePos(2,2),55,[.1 .9 .2],'^','filled', ...
            'HandleVisibility','off');
    end
    scatter(ax,A.imageCenter(1),A.imageCenter(2),70,'m','x','LineWidth',2, ...
        'HandleVisibility','off');
    for c=1:3
        q=F(i).labels==C.ClassOrder(c); id=F(i).originalIds(q);
        scatter(ax,A.xyRawAll(id,1),A.xyRawAll(id,2),26,C.ClassColors(c,:), ...
            'filled','MarkerEdgeColor','k','DisplayName',C.ClassOrder(c));
    end
    title(ax,'Raw image, ROI centroids and caudal-to-rostral line');
    legend(ax,'Location','best');

    ax=nexttile(t); hold(ax,'on');
    scatter(ax,A.xAll,A.yAll,7,[.82 .82 .82],'filled','HandleVisibility','off');
    for c=1:3
        q=F(i).labels==C.ClassOrder(c);
        scatter(ax,F(i).x(q),F(i).y(q),30,C.ClassColors(c,:), ...
            'filled','MarkerEdgeColor','k','DisplayName',C.ClassOrder(c));
    end
    xline(ax,0,'k:','HandleVisibility','off');
    yline(ax,0,'k:','HandleVisibility','off');
    axis(ax,'equal'); grid(ax,'on'); box(ax,'off');
    xlabel(ax,C.AnatomyXLabel); ylabel(ax,C.AnatomyYLabel);
    title(ax,sprintf(['Aligned coordinates; FOV angle %.3f deg\n' ...
        '%s; median mapping r = %.4f'],A.fovAngleDeg,F(i).mappingMethod, ...
        median(F(i).mappingCorrelation,'omitnan')),'Interpreter','none');
    legend(ax,'Location','best');

    title(t,sprintf('Anatomy alignment QC | %s | %s',F(i).morph,F(i).fish));
    paths{i}=saveBoth(fig,out,sprintf('QC_%s_%s',F(i).morph,F(i).fish));
end
end

function path = plotAHVProfiles(P,C,centers)
fig=figure('Color','w','Visible',C.Visible,'Position',[100 40 540 950]);
t=tiledlayout(fig,3,1,'TileSpacing','compact');
for m=1:3
    ax=nexttile(t,m); hold(ax,'on');
    for c=1:3
        q=P.Morph==C.MorphOrder(m)&P.Class==C.ClassOrder(c);
        fish=unique(P.Fish(q)); M=nan(numel(fish),numel(centers));
        for i=1:numel(fish)
            u=P(q&P.Fish==fish(i),:); [~,o]=sort(u.AHVBin); M(i,:)=u.DeltaFoF(o);
        end
        mu=mean(M,1,'omitnan'); se=std(M,0,1,'omitnan')/sqrt(size(M,1));
        fill(ax,[centers fliplr(centers)],[mu-se fliplr(mu+se)], ...
            C.ClassColors(c,:),'FaceAlpha',.16,'EdgeColor','none','HandleVisibility','off');
        plot(ax,centers,mu,'Color',C.ClassColors(c,:),'LineWidth',2, ...
            'DisplayName',C.ClassOrder(c));
    end
    xline(ax,0,'k:','HandleVisibility','off');
    title(ax,C.MorphOrder(m)); xlabel(ax,'Angular head velocity (deg/s; CCW +)');
    ylabel(ax,'\DeltaF/F_0'); grid(ax,'on');
    pbaspect(ax,[1.35 1 1]);
    if m==1, legend(ax,'Location','best'); end
end
title(t,'AHV response profiles in raw \DeltaF/F_0');
path=saveBoth(fig,C.OutputDir,'04_ahv_response_profiles');
end

function path = plotSpeed(S,C)
fig=figure('Color','w','Visible',C.Visible,'Position',[50 80 1200 420]);
t=tiledlayout(fig,1,3);
speedColors=[.55 .55 .55; .25 .55 .85; .85 .35 .25];
for c=1:3
    ax=nexttile(t); M=nan(3,3);
    for m=1:3
        for s=1:3
            q=S.DirectionClass==C.ClassOrder(c)&S.Morph==C.MorphOrder(m)& ...
                S.SpeedClass==C.SpeedOrder(s);
            M(m,s)=mean(S.ProportionAllCells(q),'omitnan');
        end
        rowTotal=sum(M(m,:),'omitnan');
        if rowTotal>0, M(m,:)=M(m,:)/rowTotal; end
    end
    b=bar(ax,M,'stacked');
    for s=1:3, b(s).FaceColor=speedColors(s,:); end
    ylim(ax,[0 1]); set(ax,'XTickLabel',C.MorphOrder);
    ylabel(ax,'Proportion within direction class');
    title(ax,C.ClassOrder(c)); grid(ax,'on');
    if c==3, legend(ax,C.SpeedOrder,'Location','eastoutside'); end
end
title(t,'Direction class x speed-modulation class');
path=saveBoth(fig,C.OutputDir,'05_speed_modulation');
end

function path = plotPreferred(T,C)
previousRng=rng; rngCleanup=onCleanup(@()rng(previousRng));
rng(C.RandomSeed,"twister");
fig=figure("Color","w","Visible",C.Visible,"Position",[50 80 1350 560]);
tl=tiledlayout(fig,3,3,"TileSpacing","compact","Padding","compact");
edges=linspace(-pi,pi,19);
centers=rad2deg((edges(1:end-1)+edges(2:end))/2);
for m=1:3
    for c=1:3
        ax=nexttile(tl,(m-1)*3+c); hold(ax,"on");
        q=T.Morph==C.MorphOrder(m)&T.Class==C.ClassOrder(c)&isfinite(T.PreferredRad);
        fish=unique(T.Fish(q),"stable");
        fishPercent=nan(numel(fish),numel(centers));
        for f=1:numel(fish)
            angles=T.PreferredRad(q&T.Fish==fish(f));
            counts=histcounts(angles,edges);
            if sum(counts)>0, fishPercent(f,:)=100*counts/sum(counts); end
        end
        estimate=mean(fishPercent,1,"omitnan");
        ci=nan(numel(centers),2);
        for b=1:numel(centers)
            ci(b,:)=bootstrapMeanCI(fishPercent(:,b),2000);
        end
        bar(ax,centers,estimate,1,"FaceColor",C.ClassColors(c,:), ...
            "EdgeColor","none","FaceAlpha",.75);
        lowerError=estimate-transpose(ci(:,1));
        upperError=transpose(ci(:,2))-estimate;
        errorbar(ax,centers,estimate,lowerError,upperError,"k", ...
            "LineStyle","none","LineWidth",.8,"CapSize",2);
        xlim(ax,[-180 180]); xticks(ax,[-180 0 180]);
        title(ax,sprintf("%s | %s | n=%d fish",C.MorphOrder(m),C.ClassOrder(c),numel(fish)));
        xlabel(ax,"Preferred direction (deg)"); ylabel(ax,"Cells (%)");
        grid(ax,"on"); pbaspect(ax,[2.2 1 1]);
    end
end
title(tl,{"Preferred head direction by class and morph", ...
    "Fish-weighted bin percentages with fish-bootstrap 95% CIs"});
path=saveBoth(fig,C.OutputDir,"07_preferred_head_direction");
end

function addClassFeatureKWBrackets(ax,sub,morphOrder)
% Same bracket geometry and significance mapping as the v3 model comparison.
if isempty(sub)||~ismember("pAdjusted",string(sub.Properties.VariableNames)), return; end
sub=sub(isfinite(sub.pAdjusted),:); if isempty(sub), return; end
limits=ylim(ax); dataRange=limits(2)-limits(1);
if ~isfinite(dataRange)||dataRange<=0, dataRange=max(1,abs(limits(2))); end
base=limits(2)+.08*dataRange; step=.13*dataRange; tick=.035*dataRange;
for j=1:height(sub)
 x1=find(strcmpi(string(morphOrder),sub.groupA(j)),1); x2=find(strcmpi(string(morphOrder),sub.groupB(j)),1);
 if isempty(x1)||isempty(x2), continue; end
 y=base+(j-1)*step;
 plot(ax,[x1 x1 x2 x2],[y-tick y y y-tick],"k-","LineWidth",1.15, ...
  "HandleVisibility","off","Clipping","off");
 text(ax,mean([x1 x2]),y+.02*dataRange,classFeatureSignificanceLabel(sub.pAdjusted(j)), ...
  "HorizontalAlignment","center","VerticalAlignment","bottom","FontWeight","bold");
end
ylim(ax,[limits(1) base+max(height(sub)-1,0)*step+.16*dataRange]);
pieces=strings(height(sub),1);
for j=1:height(sub)
 pieces(j)=sub.groupA(j)+"-"+sub.groupB(j)+": p="+string(sprintf("%.3g",sub.pAdjusted(j)));
end
text(ax,.02,.98,strjoin(cellstr(pieces),newline),"Units","normalized", ...
 "HorizontalAlignment","left","VerticalAlignment","top","FontSize",8, ...
 "Interpreter","none","BackgroundColor","w","Margin",2,"Clipping","off");
end

function label=classFeatureSignificanceLabel(p)
if p<.001, label="***"; elseif p<.01, label="**"; elseif p<.05, label="*"; else, label="ns"; end
end

function fishBoxplots(ax,T,v,C)
hold(ax,'on');
for m=1:3
    x=v(T.Morph==C.MorphOrder(m)); x=x(isfinite(x));
    if isempty(x), continue; end
    boxchart(ax,repmat(m,numel(x),1),x,'BoxFaceColor',C.MorphColors(m,:), ...
        'MarkerStyle','none','BoxWidth',.55);
    xx=m+linspace(-.10,.10,numel(x))';
    scatter(ax,xx,x,26,.4+(1-.4)*C.MorphColors(m,:),'filled', ...
        'MarkerEdgeColor','k','HandleVisibility','off');
end
set(ax,'XTick',1:3,'XTickLabel',C.MorphOrder,'XTickLabelRotation',15, ...
    'TickDir','out','Box','off');
grid(ax,'on');
end

function s = pValueText(p)
if isempty(p) || ~isfinite(p)
    s = "NA";
elseif p < 1e-4
    s = string(sprintf("%.1e",p));
else
    s = string(sprintf("%.4f",p));
end
end

function p = saveBoth(fig,out,name)
out=char(string(out)); name=char(string(name));
p=struct('png',fullfile(out,[name '.png']),'svg',fullfile(out,[name '.svg']));
exportgraphics(fig,p.png,'Resolution',200,'ContentType','image');
% Keep the SVG as real vector artwork so plot objects remain separately
% selectable/editable in applications such as Illustrator or Inkscape.
% ContentType='image' rasterizes the complete figure into one embedded image.
set(fig,'Renderer','painters');
exportgraphics(fig,p.svg,'ContentType','vector');
close(fig);
end
