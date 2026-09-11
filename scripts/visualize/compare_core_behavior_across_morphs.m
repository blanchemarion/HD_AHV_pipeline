function results = compare_core_behavior_across_morphs(cfg)
%COMPARE_CORE_BEHAVIOR_ACROSS_MORPHS Original pass2 behavior, one row per fish.
% Usage (from the pipeline root):
%   addpath('scripts','scripts/visualize');
%   results = compare_core_behavior_across_morphs(pipeline_config());
%
% Inclusion: nonempty successful FishSummary rows (nFitNeurons > 0) saved by
% fit_neuron_HD_AHV_three_models_behavior_tuned_cells_only_ROC. No additional
% visualization-specific neuron-count, model-R2, or rotation-count filter.
% Each morph/session is treated as one fish, as in the source session list.
% Analyze the entire ORIGINAL *_swimResults_pass2.mat, not a cropped variant.
% Eligible turns use abs(dtheta) > the saved P.ahvMinAbsBoutAngleRad
% (strict, as upstream). Frequency = eligible turns / full recording minutes.
% Quiet intervals span consecutive eligible turns: next eligible onset minus
% preceding eligible end. Invalid/overlapping pairs are discarded.
% No leading/trailing censored interval is included. Angles follow the source
% degree-detection and [-pi,pi) wrapping convention.
%
% Bars: morph medians, fish-bootstrap 95%% CIs, individual fish dots.
% Tests: planned two-group Kruskal-Wallis, Surface-Molino and Surface-Pachon,
% uncorrected, minimum 3 fish/group; same settings as the v3 model comparison.
% Four panels cover three behavior summaries (quiet median and P90 separately).

if nargin < 1 || isempty(cfg)
    addpath(fileparts(fileparts(mfilename('fullpath'))));
    cfg = pipeline_config();
end
P = struct('morphOrder',{{'Surface','Molino','Pachon'}}, ...
    'colors',[.85 .20 .20; .20 .65 .30; .95 .72 .10], ...
    'nBootstrap',10000,'randomSeed',7,'minFishPerMorph',3,'alpha',.05);
P.figure = struct('pointSize',46,'jitterWidth',.11,'morphBarWidth',.46, ...
    'morphXLimits',[.72 3.28],'dpi',300);
P.inclusion = 'Saved FishSummary with nFitNeurons > 0; one session per fish';
P.quietIntervals = 'Consecutive eligible turns after applying the saved threshold';
P.duration = 'Last minus first frame timestamp plus median frame duration';
outDir = fullfile(cfg.OutputsDir,'core_behavior_comparison');
if ~isfolder(outDir); mkdir(outDir); end
oldRng = rng; restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(P.randomSeed);
fishRows = {}; boutRows = {}; selectionRows = {};
for m = 1:3
    morph = string(P.morphOrder{m});
    source = fullfile(cfg.ModelDataDir,sprintf( ...
        'HD_AHV_three_models_phase_tuned_only_results_%s.mat',lower(morph)));
    L = load(source,'P','FishSummary','sessions');
    assert(all(isfield(L,{'P','FishSummary','sessions'})), ...
        'Missing saved fitting metadata in %s.',source);
    threshold = L.P.ahvMinAbsBoutAngleRad;
    assert(isscalar(threshold) && isfinite(threshold) && threshold>=0);
    F = L.FishSummary;
    names = string(F.session);
    assert(numel(unique(names))==numel(names),'Duplicate fish-summary sessions.');
    included = names(F.nFitNeurons>0);
    discovered = string({L.sessions.name})';
    assert(all(ismember(included,discovered)),'Included session absent from metadata.');
    selectionRows{m} = table(repmat(morph,numel(discovered),1),discovered, ...
        ismember(discovered,included),repmat(string(source),numel(discovered),1), ...
        'VariableNames',{'morph','session','included','modelResultFile'}); %#ok<AGROW>
    fprintf('%s: %d included sessions; saved turn threshold %.15g rad\n', ...
        morph,numel(included),threshold);
    for s = 1:numel(included)
        session = included(s);
        folder = fullfile(cfg.DataRoot,lower(morph),session);
        files = dir(fullfile(folder,'*_swimResults_pass2.mat'));
        assert(numel(files)==1, ...
            'Expected exactly one original *_swimResults_pass2.mat in %s; found %d.', ...
            folder,numel(files));
        path = fullfile(files(1).folder,files(1).name);
        S = load(path,'pass2FileResult');
        assert(isfield(S,'pass2FileResult'),'Missing pass2FileResult in %s.',path);
        [row,bouts] = summarizeBehavior(S.pass2FileResult,threshold);
        row.morph = morph; row.session = session;
        row.behaviorFile = string(path); row.modelResultFile = string(source);
        bouts.morph = repmat(morph,height(bouts),1);
        bouts.session = repmat(session,height(bouts),1);
        fishRows{end+1} = row; boutRows{end+1} = bouts; %#ok<AGROW>
    end
end
assert(~isempty(fishRows),'No included fitted sessions.');
Fish = vertcat(fishRows{:}); Bouts = vertcat(boutRows{:});
Selection = vertcat(selectionRows{:});
metrics = {'eligibleTurnFrequencyPerMin','medianAbsTurnRad','medianQuietSec','p90QuietSec'};
labels = {'Eligible-turn frequency (turns/min)','Median |dtheta| of eligible turns (rad)', ...
    'Median quiet interval (s)','90th percentile quiet interval (s)'};
estimateRows = {}; testRows = {};
fig = figure('Color','w','Visible',cfg.FigureVisible,'Position',[60 60 720 720], ...
    'Renderer','painters');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for k = 1:numel(metrics)
    ax = nexttile(tl); hold(ax,'on'); values = cell(3,1);
    for m = 1:3
        x = Fish.(metrics{k})(Fish.morph==P.morphOrder{m});
        x = x(isfinite(x)); values{m} = x;
        [center,ci] = medianCI(x,P.nBootstrap);
        estimateRows{end+1} = table(string(metrics{k}),string(P.morphOrder{m}), ...
            numel(x),center,ci(1),ci(2),'VariableNames', ...
            {'metric','morph','nFish','median','ciLow','ciHigh'}); %#ok<AGROW>
        bar(ax,m,center,P.figure.morphBarWidth,'FaceColor',P.colors(m,:), ...
            'FaceAlpha',.78,'EdgeColor','k','LineWidth',.8, ...
            'HandleVisibility','off');
        errorbar(ax,m,center,center-ci(1),ci(2)-center,'k', ...
            'LineStyle','none','LineWidth',1.7,'CapSize',8, ...
            'HandleVisibility','off');
        scatter(ax,m+deterministicJitter(numel(x),P.figure.jitterWidth),x, ...
            max(18,round(P.figure.pointSize*.55)), ...
            'MarkerFaceColor',mixWithWhite(P.colors(m,:),.38), ...
            'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82, ...
            'HandleVisibility','off');
    end
    testLabels = strings(2,1);
    for j = 1:2
        a = values{1}; b = values{j+1}; p = NaN; effect = NaN; ci = [NaN NaN];
        if ~isempty(a) && ~isempty(b)
            effect = median(b)-median(a);
        end
        if numel(a)>=P.minFishPerMorph && numel(b)>=P.minFishPerMorph
            if numel(unique([a;b]))==1
                p = 1;
            else
                p = kruskalwallis([a;b],[ones(size(a));2*ones(size(b))],'off');
            end
            ba = median(a(randi(numel(a),numel(a),P.nBootstrap)),1);
            bb = median(b(randi(numel(b),numel(b),P.nBootstrap)),1);
            ci = prctile(bb-ba,[2.5 97.5]);
        end
        testRows{end+1} = table(string(metrics{k}),"Surface", ...
            string(P.morphOrder{j+1}),numel(a),numel(b),p,effect,ci(1),ci(2), ...
            'VariableNames',{'metric','groupA','groupB','nA','nB','pRaw', ...
            'medianDifferenceBminusA','differenceCILow','differenceCIHigh'}); %#ok<AGROW>
        testLabels(j) = sprintf('Surf-%s p=%s', ...
            shortMorph(P.morphOrder{j+1}),formatP(p));
    end
    tickLabels = arrayfun(@(m) sprintf('%s (n=%d fish)', ...
        P.morphOrder{m},numel(values{m})),1:3,'UniformOutput',false);
    set(ax,'XTick',1:3,'XTickLabel',tickLabels,'TickDir','out','Box','off');
    xlim(ax,P.figure.morphXLimits); ylabel(ax,labels{k},'Interpreter','tex');
    grid(ax,'on');
    title(ax,['Two-group KW: ' strjoin(cellstr(testLabels),'; ')], ...
        'FontWeight','normal');
    metricTests = vertcat(testRows{end-1:end});
    addPairwiseBrackets(ax,metricTests,P.morphOrder);
end
title(tl,'Original behavior: fish medians and bootstrap 95% CIs', ...
    'Interpreter','none','FontWeight','bold');
Estimates = vertcat(estimateRows{:}); Tests = vertcat(testRows{:});
Tests.test = repmat("Two-group Kruskal-Wallis",height(Tests),1);
Tests.multipleComparisonCorrection = repmat("none",height(Tests),1);
writetable(Fish,fullfile(outDir,'fish_behavior_summary.csv'));
writetable(Bouts,fullfile(outDir,'bout_and_quiet_interval_audit.csv'));
writetable(Selection,fullfile(outDir,'session_selection.csv'));
writetable(Estimates,fullfile(outDir,'morph_estimates.csv'));
writetable(Tests,fullfile(outDir,'planned_pairwise_statistics.csv'));
exportgraphics(fig,fullfile(outDir,'core_behavior_bar_charts.png'), ...
    'Resolution',P.figure.dpi);
exportgraphics(fig,fullfile(outDir,'core_behavior_bar_charts.svg'), ...
    'ContentType','vector');
legacyFig = fullfile(outDir,'core_behavior_bar_charts.fig');
if isfile(legacyFig); delete(legacyFig); end
if strcmpi(cfg.FigureVisible,'off'); close(fig); end
results = struct('fish',Fish,'bouts',Bouts,'selection',Selection, ...
    'estimates',Estimates,'tests',Tests,'parameters',P,'outputDir',outDir);
save(fullfile(outDir,'core_behavior_results.mat'),'results');
fprintf('Saved behavior summaries for %d fish to %s\n',height(Fish),outDir);
end

function [row,B] = summarizeBehavior(beh,threshold)
fps = double(beh.fps);
assert(isscalar(fps) && isfinite(fps) && fps>0,'Invalid behavior fps.');
n = NaN;
for name = {'behaviorTimeSec','heading_est','tail_angle'}
    if isfield(beh,name{1}) && ~isempty(beh.(name{1}))
        n = numel(beh.(name{1})); break;
    end
end
if isnan(n) && isfield(beh,'nFrames'); n = double(beh.nFrames); end
assert(isscalar(n) && isfinite(n) && n>=2 && n==round(n),'Invalid behavior length.');
t = (0:n-1)'/fps; timeSource = "fps";
if isfield(beh,'behaviorTimeSec') && ~isempty(beh.behaviorTimeSec)
    t = double(beh.behaviorTimeSec(:)); timeSource = "behaviorTimeSec";
end
assert(numel(t)==n && all(isfinite(t)) && all(diff(t)>0),'Invalid timestamps.');
durationSec = t(end)-t(1)+median(diff(t));
angleSource = "dtheta";
if ~isfield(beh,'dtheta') || isempty(beh.dtheta); angleSource = "LI"; end
assert(isfield(beh,angleSource),'Neither dtheta nor LI is present.');
d = double(beh.(angleSource)(:)); finiteAngles = d(isfinite(d));
convertedDegrees = ~isempty(finiteAngles) && max(abs(finiteAngles))>2*pi+.5;
if convertedDegrees; d = deg2rad(d); end
d = mod(d+pi,2*pi)-pi;
sf = double(beh.swimFrames);
if size(sf,1)==2; starts=sf(1,:)'; ends=sf(2,:)';
elseif size(sf,2)==2; starts=sf(:,1); ends=sf(:,2);
else; error('swimFrames must be 2 x N or N x 2.'); end
assert(numel(starts)==numel(d),'Bout/angle counts differ; refusing silent truncation.');
originalRow = (1:numel(d))';
[starts,order] = sort(round(starts)); ends=round(ends(order)); d=d(order);
originalRow = originalRow(order);
valid = isfinite(starts) & starts>=1 & starts<=n & isfinite(d);
eligible = valid & abs(d)>threshold;
validEnd = isfinite(ends) & ends>=starts & ends<=n;
quiet = nan(size(d));
eligibleRows = find(eligible);
if numel(eligibleRows) >= 2
    previous = eligibleRows(1:end-1);
    following = eligibleRows(2:end);
    pairOk = validEnd(previous) & starts(following)>=ends(previous);
    previous = previous(pairOk);
    following = following(pairOk);
    quiet(previous) = t(starts(following))-t(ends(previous));
end
q = quiet(isfinite(quiet)); turns = abs(d(eligible));
row = table(durationSec,numel(d),sum(valid),sum(eligible),sum(~validEnd), ...
    numel(q),threshold,sum(eligible)/(durationSec/60),safeMedian(turns), ...
    safeMedian(q),safePercentile(q,90),timeSource,angleSource,convertedDegrees, ...
    'VariableNames',{'durationSec','nDetectedBouts','nValidBouts','nEligibleTurns', ...
    'nInvalidBoutEnds','nQuietIntervals','turnThresholdRad','eligibleTurnFrequencyPerMin', ...
    'medianAbsTurnRad','medianQuietSec','p90QuietSec','timeSource','angleSource', ...
    'anglesConvertedFromDegrees'});
B = table(originalRow,starts,ends,d,valid,eligible,validEnd,quiet, ...
    'VariableNames',{'originalBoutRow','onsetFrame','endFrame','dthetaRad', ...
    'validBout','eligibleTurn','validEnd','quietToNextEligibleTurnSec'});
end

function x = safeMedian(v)
if isempty(v); x=NaN; else; x=median(v); end
end

function x = safePercentile(v,p)
if isempty(v); x=NaN; else; x=prctile(v,p); end
end

function [center,ci] = medianCI(x,nBoot)
center = safeMedian(x); ci = [NaN NaN];
if isempty(x); return; end
boot = median(x(randi(numel(x),numel(x),nBoot)),1);
ci = prctile(boot,[2.5 97.5]);
end


function jitter = deterministicJitter(n,width)
if n <= 1; jitter = 0; else; jitter = linspace(-width,width,n)'; end
end

function color = mixWithWhite(color,fraction)
color = color.*(1-fraction) + fraction;
end

function addPairwiseBrackets(ax,T,morphOrder)
T = T(isfinite(T.pRaw),:);
if isempty(T); return; end
limits = ylim(ax);
dataRange = limits(2)-limits(1);
if ~isfinite(dataRange) || dataRange<=0; dataRange=max(1,abs(limits(2))); end
base = limits(2)+.08*dataRange;
step = .13*dataRange;
tick = .035*dataRange;
for j = 1:height(T)
    x1=find(strcmpi(morphOrder,T.groupA(j)),1);
    x2=find(strcmpi(morphOrder,T.groupB(j)),1);
    if isempty(x1) || isempty(x2); continue; end
    y=base+(j-1)*step;
    plot(ax,[x1 x1 x2 x2],[y-tick y y y-tick],'k-', ...
        'LineWidth',1.15,'HandleVisibility','off','Clipping','off');
    text(ax,mean([x1 x2]),y+.02*dataRange,significanceLabel(T.pRaw(j)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontWeight','bold');
end
ylim(ax,[limits(1) base+max(height(T)-1,0)*step+.16*dataRange]);
addPValueText(ax,T);
end

function addPValueText(ax,T)
pieces = strings(height(T),1);
for i=1:height(T)
    pieces(i)=string(T.groupA(i))+"-"+string(T.groupB(i))+ ...
        ": p="+string(sprintf('%.3g',T.pRaw(i)));
end
text(ax,.02,.98,strjoin(cellstr(pieces),newline),'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top','FontSize',8, ...
    'Interpreter','none','BackgroundColor','w','Margin',2,'Clipping','off');
end

function label = significanceLabel(p)
if p<.001; label='***'; elseif p<.01; label='**'; ...
elseif p<.05; label='*'; else; label='ns'; end
end

function s = shortMorph(morph)
switch char(morph)
    case 'Surface'; s='Surf';
    case 'Molino'; s='Mol';
    case 'Pachon'; s='Pach';
    otherwise; s=char(morph);
end
end

function s = formatP(p)
if ~isfinite(p); s='NA'; elseif p<1e-4; s=sprintf('%.1e',p); ...
else; s=sprintf('%.4f',p); end
end
