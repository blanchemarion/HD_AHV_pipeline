function [Y,provenance] = preprocessActivity(Y,doZscore)
%PREPROCESSACTIVITY Apply the fitter's finite, per-neuron timewise z-score.
Y=double(Y);
if doZscore
    for j=1:size(Y,2)
        ok=isfinite(Y(:,j));
        if any(ok)
            mu=mean(Y(ok,j)); sigma=std(Y(ok,j));
            if sigma>0; Y(ok,j)=(Y(ok,j)-mu)./sigma;
            else; Y(ok,j)=0; end
        end
    end
end
provenance=struct('zscoreActivityWithinWindow',logical(doZscore), ...
    'scope','independently within supplied time window, per neuron', ...
    'finiteHandling','mean/std calculated from finite response samples only');
end
