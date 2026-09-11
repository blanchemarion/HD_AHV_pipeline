function [labels,pAdjusted] = classifyB2(b2,pB2,familyN,alpha,positiveIsCCW,valid)
if nargin<6; valid=isfinite(b2)&isfinite(pB2); end
b2=b2(:); pB2=pB2(:); valid=logical(valid(:)); pAdjusted=nan(size(pB2)); pAdjusted(valid)=min(1,familyN.*pB2(valid)); labels=repmat("NotFit",numel(b2),1); labels(valid)="Symmetric"; sig=valid&pAdjusted<alpha&b2~=0;
if positiveIsCCW; labels(sig&b2>0)="CCW"; labels(sig&b2<0)="CW"; else; labels(sig&b2>0)="CW"; labels(sig&b2<0)="CCW"; end
end
