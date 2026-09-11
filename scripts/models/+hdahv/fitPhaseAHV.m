function [F,Audit] = fitPhaseAHV(Y,phase,ahv,sourcePrefRad,minSamples)
%FITPHASEAHV Shared fixed-preference phase+AHV OLS computational core.
N=size(Y,2); z=nan(N,1);
F=struct('betaTheta',z,'b1',z,'b2',z,'c0',z,'seBetaTheta',z,'seB1',z,'seB2',z,'seC0',z,'pBetaTheta',z,'pB1',z,'pB2',z,'pC0',z,'r2',z,'designRank',z,'conditionNumberRaw',z,'conditionNumberStandardized',z,'nFitSamples',zeros(N,1),'failureReason',repmat("NotAttempted",N,1));
F.coefficientOrder={'betaTheta','b1','b2','c0'};
F.modelFormula='y = betaTheta*cos(sourcePrefRad-networkPhase) + b1*abs(AHV) + b2*AHV + c0 + epsilon';
phase=phase(:); ahv=ahv(:); sourcePrefRad=sourcePrefRad(:);
assert(size(Y,1)==numel(phase)&&numel(phase)==numel(ahv)); assert(size(Y,2)==numel(sourcePrefRad));
if nargout>1
 Audit=struct("coefficientOrder",{F.coefficientOrder},"phase",phase, ...
  "absAHV",abs(ahv),"signedAHV",ahv,"intercept",ones(numel(phase),1), ...
  "sourcePrefRad",sourcePrefRad,"phaseRegressor",cos(sourcePrefRad(:).'-phase), ...
  "finiteFrameMask",isfinite(Y)&isfinite(phase)&isfinite(ahv)&isfinite(sourcePrefRad(:).'));
else
 Audit=struct();
end
for i=1:N
 if ~isfinite(sourcePrefRad(i)); F.failureReason(i)="InvalidSourcePref"; continue; end
 X=[cos(sourcePrefRad(i)-phase),abs(ahv),ahv,ones(numel(phase),1)]; ok=all(isfinite(X),2)&isfinite(Y(:,i)); n=sum(ok); F.nFitSamples(i)=n;
 if n<minSamples; F.failureReason(i)="TooFewSamples"; continue; end
 ahvValid=ahv(ok); tolerance=1e-12;
 if ~any(ahvValid>tolerance)||~any(ahvValid<-tolerance); F.failureReason(i)="InsufficientDirectionalSupport"; continue; end
 X=X(ok,:); y=Y(ok,i); F.designRank(i)=rank(X); F.conditionNumberRaw(i)=cond(X);
 Z=X(:,1:3); mu=mean(Z,1); sd=std(Z,0,1); if any(sd==0); F.conditionNumberStandardized(i)=Inf; else; F.conditionNumberStandardized(i)=cond([(Z-mu)./sd,ones(n,1)]); end
 if F.designRank(i)<size(X,2); F.failureReason(i)="RankDeficientDesign"; continue; end
 beta=X\y;
 resid=y-X*beta; dof=n-size(X,2); if dof<=0; F.failureReason(i)="NoResidualDegreesOfFreedom"; continue; end
 mse=sum(resid.^2)/dof; se=sqrt(max(0,diag(mse*((X'*X)\eye(size(X,2)))))); p=2*(1-tcdf(abs(beta./se),dof));
 F.betaTheta(i)=beta(1); F.b1(i)=beta(2); F.b2(i)=beta(3); F.c0(i)=beta(4); F.seBetaTheta(i)=se(1); F.seB1(i)=se(2); F.seB2(i)=se(3); F.seC0(i)=se(4); F.pBetaTheta(i)=p(1); F.pB1(i)=p(2); F.pB2(i)=p(3); F.pC0(i)=p(4);
 sst=sum((y-mean(y)).^2); if sst>0; F.r2(i)=1-sum(resid.^2)/sst; end; F.failureReason(i)="";
end
end
