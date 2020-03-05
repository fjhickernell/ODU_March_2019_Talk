function [KDataData, KPlotData, coeff, fAppPlot, RMSPE, xBad, fBad, ...
   maxRMSPE, whMiss, Miss, thetaout] = ...
   RKHS(kname, f, xData, fData, xPlot, fPlot, theta)
%RKHS sets up stuff for approximation

d = size(xData,2);

%% Infer parameters
if strcmp(kname{1},{'MaternTh'})
   Ktheta = @(logth) kernel(kname,xData,xData,exp(logth));
   objective = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
   if d == 1
      logthopt = fminbnd(@(logth) objective(Ktheta(logth),fData),-5,5);
   else
      logthopt = fminsearch(@(logth) objective(Ktheta(logth),fData),zeros(1,d));
   end
   thetaout = exp(logthopt);
elseif strcmp(kname{1},{'MaternMod'})
   Ktheta = @(theta) kernel(kname,xData,xData,[exp(theta(1:d)), theta(d+1:2*d)]);
   objective = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
   [bopt,~] = fminsearch(@(b) objective(Ktheta(b),fData), zeros(1,2*d));
   thetaout = [exp(bopt(1:d)), bopt(d+1:2*d)];
else
   thetaout = theta;
end

%% Evaluate at parameters
KDataData = kernel(kname,xData,xData,thetaout);
KPlotData = kernel(kname,xPlot,xData,thetaout);
coeff = KDataData\fData;
fAppPlot = KPlotData*coeff;
normf = sqrt(coeff'*fData);

%% Compute RMSPE
if strcmp(kname{2},'Afix')
   A = kname{3};
end
RMSPE = real(sqrt(kernDiag(kname, xPlot, thetaout) - ...
   sum(KPlotData.*(KDataData\KPlotData')',2))) .* (A*normf);
[maxRMSPE,whBad] = max(RMSPE);
xBad = xPlot(whBad,:);
fBad = f(xBad);

whMiss = find((fPlot > fAppPlot + RMSPE + 1e6*eps) | ...
   (fPlot < fAppPlot - RMSPE - 1e5*eps));
Miss = [xPlot(whMiss,:) fPlot(whMiss) fAppPlot(whMiss) + [-1 1].*RMSPE(whMiss)];


end

function out = kernel(kname, x, y, theta)
   d = size(x,2);
   dist = @(x,y,theta) sqrt(sum((reshape(theta(1:d).*x,[size(x,1),1,d]) - ...
      reshape(theta.*y,[1,size(y,1),d])).^2,3));
   kdist = @(z) (1 + z).*exp(-z);
   if any(strcmp(kname{1}, {'Matern', 'MaternTh'}))
      out = kdist(dist(x,y,theta(1:d)));
   elseif strcmp(kname{1}, 'MaternMod')
      out = exp((x*theta(d+1:2*d)'+theta(d+1:2*d)*y')).* kdist(dist(x,y,theta(1:d)));
   else
      disp('Invalid')
   end
end

function out = kernDiag(kname, x, theta)
   d = size(x,2);
   if any(strcmp(kname{1}, {'Matern', 'MaternTh'}))
      out = ones(size(x,1),1);
   elseif strcmp(kname{1}, 'MaternMod')
      out = exp(x*(2*theta(d+1:2*d)'));
   else
      disp('Invalid')
   end
end

