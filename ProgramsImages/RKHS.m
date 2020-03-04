function [KDataData, KPlotData, coeff, fAppPlot, RMSPE, xBad, fBad, ...
   maxRMSPE, whMiss, Miss, varargin] = ...
   RKHS(kname, f, xData, fData, xPlot, fPlot, varargin)
%RKHS sets up stuff for approximation

%% Infer parameters
if strcmp(kname{1},{'MaternTh'})
   Ktheta = @(logth) kernel(kname,xData,xData,exp(logth));
   objective = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
   logthopt = fminbnd(@(logth) objective(Ktheta(logth),fData),-5,5);
   varargin = {exp(logthopt)};
elseif strcmp(kname{1},{'MaternMod'})
   Ktheta = @(logth,a) kernel(kname,xData,xData,exp(logth),a);
   objective = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
   [bopt,objmin] = fminsearch(@(b) objective(Ktheta(b(1),b(2)),fData),[2,-10]);
   varargin = {exp(bopt(1)), bopt(2)};
end

%% Evaluate at parameters
KDataData = kernel(kname,xData,xData,varargin{:});
KPlotData = kernel(kname,xPlot,xData,varargin{:});
coeff = KDataData\fData;
fAppPlot = KPlotData*coeff;
normf = sqrt(coeff'*fData);

%% Compute RMSPE
if strcmp(kname{2},'Afix')
   A = kname{3};
end
RMSPE = real(sqrt(kernDiag(kname, xPlot, varargin{:}) - ...
   sum(KPlotData.*(KDataData\KPlotData')',2))) .* (A*normf);
[maxRMSPE,whBad] = max(RMSPE);
xBad = xPlot(whBad);
fBad = f(xBad);

whMiss = find((fPlot > fAppPlot + RMSPE + 1e5*eps) | ...
   (fPlot < fAppPlot - RMSPE - 1e5*eps));
Miss = [xPlot(whMiss) fPlot(whMiss) fAppPlot(whMiss) + [-1 1].*RMSPE(whMiss)];


end

function out = kernel(kname, x, y, varargin)
   dist = @(x,y) abs(x - y');
   kdist = @(z,theta) (1 + theta*z).*exp(-theta*z);
   if any(strcmp(kname{1}, {'Matern', 'MaternTh'}))
      out = kdist(dist(x,y),varargin{1});
   elseif strcmp(kname{1}, 'MaternMod')
      out = exp(varargin{2}*(x+y')).* kdist(dist(x,y),varargin{1});
   else
      disp('Invalid')
   end
end

function out = kernDiag(kname, x, varargin)
   if any(strcmp(kname{1}, {'Matern', 'MaternTh'}))
      out = ones(size(x,1),1);
   elseif strcmp(kname{1}, 'MaternMod')
      out = exp(varargin{2}*(2*x));
   else
      disp('Invalid')
   end
end

