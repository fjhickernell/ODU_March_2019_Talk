%% Adaptive sampling and approximation example
close all
clearvars
gail.InitializeDisplay

%% Cheng's function
d = 2;
f = @(x) cos(sum(x,2)).*exp(prod(x,2));
%f = @(x) cos(sum(x,2)).*exp(prod(1-x,2));
m0 = 4;
delta = 0.1;
xlat = 1 - (gail.lattice_gen(1,2^m0,d) + 2^(-m0-1));
xData = xlat;
fData = f(xData);
xx = (0:0.005:1)';
[XX, YY] = meshgrid(xx);
xPlot = [XX(:) YY(:)];
fPlot = f(xPlot);
ZZ = reshape(fPlot,size(XX));
gail.Remove3DPlotAxes
surf(XX,YY,ZZ,'FaceColor','Interp','EdgeColor','None')
hold on
plot3(xData(:,1), xData(:,2), fData + delta, '.', 'color', MATLABOrange)
colorbar
%title('\(f(\textbf{\textit{y}}) = \cos(x_1 + x_2) \exp(x_1 * x_2)\)')
axis equal
squeezeWhite
print('-depsc','ChengFunData.eps')


%% No inference
kname = {'Matern','Afix',1}
theta = ones(1,d); %shape parameter
[KDataData, KPlotData, coeff, fAppPlot, RMSPE, xBad, fBad, ...
   maxRMSPE, whMiss, Miss, thetaOut] = ...
   RKHS(kname, f, xData, fData, xPlot, fPlot, theta);

gail.Remove3DPlotAxes
EE = reshape(RMSPE,size(XX));
surf(XX,YY,EE,'FaceColor','Interp','EdgeColor','None')
hold on
plot3(xData(:,1), xData(:,2), maxRMSPE*ones(size(xData,1),1), '.', 'color', MATLABOrange)
scatter3(xBad(:,1), xBad(:,2), maxRMSPE+delta, 200,MATLABPurple,'filled','d')
colorbar
axis equal
squeezeWhite
print('-depsc','ChengFunErrNoOpt.eps')
clax = caxis; %want to keep this color axis
whMiss, Miss



%% Inference with mod Matern
kname = {'MaternMod','Afix',1}


[KDataData, KPlotData, coeff, fAppPlot, RMSPEOpty, xBadOpty, fBadOpty, ...
   maxRMSPEOpty, whMissOpty, MissOpty, thetaOptyOut] = ...
   RKHS(kname, f, xData, fData, xPlot, fPlot, theta);

gail.Remove3DPlotAxes
EE = reshape(RMSPEOpty,size(XX));
surf(XX,YY,EE,'FaceColor','Interp','EdgeColor','None')
caxis(clax);
hold on
plot3(xData(:,1), xData(:,2), maxRMSPEOpty*ones(size(xData,1),1), '.', 'color', MATLABOrange)
scatter3(xBad(:,1), xBad(:,2), maxRMSPEOpty+delta, 200,MATLABPurple,'filled','d')
colorbar
axis equal
squeezeWhite
print('-depsc','ChengFunErrOpty.eps')

whMissOpty, MissOpty

return

%% Finish iterating
n0 = length(xData);
tol = 0.05; %tolerance
thetaOptyOutVec = [];
while maxRMSPEOpty > tol
   xData = [xData; xBadOpty];
   fData = [fData; fBadOpty];
   thetaOptyOutVec = [thetaOptyOutVec; thetaOptyOut];
   [KDataData, KPlotData, coeff, fOptyAppPlot, RMSPEOpty, xBadOpty, fBadOpty, ...
      maxRMSPEOpty, whMissOpty, MissOpty, thetaOptyOut] = ...
      RKHS(kname, f, xData, fData, xPlot, fPlot, theta);
   maxRMSPEOpty
end
thetaOptyOutVec = [thetaOptyOutVec; thetaOptyOut];
whMissOpty, MissOpty

%% Plot final results
gail.Remove3DPlotAxes
n1 = size(xData,1);
EE = reshape(RMSPEOpty,size(XX));
surf(XX,YY,EE,'FaceColor','Interp','EdgeColor','None')
hold on
plot3(xData(1:n0,1), xData(1:n0,2), maxRMSPEOpty*ones(n0,1), '.', 'color', MATLABOrange)
scatter3(xData(n0+1:n1,1), xData(n0+1:n1,2), maxRMSPEOpty*ones(n1-n0,1), 200,MATLABPurple,'filled','d')
colorbar
axis equal
squeezeWhite
print('-depsc','ChengFunErrOptyFinal.eps')

