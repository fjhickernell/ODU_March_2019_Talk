%% Adaptive sampling and approximation example
kernel = @(a,b,t,x) exp(WAddMatrix(t,x,b)).*(1+WDistanceMatrix(t,x,a))...
    .*exp(-WDistanceMatrix(t,x,a));

f = @(x1,x2) cos(x1 + x2).*exp(x1.*x2);
d = 2;
n = 5;
xlat = gail.lattice_gen(1,2^n,d);
%xdata = [xlat; [1,1]];
fData = f(xlat(:,1),xlat(:,2));
%find optimal shape parameter
%if change function f, need to modify in minfun
optimv = fminsearch(@minfun,[2;2;-2;-2]);

fprintf('The optimal shape parameter is %.4f\n', optimv(1:2));
fprintf('The optimal stretch parameter is %.4f\n', optimv(3:4));
optima = optimv(1:2);
optimb = optimv(3:4);
KDataData = kernel(optima,optimb,xlat,xlat);
    coeff = KDataData\fData;
xPlot = (0:0.002:1)';
[XPlot, YPlot] = meshgrid(xPlot,xPlot);
xplot = XPlot(:);
yplot = YPlot(:);
fPlot = f(xplot,yplot);
KPlotData = kernel(optima,optimb,[xplot, yplot],xlat);
 fAppPlot  = KPlotData*coeff;
  yvalue  = fData'*coeff;
temp = kernel(optima, optimb,[xPlot, xPlot],[xPlot, xPlot]);
M = sum(KPlotData.*(KDataData\KPlotData')',2);
value = (temp(:)-M).*yvalue;
%index = find(value>0.995*max(value));
index = find(value == max(value));
disp('The sample should be:  ');
disp([xplot(index), yplot(index)]);

error = abs(fPlot-fAppPlot);
errorplot = reshape(error,length(xPlot),length(xPlot));

%Plot the graph
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24, ... %make font larger
      'defaultLineLineWidth',3, ... %thick lines
      'defaultLineMarkerSize',18) %big dots
figure;
mesh(XPlot,YPlot,errorplot,errorplot);
colorbar;
figure;
contour(XPlot,YPlot,errorplot);
colorbar;
hold on;
plot(xlat(:,1),xlat(:,2),'r*');
hold on;
plot(xplot(index),yplot(index),'ko')
legend('error', 'datasites', 'nextsample')

%% Using same datasites with a different function
g = @(x1,x2) 1/6*((30+ 5* x1.*sin(5*x1)).* (4+exp(-5* x2))-100);
optimv = fminsearch(@ming,[1;1;-1;-1]);
optimag = optimv(1:2);
optimbg = optimv(3:4);
fprintf('The optimal shape parameter is %.4f\n', optimag);
fprintf('The optimal stretch parameter is %.4f\n', optimbg);
KDataData = kernel(optimag,optimbg,xlat,xlat);
gData = g(xlat(:,1),xlat(:,2));
coeff = KDataData\gData;
gPlot = g(xplot,yplot);
KPlotData = kernel(optimag,optimbg,[xplot, yplot],xlat);
 gAppPlot  = KPlotData*coeff;
  yvalue  = gData'*coeff;
temp = kernel(optimag, optimbg,[xPlot, xPlot],[xPlot, xPlot]);
M = sum(KPlotData.*(KDataData\KPlotData')',2);
value = (temp(:)-M).*yvalue;
%indexg = find(value>0.995*max(value));
indexg = find(value == max(value));
disp('The sample should be:  ');
disp([xplot(indexg), yplot(indexg)]);
error = abs(gPlot-gAppPlot);
errorplot = reshape(error,length(xPlot),length(xPlot));

%Plot the graph
figure;
mesh(XPlot,YPlot,errorplot,errorplot);
colorbar;
figure;
contour(XPlot,YPlot,errorplot);
colorbar;
hold on;
plot(xlat(:,1),xlat(:,2),'r*');
hold on;
plot(xplot(indexg),yplot(indexg),'ko')
legend('error', 'datasites', 'nextsample')