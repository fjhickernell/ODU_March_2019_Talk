%% Adaptive sampling and approximation example
clearvars
gail.InitializeDisplay
f = @(x) exp(-6*x).*sin(8*x+0.1) - 0.1;
%f = @(x) f0(x);
axisBox = [0 1 -0.3 0.5];

%% Plot function
xData = [0:0.1:0.6 0.8:0.1:1]';
fData = f(xData);
xPlot = (0:0.002:1)';
fPlot = f(xPlot);
kname = {'Matern','Afix',1}
theta = 1; %shape parameter
[KDataData, KPlotData, coeff, fAppPlot, RMSPE, xBad, fBad, ...
   maxRMSPE, whMiss, Miss, thetaOut] = ...
   RKHS(kname, f, xData, fData, xPlot, fPlot, theta);

meanf = mean(fPlot)


figure(1)
plot(xPlot,fPlot,xData,fData,'.')
xlabel('\(x\)')
ylabel('\(f(x)\)')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandData.eps')

%% Compute and plot approximation
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot);
set(h(1),'color',MATLABBlue)
set(h(2),'color',MATLABOrange)
set(h(3),'color',MATLABGreen)
xlabel('\(x\)')
[lgd,~] = legend(h,{'\(f\)','initial data', ...
   'APP\((\textsf{X},\textbf{\textit{y}})\)'}, ...
   'location', 'north','box','off', 'orientation','horizontal');
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppx.eps')

%% Next data point based on prediction error, no inference
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot, ...
   xPlot,fAppPlot + [-1,1].*RMSPE);
hold on
h = [h; scatter(xBad,fBad,200,MATLABPurple,'filled','d')];
set(h(1),'color',MATLABBlue)
set(h(2),'color',MATLABOrange)
set(h(3),'color',MATLABGreen)
set(h(4:5),'color',MATLABMaroon)
xlabel('\(x\)')
[~,icons] = legend(h([1:4 6]),{'\(f(x)\)','initial data', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x) \pm \)ERR\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   '\(\bigl(x_{11},f(x_{11})\bigr)\)'}, ...
   'box', 'off');
icons(14).Children.MarkerSize = 15;
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPE.eps')

whMiss, Miss


%% Infer theta using empirical Bayes
kname = {'MaternTh','Afix',1}
[KDataData, KPlotData, coeff, fOptAppPlot, RMSPEOpt, xBadOpt, fBadOpt, ...
   maxRMSPEOpt, whMissOpt, MissOpt, thetaOptOut] = ...
   RKHS(kname, f, xData, fData, xPlot, fPlot);

%% Compute and plot approximation with optimal theta
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptAppPlot);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','initial data','APP\((10)(x)\)'})
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndOptAppx.eps')

%% Next data point based on prediction error for optimal theta
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptAppPlot, ...
   xPlot,fOptAppPlot + [-1,1].*RMSPEOpt);
hold on
h = [h; scatter(xBadOpt,fBadOpt,200,MATLABPurple,'filled','d')];
set(h(1),'color',MATLABBlue)
set(h(2),'color',MATLABOrange)
set(h(3),'color',MATLABGreen)
set(h(4:5),'color',MATLABMaroon)
xlabel('\(x\)')
[~,icons] = legend(h([1:4 6]),{'\(f(x)\)','initial data', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x) \pm \)ERR\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   '\(\bigl(x_{11},f(x_{11})\bigr)\)'}, ...
   'box', 'off');
%lgd.NumColumns = 2;
icons(14).Children.MarkerSize = 15;
%legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPEOpt.eps')

whMissOpt, MissOpt

%% Infer y-varying kernel using empirical Bayes
kname = {'MaternMod','Afix',1}
[KDataData, KPlotData, coeff, fOptyAppPlot, RMSPEOpty, xBadOpty, fBadOpty, ...
   maxRMSPEOpty, whMissOpty, MissOpty, thetaOptyOut] = ...
   RKHS(kname, f, xData, fData, xPlot, fPlot);

%% Compute and plot approximation with optimal y-varying
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptyAppPlot);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','initial data','APP\((n)(x)\)'})
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndOptyAppx.eps')

%% Next data point based on prediction error for y-varying
n = length(xData);
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptyAppPlot, ...
   xPlot,fOptyAppPlot + [-1,1].*RMSPEOpty);
hold on
h = [h; scatter(xBadOpty,fBadOpty,200,MATLABPurple,'filled','d')];
set(h(1),'color',MATLABBlue)
set(h(2),'color',MATLABOrange)
set(h(3),'color',MATLABGreen)
set(h(4:5),'color',MATLABMaroon)
xlabel('\(x\)')
[~,icons] = legend(h([1:4 6]),{'\(f(x)\)','initial data', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x) \pm \)ERR\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   '\(\bigl(x_{11},f(x_{11})\bigr)\)'}, ...
   'box', 'off');
%lgd.NumColumns = 2;
icons(14).Children.MarkerSize = 15;
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPEOpty.eps')

whMissOpty, MissOpty


%% Finish iterating until tolerance is reached
n0 = length(xData);
tol = 0.02; %tolerance
thetaOptyOutVec = [];
while maxRMSPEOpty > tol
   xData = [xData; xBadOpty];
   fData = [fData; fBadOpty];
   thetaOptyOutVec = [thetaOptyOutVec; thetaOptyOut];
   [KDataData, KPlotData, coeff, fOptyAppPlot, RMSPEOpty, xBadOpty, fBadOpty, ...
      maxRMSPEOpty, whMissOpty, MissOpty, thetaOptyOut] = ...
      RKHS(kname, f, xData, fData, xPlot, fPlot, theta);
end
thetaOptyOutVec = [thetaOptyOutVec; thetaOptyOut];

%% Plot final result
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData(1:n0),fData(1:n0),'.',xPlot,fOptyAppPlot, ...
   xPlot,fOptyAppPlot + [-1,1].*RMSPEOpty);
hold on
h = [h; scatter(xData(n0+1:end),fData(n0+1:end),200,MATLABPurple,'filled','d')];
set(h(1),'color',MATLABBlue)
set(h(2),'color',MATLABOrange)
set(h(3),'color',MATLABGreen)
set(h(4:5),'color',MATLABMaroon)
xlabel('\(x\)')
[~,icons] = legend(h([1:4 6]),{'\(f(x)\)','initial data', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   'APP\((\textsf{X},\textbf{\textit{y}})(x) \pm \)ERR\((\textsf{X},\textbf{\textit{y}})(x)\)', ...
   ['additional data for \(\varepsilon = ' num2str(tol) '\)']}, ...
   'box', 'off');
%lgd.NumColumns = 2;
icons(14).Children.MarkerSize = 15;
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc',['fandDataAndAppxAndRMSPEOptyFinal' int2str(1000*tol) '.eps'])

whMissOpty, MissOpty
