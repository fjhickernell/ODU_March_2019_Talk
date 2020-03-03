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
s = 1; %scale parameter
theta = 1; %shape parameter
dist = @(x,y) abs(x - y');
kernel = @(z,s,theta) s*(1 + theta*z).*exp(-theta*z);
KDataData = kernel(dist(xData,xData),s,theta);
coeff = KDataData\fData;
KPlotData = kernel(dist(xPlot,xData),s,theta);
fAppPlot = KPlotData*coeff;
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot);
set(h(1),'color',MATLABBlue)
set(h(2),'color',MATLABOrange)
set(h(3),'color',MATLABGreen)
xlabel('\(x\)')
[lgd,icons] = legend(h,{'\(f\)','data', ...
   'APP\((\textsf{X},\textbf{\textit{y}})\)'}, ...
   'location', 'north','box','off', 'orientation','horizontal');
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppx.eps')

% %% Compute and plot approximation also using small design
% xSmallData = [0 0.4 0.6 1]';
% fSmallData = f(xSmallData);
% KSmallDataData = kernel(dist(xSmallData,xSmallData),s,theta);
% coeffSmall = KSmallDataData\fSmallData;
% KPlotSmallData = kernel(dist(xPlot,xSmallData),s,theta);
% fAppPlotSmall = KPlotSmallData*coeffSmall;
% figure
% h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot,xPlot,fAppPlotSmall);
% xlabel('\(x\)')
% legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((10)(x)\)', 'APP\((4)(x)\)'})
% legend('boxoff')
% axis(axisBox)
% set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
% pos = get(gcf,'Position');
% set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
% print('-depsc','fandDataAndAppxSmall.eps')

%% Next data point based on prediction error, no inference
A = 1;
normf = sqrt(coeff'*fData);
RMSPE = real(sqrt(kernel(0,s,theta) - ...
   sum(KPlotData.*(KDataData\KPlotData')',2))) .* normf;
[~,whBad] = max(RMSPE);
xBad = xPlot(whBad);
fBad = f(xBad);
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot, ...
   xPlot,fAppPlot + A*[-1,1].*RMSPE);
hold on
h = [h; scatter(xBad,fBad,200,MATLABPurple,'filled','d')];
set(h(1),'color',MATLABBlue)
set(h(2),'color',MATLABOrange)
set(h(3),'color',MATLABGreen)
set(h(4:5),'color',MATLABMaroon)
xlabel('\(x\)')
[~,icons] = legend(h([1:4 6]),{'\(f(x)\)','\(f(x_i)\)', ...
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

whMiss = find((fPlot > fAppPlot + A*RMSPE + 1000*eps) | ...
   (fPlot < fAppPlot - A*RMSPE - 1000*eps))
Miss = [xPlot(whMiss) fPlot(whMiss) fAppPlot(whMiss) + A*[-1 1].*RMSPE(whMiss)]

% %% Find next data point for optimization
% [~,whBadMin] = min(fAppPlot - A.*RMSPE);
% xBadMin = xPlot(whBadMin);
% fBadMin = f(xBadMin);
% [fDataMin,whDataMin] = min(fData);
% xDataMin = xData(whDataMin);
% figure
% h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot, ...
%    xPlot,fAppPlot + A*[-1,1].*RMSPE);
% hold on
% h = [h; scatter(xBad,fBad,200,MATLABPurple,'filled','d'); ...
%    scatter(xDataMin,fDataMin,800,MATLABGreen,'filled','p'); ...
%    scatter(xBadMin,fBadMin,200,MATLABCyan,'filled','s')];
% set(h(4:5),'color',MATLABMaroon)
% xlabel('\(x\)')
% [~,icons,plts,txt] = legend(h([1:3 7 4 6 8]), ...
%    {'\(f(x)\)','\(f(x_i)\)',...
%    'APP\((10)(x)\)', ...
%    '\(\bigl(\widehat{x}_{\mathrm{MIN}},\widehat{\mathrm{MIN}}(10) \bigr)\)', ... 
%    'APP\((10)(x) \pm \)SVU\((10)(x)\)', ...
%    '\(\bigl(x_{\mathrm{ID},11},f(x_{\mathrm{ID},11})\bigr)\)', ...
%    '\(\bigl(x_{\mathrm{MIN},11},f(x_{\mathrm{MIN}, 11})\bigr)\)'});
% legend('boxoff')
% icons(17).Children.MarkerSize = 15;
% icons(14).Children.MarkerSize = 20;
% icons(18).Children.MarkerSize = 15;
% %lgd.NumColumns = 2;
% axis(axisBox)
% set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
% pos = get(gcf,'Position');
% set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
% print('-depsc','fandDataAndAppxAndRMSPEAndMin.eps')


%% Infer theta using empirical Bayes
AOpt0 = 5;
Ktheta = @(logth) kernel(dist(xData,xData),s,exp(logth));
objective = @(K,y) mean(log(eig(K))) + log(y'*(K\y));
figure
logthetaRange = log(0.01):0.005:log(100);
semilogx(exp(logthetaRange), ...
   arrayfun(@(lgth) objective(Ktheta(lgth),fData), logthetaRange))
logthopt = fminbnd(@(logth) objective(Ktheta(logth),fData),-5,5);
thetaopt = exp(logthopt)

%% Compute and plot approximation with optimal theta
theta = thetaopt;
KOptDataData = kernel(dist(xData,xData),s,theta);
coeffOpt = KOptDataData\fData;
KOptPlotData = kernel(dist(xPlot,xData),s,theta);
fOptAppPlot = KOptPlotData*coeffOpt;
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptAppPlot);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((10)(x)\)'})
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndOptAppx.eps')

%% Next data point based on prediction error for optimal theta
%AOpt = 1;
normf = sqrt(coeffOpt'*fData);
kernelDiag = @(s,theta) s;
fxOptnorm = real(sqrt(kernelDiag(s,theta) - ...
   sum(KOptPlotData.*(KOptDataData\KOptPlotData')',2)));
AOpt = AOpt0 * max(fxOptnorm./sqrt(kernelDiag(s,theta)))
RMSPEOpt = fxOptnorm .* (AOpt * normf);
[~,whBadOpt] = max(RMSPEOpt);
xBadOpt = xPlot(whBadOpt);
fBadOpt = f(xBadOpt);
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
[~,icons] = legend(h([1:4 6]),{'\(f(x)\)','\(f(x_i)\)', ...
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

whMissOpt = find((fPlot > fOptAppPlot + AOpt*RMSPEOpt + 1000*eps) | ...
   (fPlot < fOptAppPlot - AOpt*RMSPEOpt - 1000*eps))
Miss = [xPlot(whMissOpt) fPlot(whMissOpt) fOptAppPlot(whMissOpt) + ...
   AOpt*[-1 1].*RMSPEOpt(whMissOpt)]

%% Infer y-varying kernel using empirical Bayes
kernel = @(z,s,theta,a,x,y) s*exp(a*(x+y')).*(1 + theta*z).*exp(-theta*z);
Ktheta = @(logth,a) kernel(dist(xData,xData),s,exp(logth),a,xData,xData);
objective = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
% theta = 1;
% aRange = -15:0.01:0;
% plot(aRange, ...
%    arrayfun(@(a) objective(Ktheta(log(theta),a),fData), aRange))
% aopt = fminbnd(@(a) objective(Ktheta(log(theta),a),fData),-15,0)
[bopt,objmin] = fminsearch(@(b) objective(Ktheta(b(1),b(2)),fData),[2,-10])
thetaopt = exp(bopt(1))
aopt = bopt(2)

%% Compute and plot approximation with optimal y-varying
theta = thetaopt;
a = aopt;
KOptyDataData = kernel(dist(xData,xData),s,theta,a,xData,xData);
coeffOpty = KOptyDataData\fData;
KOptyPlotData = kernel(dist(xPlot,xData),s,theta,a,xPlot,xData);
fOptyAppPlot = KOptyPlotData*coeffOpty;
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptyAppPlot);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((n)(x)\)'})
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndOptyAppx.eps')

%% Next data point based on prediction error for y-varying
n = length(xData);
%AOpty = -tinv(0.005,n)/sqrt(length(xData))
%AOpty = 1;
normf = sqrt(coeffOpty'*fData);
kernelDiag = @(s,theta,a,x) s*exp(a*(2*x));
fxOptnormy = real(sqrt(kernelDiag(s,theta,a,xPlot) - ...
   sum(KOptyPlotData.*(KOptyDataData\KOptyPlotData')',2)));
AOpty = AOpt0 * max(fxOptnormy./sqrt(kernelDiag(s,theta,a,xPlot)))
RMSPEOpty = fxOptnormy .* (AOpty * normf);
[maxRSMEOpty,whBadOpty] = max(RMSPEOpty);
xBadOpty = xPlot(whBadOpty);
fBadOpty = f(xBadOpty);
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
[~,icons] = legend(h([1:4 6]),{'\(f(x)\)','\(f(x_i)\)', ...
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

whMissOpty = find((fPlot > fOptyAppPlot + AOpty*RMSPEOpty + 1000*eps) | ...
   (fPlot < fOptyAppPlot - AOpty*RMSPEOpty - 1000*eps))
Miss = [xPlot(whMissOpty) fPlot(whMissOpty) fOptyAppPlot(whMissOpty) + ...
   AOpty*[-1 1].*RMSPEOpty(whMissOpty)]

% %% Find next data point for optimization based on prediction error for y-varying
% [~,whBadOptyMin] = min(fAppPlot - A.*RMSPEOpty);
% xBadOptyMin = xPlot(whBadOptyMin);
% fBadOptyMin = f(xBadOptyMin);
% figure
% h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptyAppPlot, ...
%    xPlot,fOptyAppPlot + AOpty*[-1,1].*RMSPEOpty);
% hold on
% h = [h; scatter(xBadOpty,fBadOpty,200,MATLABPurple,'filled','d'); ...
%    scatter(xDataMin,fDataMin,800,MATLABGreen,'filled','p'); ...
%    scatter(xBadOptyMin,fBadOptyMin,200,MATLABCyan,'filled','s')];
% set(h(1),'color',MATLABBlue)
% set(h(2),'color',MATLABOrange)
% set(h(3),'color',MATLABGreen)
% set(h(4:5),'color',MATLABMaroon)
% xlabel('\(x\)')
% [~,icons] = legend(h([1:3 7 4 6 8]),...
%    {'\(f(x)\)','\(f(x_i)\)','APP\((10)(x)\)', ...
%     '\(\bigl(\widehat{x}_{\mathrm{MIN}},\widehat{\mathrm{MIN}}(10) \bigr)\)', ... 
%   'APP\((10)(x) \pm \)SVU\((10)(x)\)', ...
%    '\(\bigl(x_{\mathrm{ID},11},f(x_{\mathrm{ID},11})\bigr)\)', ...
%    '\(\bigl(x_{\mathrm{MIN},11},f(x_{\mathrm{MIN},11})\bigr)\)'}, ...
%    'box','off');
% %lgd.NumColumns = 2;
% %legend('boxoff')
% icons(17).Children.MarkerSize = 15;
% icons(14).Children.MarkerSize = 20;
% icons(18).Children.MarkerSize = 15;
% axis(axisBox)
% set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
% pos = get(gcf,'Position');
% set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
% print('-depsc','fandDataAndAppxAndRMSPEOptyAndMin.eps')
% 


%% Finish iterating until tolerance is reached
n0 = length(xData);
tol = 0.02; %tolerance
err = AOpty * maxRSMEOpty; %error bound
while err > tol
   [bopt,objmin] = fminsearch(@(b) objective(Ktheta(b(1),b(2)),fData),[2,-10])
   thetaopt = exp(bopt(1))
   aopt = bopt(2)

   %% Compute and plot approximation with optimal y-varying
   theta = thetaopt;
   a = aopt;
   KOptyDataData = kernel(dist(xData,xData),s,theta,a,xData,xData);
   coeffOpty = KOptyDataData\fData;
   KOptyPlotData = kernel(dist(xPlot,xData),s,theta,a,xPlot,xData);
   fOptyAppPlot = KOptyPlotData*coeffOpty;

   %% Next data point based on prediction error for y-varying
   n = length(xData);
   %AOpty = 1;
   normf = sqrt(coeffOpty'*fData);
   kernelDiag = @(s,theta,a,x) s*exp(a*(2*x));
   fxOptnormy = real(sqrt(kernelDiag(s,theta,a,xPlot) - ...
      sum(KOptyPlotData.*(KOptyDataData\KOptyPlotData')',2)));
   AOpty = AOpt0 * max(fxOptnormy./sqrt(kernelDiag(s,theta,a,xPlot)))
   RMSPEOpty = fxOptnormy .* (AOpty * normf);
   [maxRSMEOpty,whBadOpty] = max(RMSPEOpty);
   err = AOpty * maxRSMEOpty; %error bound
   xBadOpty = xPlot(whBadOpty)
   fBadOpty = f(xBadOpty);
   xData = [xData; xBadOpty];
   fData = [fData; fBadOpty];
   Ktheta = @(logth,a) kernel(dist(xData,xData),s,exp(logth),a,xData,xData);
end

%% Plot final result
gail.RemovePlotAxes
h = plot(xPlot,fPlot,xData(1:n0),fData(1:n0),'.',xPlot,fOptyAppPlot, ...
   xPlot,fOptyAppPlot + AOpty*[-1,1].*RMSPEOpty);
hold on
h = [h; scatter(xData(n0+1:end-1),fData(n0+1:end-1),200,MATLABPurple,'filled','d')];
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

whMissOpt = find((fPlot > fOptAppPlot + AOpt*RMSPEOpt + 1000*eps) | ...
   (fPlot < fOptAppPlot - AOpt*RMSPEOpt - 1000*eps))
Miss = [xPlot(whMissOpt) fPlot(whMissOpt) fOptAppPlot(whMissOpt) + ...
   AOpt*[-1 1].*RMSPEOpt(whMissOpt)]
