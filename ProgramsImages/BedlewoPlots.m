%% Bedlewo Plots
gail.InitializeWorkspaceDisplay %clean up 
format long

%%
theta = 0.7;
ltblue = (1-theta)*MATLABBlue + theta*ones(1,3);
theta = 0.75;
LtMaroon = (1-theta)*MATLABMaroon + theta*ones(1,3);
theta = 0.5;
MedMaroon = (1-theta)*MATLABMaroon + theta*ones(1,3);

%%
f = @(x) x./sqrt(1 + x.^2);
fpp = @(x) -3.*x.*((1 + x.^2).^(-2.5));
%xnodes = [0 0.2 0.4 0.6 0.8 1.2];
%xnodes = [0:0.4:1.6 2.4];
xnodes = [0:0.5:2 3];
xxnodes = [2*xnodes(1) - xnodes(2) xnodes 2*xnodes(end) - xnodes(end-1)];
fnodes = f(xnodes);
xfnodes = f(xxnodes);
nnodes = numel(xnodes);
noffset = 43;
extra = 0.15;
xmin = min(xnodes) - extra;
xmax = max(xnodes) + extra;
xplot = xmin + (0:1e-3:1)*(xmax - xmin);
fplot = f(xplot);
fppplot = fpp(xplot);
fppnodes = fpp(xnodes);
fmin = min(fplot);
fmax = max(fplot);
ylabelval = fmin - 0.03;

%% Plot function and spline
gail.RemovePlotAxes
h = plot(xplot,fplot,'color',MATLABBlue);
axis([xmin xmax fmin fmax+0.3])
plot([xnodes; xnodes],[fmin*ones(size(xnodes)); fnodes],'--k')
h = [h; plot(xxnodes,xfnodes,'color',MATLABGreen)];
h = [h; plot(xnodes,fnodes,'.','color',MATLABOrange)];
for ii = 1:nnodes
   text(xnodes(ii)-0.05,ylabelval-0.02,['\(x_{' int2str(noffset + ii) '}\)'],'color',zeros(1,3))
end
pbaspect([1 0.5 1])
legend(h([1 3 2]),{'\(f\)','data','APP\((n,\textsf{X},\textbf{\textit{y}})\)'},...
   'location','northwest','orientation','horizontal')
legend boxoff
print -depsc LinearSpline.eps

%% Plot second derivative too
fdd1 = diff(diff(fnodes(1:3))./diff(xnodes(1:3)))/diff(xnodes([1 3]));
fdd2 = diff(diff(fnodes(4:6))./diff(xnodes(4:6)))/diff(xnodes([4 6]));
gail.RemovePlotAxes
plot(xplot,fplot,'color',MATLABBlue);
axis([xmin xmax fmin fmax+0.3])
plot([xnodes; xnodes],[fmin*ones(size(xnodes)); fnodes],'--k')
plot(xxnodes,xfnodes,'color',MATLABGreen);
plot(xnodes,fnodes,'.','color',MATLABOrange);
h = plot(xplot,abs(fppplot),'-','color',MATLABPurple);
h = [h; plot(xnodes([1 3]),2*abs(fdd1)*ones(1,2), ...
   xnodes([4 6]),2*abs(fdd2)*ones(1,2), ...
   'color',MedMaroon)];
h = [h; plot(xnodes([1 3]),abs(fppnodes(1))*ones(1,2), ...
   xnodes([4 6]),abs(fppnodes(6))*ones(1,2), ...
   'color',MATLABMaroon)];
h = [h; plot(xnodes([3 4]),abs(fppnodes(3))*ones(1,2), ...
   xnodes([1 3]),abs(fpp(0.5))*ones(1,2), ...
   xnodes([4 6]),abs(fppnodes(4))*ones(1,2), ...
   'color',LtMaroon)];

for ii = 1:nnodes
   text(xnodes(ii)-0.05,ylabelval-0.02,['\(x_{' int2str(noffset + ii) '}\)'],'color',zeros(1,3))
end
pbaspect([1 0.5 1])
legend(h([1 6 2 4]),{'\(|f''''|\)','\(||f''''||_{\infty}\)','\(2|f[x_{i-1},x_i,x_{i+1}]|\)','\(||f''''||_{-\infty}\)'},...
   'location','northwest','orientation','horizontal')
legend boxoff
print -depsc LinearSplineSecDeriv.eps

