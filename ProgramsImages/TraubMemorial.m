% function TraubMemorial
% TRAUBMEMORIAL create some figures for Traub memorial issue paper
% This script creates figures for 
%
% S.-C. T. Choi, Y. Ding, F. J. Hickernell, X. Tong, " ...

gail.InitializeDisplay %initialize the workspace and the display parameters

%% Sample functions with wildy oscillating second derivatives
f1 = @(x) x.^4 .* sin(1./((x==0)+x));
f1pp = @(x) (12*x.^2 - 1) .* sin(1./x) - 6*x .* cos(1./x);
f2 = @(x) f1(x) + 10.*x.^2;
f2pp = @(x) f1pp(x) + 20;
xplot = (-1:.001:1);


%% Sample functions with piecwise constant second derivatives
f3param = @(x,delta,c) (1./(2*delta.^2))*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
   - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta);
f3ppparam = @(x,delta,c) (1./(delta.^2))*(1 + sign(x-c-delta) - sign(x-c+delta)) ...
   .*(abs(x-c) <= 2*delta);

%% Sampling of hump functions for funappx_g
in_param.nlo = 10;
in_param.nhi = in_param.nlo;
in_param.abstol = 0.01;
in_param.a=-1;
in_param.b=1;
in_param.output_x = 1;
c = -0.2;
delta = 0.3;
f3 = @(x) f3param(x,delta,c);
gail.RemovePlotAxes
[~,fappxout] = funappx_g(@(x) -f3(x),in_param);
disp(['funappx_g used ' int2str(fappxout.npoints) ' points and ' ...
   int2str(fappxout.iter) ' iterations'])
h = plot(xplot,-f3(xplot),'-', ...
   fappxout.x, fappxout.y,'.');
axis([-1 1 -1. 0.3])
legend(h,{'\(f(x)\)','\((x_i,f(x_i))\)'}, ...
    'location', 'north','box','off', 'orientation','horizontal')
pbaspect([1 0.5 1])
print -depsc sampling-funappxg.eps

%% Sampling of hump function for funmin_g
gail.RemovePlotAxes
fminex = @(x) -f3(x) - 0.5*f3param(x,0.1,0.3);
[~,fminout] = funmin_g(fminex,in_param);
disp(['funmin_g used ' int2str(fminout.npoints) ' points and ' ...
   int2str(fminout.iter) ' iterations'])
h = plot(xplot,fminex(xplot),'-', ...
   fminout.x, fminout.y,'.');
axis([-1 1 -1. 0.3])
legend(h,{'\(f(x)\)','\((x_i,f(x_i))\)'}, ...
    'location', 'north','box','off', 'orientation','horizontal')
pbaspect([1 0.5 1])
print -depsc sampling-funming



