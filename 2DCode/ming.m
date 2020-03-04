function value = ming(x) 
kernel = @(a,b,t,xinp) exp(WAddMatrix(t,xinp,b)).*(1+WDistanceMatrix(t,xinp,a))...
    .*exp(-WDistanceMatrix(t,xinp,a));
d = 2;
n = 5;
xlat = gail.lattice_gen(1,2^n,d);
%xdata = [xlat; [1,1]];
xdata  = xlat;
g = @(x1,x2) 1/6*((30+ 5* x1.*sin(5*x1)).* (4+exp(-5* x2))-100);
gData = g(xdata(:,1),xdata(:,2));
a = x(1:2);
b = x(3:4);
KDataData = kernel(a,b,xdata,xdata);
coeff = KDataData\gData;
value= 1/(2^n)*(log(det(KDataData)) + ...
                    log(gData'*coeff));
                
              