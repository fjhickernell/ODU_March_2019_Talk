function value = minfun(x) 
kernel = @(a,b,t,xinp) exp(WAddMatrix(t,xinp,b)).*(1+WDistanceMatrix(t,xinp,a))...
    .*exp(-WDistanceMatrix(t,xinp,a));
f = @(x1,x2) cos(x1 + x2).*exp(x1.*x2);
d = 2;
n = 5;
xlat = gail.lattice_gen(1,2^n,d);
%xdata = [xlat; [1,1]];
fData = f(xlat(:,1),xlat(:,2));
a = x(1:2);
b = x(3:4);
KDataData = kernel(a,b,xlat,xlat);
coeff = KDataData\fData;
value= 1/(2^n)*(log(det(KDataData)) + ...
                    log(fData'*coeff));
                
              