
  function DM = WAddMatrix(dsites,ctrs,b)
  M = size(dsites,1); N = size(ctrs,1);
% Algorithm is based on expanding the terms and computing each term
% explicitly, i.e.  
%         (x1 - x2)^2 = x1.^2 + x2.^2 - 2*x1*x2;
  DM = zeros(M,N);
  for i = 1:M
      temp = (dsites(i,:)+ ctrs)*b;
      DM(i,:) = temp';
  end
 