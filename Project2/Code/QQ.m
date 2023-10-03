function q=QQ(x)
%function q=QQ(x)
%  Performs the Gaussian error function, with mean 0 and unity variance
%    Q(x)=integral from x to infinity of exp(-t^2/2) dt
%
%    Method uses the MATLAB erfc
%       QQ(x)=0.5*erfc(x/sqrt(2))
%
q=0.5*erfc(x/sqrt(2));

