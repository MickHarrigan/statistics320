function rexp = randx(n,k,lambda)
%
%function rexp = randx(n,k,lambda)
% Generates samples of an exponentially distributed random variable
%
Z = rand(n,k); 
rexp=zeros(n,k);
rexp=-log(1-Z)/lambda;