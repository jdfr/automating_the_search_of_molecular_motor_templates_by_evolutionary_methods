%return the number for which the distance from 1 to the Weibull cumulative
%function (multiplied by a given factor) becomes equal or smaller than a
%given tolerance 
function d = weibullMaxDist(l,k,factor,tol)
d=l*(-log(tol/factor))^(1/k);
