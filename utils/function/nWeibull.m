%normalized weibull (mode is 1)
function y = nWeibull(l,k,x)
C = (k-1)/k;
Xl = x/l;
y= ((Xl/(C^(1/k))).^(k-1)).*exp(C-Xl.^k);