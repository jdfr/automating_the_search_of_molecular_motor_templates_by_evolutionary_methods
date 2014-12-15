%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%given a precision (like 0.001), rounds the number to conform to that
%precision
function n = capPrecision(n, p)
m = mod(n, p);
n = n-m+round(m/p)*p;
end
