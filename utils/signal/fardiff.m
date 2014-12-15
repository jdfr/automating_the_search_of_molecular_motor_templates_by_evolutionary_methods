function dxn = fardiff(x, n)

dxn = x(1:(end-n))-x((n+1):end);
%dxn=-diff(x);

