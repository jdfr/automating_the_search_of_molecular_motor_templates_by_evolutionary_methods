function Y = pdist1D(X)
%customized version of PDIST for 1D data, faster than pdist, of course
%X is a column vector of n elements, Y is a column vector containing
%pairwise distances, ordered

n = length(X);
Y = zeros(1,n*(n-1)./2, class(X));

k = 1;
for i = 1:n-1
  Y(k:(k+n-i-1)) = X((i+1):n,:)-X(i,:);
  k = k + (n-i);
end
Y = abs(Y);
