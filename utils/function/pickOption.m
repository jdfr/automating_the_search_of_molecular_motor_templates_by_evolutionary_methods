function choices = pickOption(options, weights, n)
%given a table of relative probabilities, pick options according to the
%probabilities

if nargin<3
  n=1;
end

if (numel(options)~=numel(weights)) && (numel(weights)>0)
  error('The amount of options must match with the amount of weights (%g, %g)!!!', numel(options), numel(weights));
end

rnds = rand(n,1);

if isempty(weights)
  choices = options(ceil(rnds*numel(options)));
  return
end

%wipe null options out
notnulls = weights~=0;
if ~all(notnulls)
  weights = weights(notnulls);
  options = options(notnulls);
end

%cumulative distribution out of weights
cumulative = cumsum(weights)/sum(weights);

% if numel(rnds)==1
%   if numel(find(rnds<cumulative, 1, 'first'))~=1
%     rnds=rnds;
%   end
% end

if n==1
  choices = options(find(rnds<cumulative, 1));%IMPLIED:, 'first');
else
  choices = options(arrayfun(@(rnd) find(rnd<cumulative, 1), rnds)); %IMPLIED: , 'first'
end
