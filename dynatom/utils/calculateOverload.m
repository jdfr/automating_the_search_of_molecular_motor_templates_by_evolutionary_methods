%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [overload ndatoms] = calculateOverload(evparams, ss)
%calculate statistics of superposition
[counts ndatoms] = evaluateSuperimposition(evparams, ss);
%discard first number which is the count of zeros
counts           = counts(2:end);
ncounts          = 1:numel(counts);
%global measure of superimposition
overload         = sum(counts.*exp(ncounts-1))/sum(counts);
%overload         = sum(counts.*(realpow(2, ncounts-1))/sum(counts);
%overload         = sum(counts.*(realpow(ncounts, 2))/sum(counts);
