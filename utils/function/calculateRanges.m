function [ranges rangesizes] = calculateRanges(popsize, numRanges)
%given a number of objects, and a number of ranges, create indexes to distribute evenly the objects in the given number of ranges
rangesizes  = diff(floor(linspace(0, popsize, numRanges+1)));
rangesizes  = rangesizes(rangesizes~=0);
ranges      = [0; cumsum(rangesizes)'];
ranges      = [ranges(1:(end-1))+1, ranges(2:end)];
