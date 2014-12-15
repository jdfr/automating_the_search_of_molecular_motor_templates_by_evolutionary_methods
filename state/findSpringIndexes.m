%given a set of springs specified by their ends, pinpoint their indexes in
%the state ss
function [indexes members] = findSpringIndexes(ss, springs, doNotSort)
if (nargin>2) && doNotSort
  [members indexes] = ismember(springs, ss.springsEnds, 'rows');
else
  [members indexes] = ismember(sort(springs,2), sort(ss.springEnds,2), 'rows');
end
