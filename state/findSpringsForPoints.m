%given a set of points, return their incident springs.
function indexSprings = findSpringsForPoints(ss, indexPoints, func)
  aux = ismember(ss.springEnds, indexPoints);
  if nargin<3
    indexSprings = find(aux(:,1)|aux(:,2));
  else
    indexSprings = find(func(aux(:,1), aux(:,2)));
  end
end
