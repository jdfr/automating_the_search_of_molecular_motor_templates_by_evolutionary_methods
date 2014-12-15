%to be used before simulating if the number of dimensions, points and
%springs has changed
function ss= makeSimulationAmounts(ss)
  ss.ndims = size(ss.pos,2);
  ss.npoints = size(ss.pos,1);
  ss.nsprings = size(ss.k,1);
end
