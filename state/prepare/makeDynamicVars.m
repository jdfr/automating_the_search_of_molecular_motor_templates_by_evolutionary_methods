%to be used before simulating if the number of points or spring changes, or
%any flag in the dynamic_X vectors change 
function ss = makeDynamicVars(ss)
  %for each vector of dynflags, calculate amount of trues and cumsums
  for k=1:numel(ss.flagVars)
    campo = ss.flagVars{k};
    ss.(ss.numFlagVars{k})    = sum(ss.(campo));
    ss.(ss.cumsumFlagVars{k}) = cumsum(ss.(campo));
  end
  %selectors for positional variables X,Y,Z,...
  ss.selector.pos = [1 ss.numdyn_p*ss.ndims];
  %selectors for velocity variables dX,dY,dZ,...
  ss.selector.vel = ss.selector.pos(end) + ss.selector.pos;
  %selectors for other dynamical variables
  offset = ss.selector.vel(end);
  for k=3:numel(ss.dynVars) %skipping the first two: pos and vel
    ss.selector.(ss.dynVars{k}) = offset + [1 ss.(ss.doubledNumFlagVars{k})];
    if (ss.(ss.doubledNumFlagVars{k})~=0);
      offset = ss.selector.(ss.dynVars{k})(end);
    end
  end
end