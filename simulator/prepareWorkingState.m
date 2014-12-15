%dump 'state' into 'ss' in order to have properly arranged
%vars. Again, this function may be a very simple and generic loop, but
%it is critical to improve performance, and using dynamic fields to
%access to properties is very slow
function [ss dynamicVel] = prepareWorkingState(ss, state)
  switch ss.ndims
    case 2
      np = ss.npoints;
      dynamicVel = [state((2*np+1):(3*np)), state((3*np+1):(4*np))];
      ss.pos(ss.dynamical_p,:) = [state(1:np), state((np+1):(2*np))];
      ss.vel(ss.dynamical_p,:) = dynamicVel;
    otherwise
      dynamicVel = reshape(state(ss.selector.vel(1,1):ss.selector.vel(end,2)), [], ss.ndims);
      ss.pos(ss.dynamical_p,:) = reshape(state(ss.selector.pos(1,1):ss.selector.pos(end,2)), [], ss.ndims);
      ss.vel(ss.dynamical_p,:) = dynamicVel;
  end
  if ss.numdyn_m>0;
    ss.m(ss.dynamical_m,:) = state(ss.selector.m(1):ss.selector.m(end));
  end
  if ss.numdyn_k>0;
    ss.k(ss.dynamical_k,:) = state(ss.selector.k(1):ss.selector.k(end));
  end
  if ss.numdyn_r>0;
    ss.r(ss.dynamical_r,:) = state(ss.selector.r(1):ss.selector.r(end));
  end
  if ss.numdyn_c>0;
    ss.c(ss.dynamical_c,:) = state(ss.selector.c(1):ss.selector.c(end));
  end
  if ss.numdyn_u>0;
    ss.u(ss.dynamical_u,:) = state(ss.selector.u(1));
  end
end
