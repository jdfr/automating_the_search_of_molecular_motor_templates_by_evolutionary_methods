%forces a velocity by forcing both the velocity and the acceleration
% updateFlagVars: flag to say that affected flag vars should be set to true
function ss = forceVelocity(ss, updateFlagVars, indexes, startTs, endTs, vel)
  %force velocities to the given values
  ss = forceDerivative(ss, updateFlagVars, 'pos', indexes, startTs, endTs, vel);
  %force accelerations to 0
  ss = forceDerivative(ss, updateFlagVars, 'vel', indexes, startTs, endTs, zeros(size(vel)));
end
