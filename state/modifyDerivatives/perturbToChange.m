%adds a perturbation to the derivative of a dynamical variable, to drive it
%to a new value in a given amount of time
% updateFlagVars: flag to say that affected flag vars should be set to true
function ss = perturbToChange(ss, updateFlagVars, dynvartype, indexes, startTs, deltaTs, newVals)
  initVals = ss.(dynvartype)(indexes,:);
  ratio = (newVals-initVals)./deltaTs;
  ss = perturbDerivative(ss, updateFlagVars, dynvartype, indexes, startTs, startTs+deltaTs, ratio);
end
