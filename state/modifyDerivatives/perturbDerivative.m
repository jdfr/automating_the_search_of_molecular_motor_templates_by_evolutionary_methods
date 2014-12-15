%adds a perturbation to the derivative of a dynamical variable
% updateFlagVars: flag to say that affected flag vars should be set to true
function ss = perturbDerivative(ss, updateFlagVars, varargin) %varargin = {dynvartype, indexes, startTs, endTs, modification}
  ss = modifyDynamics(ss, updateFlagVars, 'perturbations', varargin{:});
end
