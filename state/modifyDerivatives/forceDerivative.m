%forces the value of the derivative of a dynamical variable
% updateFlagVars: flag to say that affected flag vars should be set to true
function ss = forceDerivative(ss, updateFlagVars, varargin) %varargin = {dynvartype, indexes, startTs, endTs, modification}
  ss = modifyDynamics(ss, updateFlagVars, 'enforcements', varargin{:});
end
