%adds springs; optionally returns the indexes to access their parameters
function [ss varargout] = addSprings(ss, varargin) %varargin = {springEnds, k, r, c, dynamical_k, dynamical_r, dynamical_c}) = springVars
  numss = size(varargin{1},1);
  %output indexes
  if nargout>1; varargout{1} = size(ss.springEnds,1)+(1:numss)'; end %indexes
  ss = auxAssignment(ss, numss, 1, varargin, ss.springDefaultVals, ss.springVars, 'spring');
