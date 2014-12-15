%add points; optionally returns the indexes to access their parameters
function [ss varargout] = addPoints(ss, varargin) %varargin = {pos, vel, m, dynamical_p, dynamical_m} = pointVars
  numps = size(varargin{1},1);
  %output indexes
  if nargout>1; varargout{1} = size(ss.pos,1)+(1:numps)'; end %indexes
  ndims = size(varargin{1},2);
  if ndims ~= size(ss.pos,2); error('incorrect number of dimensions'); end
  ss = auxAssignment(ss, numps, 1, varargin, ss.pointDefaultVals, ss.pointVars, 'point');
