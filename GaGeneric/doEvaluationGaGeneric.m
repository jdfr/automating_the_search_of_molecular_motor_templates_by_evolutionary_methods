%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeSpent varargout] = doEvaluationGaGeneric(evparams, taskinfo, multiEvFun, varargin)
timeStart   = cputime;
[varargout{1:nargout-1}] = feval(multiEvFun, varargin{:}, evparams, taskinfo);
timeSpent = cputime-timeStart;
