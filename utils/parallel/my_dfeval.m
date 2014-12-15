function [statusOK outputs] = my_dfeval(dfcn, numArgOut, varargin)
%customized dfeval

% %% Check arguments
% error(nargoutchk(1,inf,nargout, 'struct'))
% error(nargchk(2,inf,nargin, 'struct'))% for argOutLp

%%% Create, submit and return

numTasks = numel(varargin{1});

timeout = 60;

jobObj = my_dfevalasync(dfcn,numArgOut,varargin{:});
try
%   waitForState(jobObj,'finished');
  
  goon = true;
  while goon
    waitForState(jobObj,'finished', timeout);
    [pending running finished] = findTask(jobObj);
    goon = numel(finished)<numTasks;
    if goon
      fprintf('%g total tasks, %g done, %g running, %g pending\n', numTasks, numel(finished), numel(running), numel(pending));
    end
  end

  [statusOK outputs] = reapAllResults(jobObj);
catch ME
  fprintf('Destroying job before it is finished!!\n');
  jobObj.destroy;
  rethrow(ME);
end
jobObj.destroy;

%destroy the object, release resources


