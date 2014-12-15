function [statusOK outputs] = reapAllResults(job)
%based on getAllOutputArguments, merges errors and output vars alike

% Ensure that only one job has been passed in
if numel(job) > 1
    error('distcomp:job:InvalidArgument', 'The function GetAllOutputArguments requires a single job input');
end
% Get all the tasks from this job
tasks = job.Tasks(:);
%get all error messages
errorMessages = get(job.Tasks, {'ErrorMessage'});
errors = get(job.Tasks, {'Error'});
% Create a cell array to hold the output
outputs   = cell(numel(tasks), 1);
statusOK  = cellfun(@isempty, errorMessages);
outputs(~statusOK) = errors(~statusOK);
oks = find(statusOK);
% Loop over each task an get the output
for i = 1:numel(oks)
    tsk = oks(i);
    outputs{tsk} = tasks(tsk).OutputArguments;
end