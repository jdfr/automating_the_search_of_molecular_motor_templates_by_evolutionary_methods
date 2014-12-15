function [ended states legal] = waitForJobsGeneric(ftst, jobs, decisionFunction, terclus, params) %#ok<INUSD>
%check whether jobs are finished or not. Try it until custom condition is
%met

  fprintf('WAITING FOR JOBS... ');

  if isempty(jobs)
    ended  = false(size(jobs));
    states =  cell(size(jobs));
    legal  =  true(size(jobs));
    return
  end
tempseed = rand(params.randmethod);
  
  legal = cellfun(@(x) isa(x.job, 'distcomp.job') || ...
                       isa(x.job, 'distcomp.simplejob') || ...
                          (isstruct(x.job) && ...
                           ( (isfield(x.job, 'fakeJob') && x.job.fakeJob) || ...
                             (isfield(x.job, 'picassoJob') && x.job.picassoJob) ) ...
                           ), ...
                  jobs);
  if ~all(legal)
    notlegal = jobs(~legal);
    cls      = cellfun(@(x)[class(x) ', '], notlegal, 'uniformoutput', false);
    cls{end} = cls{end}(1:(end-2));
    cls      = horzcat(cls{:});
    mprintf(ftst, 'ERROR!!!! Some jobs are not ''distcomp.job'' but: %s\n', cls);
    st       = dbstack;
    cls      = cell(size(st));
    for k=1:numel(st)
      cls{k} = sprintf('   FILE: %s\n   NAME: %s\n   LINE: %d\n', st(k).file, st(k).name, st(k).line);
    end
    cls      = horzcat(cls{:});
    mprintf(ftst, '   THE STACK FOR THIS WARNING IS:\n%s', cls);
    if numel(jobs(legal))==0
      ended  = true(size(jobs));
      states = repmat({'notjob'}, size(jobs));
      rand(params.randmethod, tempseed);
      return;
    end
  end
  %timeStart = clock;
  legalc =  array2cellMX(legal);% arrayfun(@(x)x, legal, 'uniformoutput', false);
  safeQuery = safeQueryStateFUN(terclus, params);
  while true 
    states  = cellfun(safeQuery, jobs, legalc, 'uniformoutput', false);
%     statest = cell(numel(jobs{1}.job.Tasks), 1); for k=1:numel(statest); statest{k} = jobs{1}.job.Tasks(k).State; end;
%     hayfailed = cellfun(@(x)strcmp(x, 'failed'), statest);
%     hayrunning = cellfun(@(x)strcmp(x, 'running'), statest);
%     hayfinished = cellfun(@(x)strcmp(x, 'finished'), statest);
%     haypending = cellfun(@(x)strcmp(x, 'pending'), statest);
%     hayunavailable = cellfun(@(x)strcmp(x, 'unavailable'), statest);
%     
%     mprintf(ftst, '------------HEY, %d TASKS: TASKS %s FAILED, TASKS %s RUNNING, TASKS %s FINISHED, TASKS %s PENDING, TASKS %s UNAVAILABLE, JOB IS %s!!!!!!\n', numel(statest), mat2str(find(hayfailed)), mat2str(find(hayrunning)), mat2str(find(hayfinished)), mat2str(find(haypending)), mat2str(find(hayunavailable)), jobs{1}.job.State);
    
    ended   = cellfun(@stateEnd, states);
    if decisionFunction(ended) %do this until custom condition is met
      break;
    else
%       pause(1);
      if (terclus.numRanges==1) && (~isstruct(jobs{1}.job))
        lastOK = waitForState(jobs{1}.job, 'finished', terclus.pauseTime); %#ok<NASGU>
      else
        pause(terclus.pauseTime); %DAN: better than pooling the
%         %pause(terclus.pauseTime); %DAN: better than pooling the
%         %cluster for the job state is to leave the cluster to awake you
%         for i=1:size(jobs, 2)
%           waitForState(jobs{i}.job, 'finished', 300); wait
%         end
      end
    end
  end
rand(params.randmethod, tempseed);
  
fprintf('DONE!\n');
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = safeQueryStateFUN(terclus, params)
fun = @(x, legal) safeQueryState(terclus, params, x, legal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = safeQueryState(terclus, params, x, legal)
if legal
  if isfield(x.job, 'picassoJob') && x.job.picassoJob
    s = getPicassoJobState(x.job.jobID, terclus, params);
  else
    s = x.job.State;
  end
else
  s = 'notjob';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ended = stateEnd(state)
  switch state
    case {'finished', 'failed', 'destroyed', 'unavailable', 'notjob'}
      ended = true;
    otherwise
      ended = false;
  end
