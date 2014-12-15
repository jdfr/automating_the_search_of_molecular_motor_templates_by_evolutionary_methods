function subnomdir = GAGeneric(params)
startTime = clock;
initcputime = cputime;

%default values for several parameters
%params = putDefaultParameters(params);

terclus = params.terclus;
params  = rmfield(params, 'terclus');

params.evparams.logInEval = params.logInEval;

params.basedir = params.nomdir;

if filesep=='/' %make sure that file separators are OK
  params.nomdir = strrep(params.nomdir, '\', filesep);
elseif filesep=='\'
  params.nomdir = strrep(params.nomdir, '/', filesep);
end
params.nomdir = [params.nomdir,filesep,num2str(params.num_iteracion),'_',num2str(params.presion),'_',num2str(params.ensayo)];
mkdir(params.nomdir);
subnomdir = params.nomdir;
params.nomdir = [params.nomdir,filesep];
params.absdir = [pwd filesep params.nomdir];

if ~isfield(params, 'getRasters')
  params.getRasters = true;
end
if ~isfield(params, 'rewindNans')
  params.rewindNans = false;
end
if ~isfield(params, 'nansReplaced')
  params.nansReplaced = false;
end
if ~isfield(params, 'newsQuotaMode')
  params.newsQuotaMode    = 'none';
end
if (~isfield(params, 'selection')) || (~isfield(params.selection, 'mode'))
  params.selection.mode = 'PACO';
end

if (terclus.numRanges==1) && (~ischar(terclus.jobMgr))
  terclus.pauseTime = Inf;
end
if ~isfield(params, 'saveRnd')
  params.saveRnd = true;
end


params.fpob = fopen([params.nomdir,'poblacion.txt'],'a');
params.farb = fopen([params.nomdir,'arbol.txt'],'a');
params.fres = fopen([params.nomdir,'resultado.txt'],'a');
params.ftst = fopen([params.nomdir,'timestats.txt'],'a');
params.allfs = {params.fpob; params.farb; params.fres; params.ftst};

usePicasso = ischar(terclus.jobMgr) && strcmpi(terclus.jobMgr, 'picasso');
if usePicasso && terclus.remoteSSH
  terclus.scriptBaseName = terclus.scriptFileName;
  terclus.scriptFileName = [terclus.scriptFileName params.ahora '.sh'];
  terclus.ssh.nomdir = strrep(params.nomdir, '\', '/');
  channel = createSSHChannel(params, terclus);

  dirs=cellfun(@(x)x{1}, regexp(terclus.ssh.nomdir, '([^/]+)', 'tokens'), 'uniformoutput', false);
  mydir = '';
  for k=1:numel(dirs)
    [channel res] = executeCommand(channel, ['cd ' mydir '; ls ' dirs{k}]);
    if isempty(res) || all(isspace(res))
      channel = executeCommand(channel, ['cd ' mydir '; mkdir ' dirs{k}]);
    end
    mydir = [mydir dirs{k} '/']; %#ok<AGROW>
  end
  
  channel.close;
  
  clear channel;
end

terclus.elapsed = 0;

params.writeMejor            = params.writeMejorFun(params.fres, params.genome2str);
params.writeIndividual       = params.writeIndividualFun(params.fpob, params.genome2str);
params.writeGenealogy        = params.writeGenealogyFun(params.farb);

if params.plotting && (~isfield(params, 'writeGeneration'))
  params.writeGeneration     = @writePNGs;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fun = writeMejorFun(fres, switchGenomeStr)
%   fun = @(mejorEnConserva) ...
%     fprintf(fres,'%g %g %g %s %s %s\n',mejorEnConserva.gen, mejorEnConserva.rind, mejorEnConserva.num, num2hex(mejorEnConserva.data.rndSeed), mat2str(mejorEnConserva.data.fitness), switchGenomeStr(mejorEnConserva.data.genome{1})); % actualizamos el fichero de resultados
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fun = writeIndividualFun(fpob, switchGenomeStr)
%   fun = @(generation, rangeind, index, pop) ...
%     fprintf(fpob,'%d %d %d %s %s %s %-06.3f %-06.3f %-06.3f %s\n',generation,rangeind,index,num2hex(pop.rndSeed(index)),mat2str(pop.fitness(index)), mat2str(pop.displacement(index)), pop.devTimeSpent(index), pop.timeSpent(index), pop.cpuTimeSpent(index), switchGenomeStr(pop.genome{index}));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fun = writeGenealogyFun(farb)
%   fun = @(generationNew, rangeindNew, indexNew, changeStr, generationOld, rangeindOld, indexOld) ...
%     fprintf(farb,'%d %d %d %s %d %d %d\n',generationNew,rangeindNew,indexNew,changeStr,generationOld, rangeindOld, indexOld);

%this struct holds a destacated first field: genome, that identifies the
%individual, and other fields, that are the individual's attributes
params.indivStruct(2,:) = cellfun(@(x){x}, params.indivStruct(2,:), 'uniformoutput', false);
params.indivStruct = struct(params.indivStruct{:});
if ~isfield(params, 'fieldsToFill')
  params.fieldsToFill = 2:numel(fieldnames(params.indivStruct)); %skip first field 'genome' when reaping results
end
%it is convenient to pre-fabricate this structure
jobStruct = {'G' 'indexNewG' 'job' 'rangeid' 'generation'};
% legend: G==individuals of this range
%         indexNewG== indexes in G of truly new individuals (only they will be evaluated)
%         job== job handle in charge of evaluating this range
%         rangeid== the id if this range
%         generation== the generation
jobStruct = reshape([jobStruct; cell(size(jobStruct))], 1, []);
params.jobStruct = struct(jobStruct{:}); clear jobStruct;

if ~isempty(params.makeupParams)
  params = params.makeupParams(params);
end

%here, we will store the best individual of the simulation. It cannot be
%safely stored in the population, since it might be overwritten
mejorEnConserva = struct('gen', [], 'rind', [], 'num', [], 'data',  {params.indivStruct});
jobsActive = cell(0,1);
try
    maxActive  = terclus.concJobs; %#ok<NODEF>
    terclus.rangeSizes = diff(terclus.ranges,1,2)+1;

    P = repmatStruct(params.indivStruct, [params.tam_poblacion,1]);
    
    [P params terclus] = prepareInitialPopulation(P, params, terclus);
    
    if params.rewindNans
      %if we are rewinding individuals with nan fitness, initialize this
      %field to all false. This way, when we reach evolucionGeneric, no
      %individual will be reverted, since the initial population, having no
      %ancestors, cannot be reverted
      params.Preverted = false(size(P.fitness));
    end
    
    i=1;
    for na=1:numel(terclus.rangeSizes)
        for nb=1:terclus.rangeSizes(na)
            params.writeGenealogy(params.initialGeneration, na, i, '=', 0, 0, 0);
            params.writeIndividual(params.initialGeneration, na, i, P);
            i=i+1;
        end
    end
    %for each range, the generation it was last updated
    generationInRange = repmat(params.initialGeneration, (size(terclus.rangeSizes)));
    %find the best in the population
    mejorEnConserva = findMejor(P, mejorEnConserva, terclus.ranges, generationInRange);
    fprintf('TIME ELAPSED: %g\n', etime(clock, startTime));

    doTermination = false;
    lastTime = clock;
    for cont=(params.initialGeneration+1):(params.initialGeneration+params.generaciones) %for each generation
        for rind=1:size(terclus.ranges,1) %for each range

            %if max amount of concurrent jobs is reached, wait for some to end,
            %and store it (or them, if several end at once)
            if numel(jobsActive)>=maxActive
                decisionFun = @any; %wait for some job to end
            else
                decisionFun = @alwaysTrue; %reap jobs if they are finished, but do not wait for them
            end
            [ended states] = waitForJobsGeneric(params.ftst, jobsActive, decisionFun, terclus, params);
            %if any job ended: reaping time!
            tooSlowEvaluation = etime(clock, lastTime)>params.tooLongEvTime; %#ok<NASGU>
            lastTime = clock;
            if any(ended)
                endedI = find(ended);
                %reap them
                [jobsActive terclus] = reapJobs(params, jobsActive, endedI, states(ended), terclus);
                reapedids = cellfun(@(x)sprintf('G=%gR=%g, ', x.generation, x.rangeid), jobsActive(ended), 'uniformoutput', false);
                reapedids = horzcat(reapedids{:});
                %store results in population
                [jobsActive, P, generationInRange, mejorEnConserva, params] = updateGA(params, terclus, endedI, jobsActive, P, generationInRange, mejorEnConserva, startTime);

%                 stopIfTooTall     = (~isempty(params.stopIfTooTall)) && any((P.dimensions(:,1)-P.offset(:,1))>params.stopIfTooTall);
%                 if stopIfTooTall || tooSlowEvaluation
%                   mprintf(params.ftst, 'We stopped prematurely!!! stopIfTooTall=%s, tooSlowEvaluation=%s\n', any2str(stopIfTooTall), any2str(tooSlowEvaluation));
%                   doTermination=true;
%                   break;
%                 end
                
                mprintf(params.ftst, 'REAPED  JOBS FOR %s\n   TIME ELAPSED: %f (CPU: %f), TIME IN CLUSTER: %f\n', reapedids, etime(clock, startTime), cputime-initcputime, terclus.elapsed);
                % MANDAR AL CLUSTER PARA PINTAR GRAFOS PIXEL
            end
            %launch a new job!!!!
            nj = params.jobStruct;
            [nj.G nj.indexNewG terclus P parentdata]  = evolucionGeneric(params, P,cont,rind,generationInRange,terclus,mejorEnConserva); % proxima generacion
            nj.generation = cont;
            nj.rangeid    = rind;
            nj.job        = evaluarRemote(terclus, nj.G.genome(nj.indexNewG), params, nj.generation, nj.rangeid, nj.indexNewG);
            nj.parentdata = parentdata; clear parentdata;
            longs         = cellfun('prodofsize', nj.G.genome(nj.indexNewG));
            jobsActive{end+1} = nj; %#ok<AGROW>
            mprintf(params.ftst, 'SPAWNED GEN=%g, RANGE=%g, Nº OF STRINGS: %g\n   TIME ELAPSED=%f (CPU: %f), TIME IN CLUSTER=%f\n', cont, rind, numel(nj.indexNewG), etime(clock, startTime), cputime-initcputime, terclus.elapsed);
            mprintf(params.ftst, '   NEW STRING LENGTH   STATS: MIN=%g, MAX=%g, MEAN=%g, MEDIAN=%g\n', min(longs), max(longs), mean(longs), median(longs));
        end
        if doTermination
          break;
        end
    end

    if ~doTermination
        mprintf(params.ftst, 'BEFORE FINAL REAPING:\n   TIME ELAPSED=%f (CPU: %f), TIME IN CLUSTER=%f\n', etime(clock, startTime), cputime-initcputime, terclus.elapsed);
        %wait for all extant jobs
        [nevermind states] = waitForJobsGeneric(params.ftst, jobsActive, @all, terclus, params); %#ok<ASGLU>
        %reap all extant jobs
        [jobsActive terclus] = reapJobs(params, jobsActive, 1:numel(jobsActive), states, terclus);
        reapedids = cellfun(@(x)sprintf('G=%gR=%g, ', x.generation, x.rangeid), jobsActive, 'uniformoutput', false);
        reapedids = horzcat(reapedids{:});
        mprintf(params.ftst, 'AT LAST, REAPED JOBS FOR %s\n   TIME ELAPSED: %g (CPU: %f), TIME IN CLUSTER: %g\n', reapedids, etime(clock, startTime), cputime-initcputime, terclus.elapsed);
        %store results in population
        [jobsActive, P, generationInRange, mejorEnConserva, params] = updateGA(params, terclus, 1:numel(jobsActive), jobsActive, P, generationInRange, mejorEnConserva, startTime); %#ok<ASGLU>
        params.writeMejor(mejorEnConserva); %write best individual for this generation
    end

catch ME
    cellfun(@fclose, params.allfs);
    fprintf('closed all files...\n');
    %  destroy(findJob(terclus.jobMgr, 'Tag', terclus.tag));
    jobsActive          = carefulDestroy(jobsActive); %#ok<NASGU>
    fprintf('destroyed all jobs...\n');
    rethrow(ME);
end
cellfun(@fclose, params.allfs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initializes population
function   [P terclus params] = startPopulation(params, terclus, P, Generation)
%get fitnesses for all

ranges          = terclus.ranges;
cj              = terclus.concJobs;
totalj          = size(ranges,1);
numconc         = 0;
jobs            = cell(1,totalj);
activeJobs      = zeros(0,1);
lP              = 0;
while (lP<totalj)
    %launch jobs
    while (lP<totalj) && (numconc<=cj)
        lP            = lP+1;
        numconc       = numconc+1;
        nj            = params.jobStruct;
        nj.indexNewG  = 1:(diff(ranges(lP,:))+1);%ranges(lP,1):ranges(lP,2);
        nj.G          = getSubStruct(P, ranges(lP,:), 'intensive');
        
        nj.rangeid    = lP;
        nj.generation = Generation;
        nj.job        = evaluarRemote(terclus, nj.G.genome, params, nj.generation, nj.rangeid, nj.indexNewG);
        jobs{lP}      = nj;
        activeJobs = [activeJobs; lP]; %#ok<AGROW>
    end

    %wait for any job to finish
    [ended states]          = waitForJobsGeneric(params.ftst, jobs(activeJobs), @any, terclus, params);
    %reap finished jobs
    [jobs terclus] = reapJobs(params, jobs, activeJobs(ended), states(ended), terclus);

    endedI = find(ended);
    for k=1:numel(endedI) %dispose non-needed fields
        jobs{activeJobs(endedI(k))}.indexNewG = [];
        jobs{activeJobs(endedI(k))}.job       = [];
    end
    numconc                 = numconc - sum(ended);
    activeJobs              = activeJobs(~ended);
    if size(activeJobs,1)>size(activeJobs,2) %this may happen sometimes
        activeJobs            = activeJobs'; %just to be safe
    end
end

%wait for all extant jobs to finish
[nevermind states]      = waitForJobsGeneric(params.ftst, jobs(activeJobs), @all, terclus, params); %#ok<ASGLU>
%reap finished jobs
[jobs terclus] = reapJobs(params, jobs, activeJobs, states, terclus);
for k=1:numel(activeJobs) %dispose non-needed fields
    jobs{activeJobs(k)}.indexNewG = [];
    jobs{activeJobs(k)}.job       = [];
end

for ei=1:numel(jobs)
  areError = find(cellfun(@(x)not(isempty(x)) && isfield(x, 'isError'), jobs{ei}.G.raster));
  if not(isempty(areError))
    errs = jobs{ei}.G.raster(areError); %#ok<NASGU>
    idxs = find(areError); %#ok<NASGU>
    save([params.nomdir sprintf('errors_%04d_%d.mat', 0, 1)], 'errs', 'idxs');
    error('There has been an error while evaluating at least one individual!!!!');
  end
end

%update population
for k=1:numel(jobs)
  P = assignSubStruct(P, ranges(k,:), 'intensive', jobs{k}.G);
  if(params.plotting)
    params.writeGeneration(params, jobs{k}.G, Generation, jobs{k}.rangeid, 1:numel(jobs{k}.G.genome), ranges(k,1));
  end
end

clear jobs;

if params.plotting
  %P.raster(1:end) = {[]};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jobstructs terclus] = reapJobs(params, jobstructs, toReap, states, terclus)
tempseed = rand(params.randmethod);
usePicasso = ischar(terclus.jobMgr) && strcmpi(terclus.jobMgr, 'picasso');
for z=1:numel(toReap)
    e = toReap(z);
    if usePicasso
      thisjob = jobstructs{e}.job;
      if terclus.remoteSSH
        channel = createSSHChannel(params, terclus);

        sp = channel.createSCPClient;
        qw=javaArray('java.lang.String', numel(thisjob.outputFiles));
        for k=1:numel(thisjob.outputFiles)
          qw(k) = java.lang.String(thisjob.outputFiles{k});
        end
        sp.get(qw, params.nomdir);

        jobID = thisjob.jobID;
        jobID = jobID((jobID>='0') & (jobID<='9'));
        channel  =  executeCommand(channel, ['rm ' terclus.scriptBaseName '*' jobID '* ; cd ' terclus.ssh.nomdir '; ' terclus.ssh.rmCommand]);
        %remove stdout and stderr, TaskInputs and TaskOutputs of this generation  
        channel.close;
      end
      %jd = {'FinishTime', 'FileDependencies', 'PathDependencies', 'JobData', 'UserData'};
      thisjob.Configuration = '';
      thisjob.Name = mat2str(thisjob.jobID);
      thisjob.ID = thisjob.jobID;
      thisjob.UserName = '';
      thisjob.Tag = '';
      thisjob.State = states{z};
      thisjob.CreateTime = '';
      thisjob.SubmitTime = '';
      thisjob.StartTime = '';
      thisjob.FinishTime = '';
      thisjob.FileDependencies = '';
      thisjob.PathDependencies = '';
      thisjob.JobData = {'picassoJob', thisjob.nomdir, thisjob.inputFileNameTemplate, thisjob.outputFileNameTemplate, thisjob.jobIDFileName};
      td = {'Configuration', 'Name', 'ID', 'Function', 'NumberOfOutputArguments', 'InputArguments', 'OutputArguments', 'CaptureCommandWindowOutput', 'CommandWindowOutput', 'State', 'Error', 'ErrorMessage', 'ErrorIdentifier', 'CreateTime', 'StartTime', 'FinishTime', 'UserData'};
      td = [td; repmat({cell(thisjob.numTasks, 1)}, size(td))]; %#ok<AGROW>
      thisjob.Tasks = struct(td{:});
      for k = 1:thisjob.numTasks
        inputFile = [params.nomdir sprintf(terclus.inputFileNameTemplate, k)];
        load(inputFile);
        thisjob.Tasks(k).ID = task.numTask;
        thisjob.Tasks(k).State = 'finished';
        outputFile = [params.nomdir sprintf(terclus.outputFileNameTemplate, k)];
        if ~exist(outputFile, 'file')
          error('The file <<%s>> is expected, but not found!!!!\n', outputFile);
        end
        load(outputFile);
        if exist('output', 'var')
          thisjob.Tasks(k).OutputArguments = output;
        elseif exist('err', 'var')
          thisjob.Tasks(k).Function = task.Function;
          thisjob.Tasks(k).NumberOfOutputArguments = task.NumberOfOutputArguments;
          thisjob.Tasks(k).InputArguments = task.InputArguments;
          thisjob.Tasks(k).generation = task.generation;
          thisjob.Tasks(k).ErrorMessage = 'There is error!!!';
          thisjob.Tasks(k).Error = err;
          thisjob.Tasks(k).UserData = task;
        else
          thisjob.Tasks(k).ErrorMessage = ['UNKNOWN PROBLEM WITH ' outputFile];
          thisjob.Tasks(k).Error = thisjob.Tasks(k).ErrorMessage;
        end
        delete(inputFile);
        delete(outputFile);
        clear task output err;
      end
      jobstructs{e}.job = thisjob; 
    end
    jrid = jobstructs{e}.rangeid;
    switch states{z}
        case 'failed'
            mprintf(params.ftst, 'JOB %s FOR RANGE %g = [%g ... %g] !!!!!\n', upper(states{e}), jrid, terclus.ranges(jrid,1), terclus.ranges(jrid,2));
            makeFailedMAT(params, jobstructs{e}, 'FailedJob');
            error('This should never happen (job %s)!!!', states{z});
            %terclus = addElapsedTime(terclus, nan);
        case {'destroyed', 'unavailable', 'notjob'}
            mprintf(params.ftst, 'JOB %s FOR RANGE %g = [%g ... %g] !!!!!\n', upper(states{e}), jrid, terclus.ranges(jrid,1), terclus.ranges(jrid,2));
            error('This should never happen (job %s)!!!', states{z});
            %terclus = addElapsedTime(terclus, nan);
        case 'finished'
            job         = jobstructs{e}.job;
            limitsByW   = job.UserData{1};
            ntasks      = numel(job.Tasks);
            indexes     = jobstructs{e}.indexNewG;
            for k=1:ntasks
                if isempty(job.Tasks(k).ErrorMessage)
                    goon = true;
                    output = job.Tasks(k).OutputArguments;
                    %there seems to be a race condition: sometimes output can be
                    %empty, supposedly because the job has not yet transferred the
                    %output. The workaround is simple, but it's stupid to have to do
                    %it
                    if isempty(output)
                        hayprob = true;
                        if isinf(terclus.pauseTime)
                          pause(1);
                        else
                          pause(terclus.pauseTime);
                        end
                        output = job.Tasks(k).OutputArguments;
                        if isempty(output)
                            mprintf(params.ftst, 'ERROR!!! THERE SEEMS TO BE A CONDITION PREVENTING TO REAP GENERATION %g RANGE %g TASK %g!!!\nERROR: %s\nOUTPUT: %s\n', jobstructs{e}.generation, jobstructs{e}.rangeid, k, showError(job.Tasks(k).Error), any2str(job.Tasks(k).OutputArguments));
                            makeFailedMAT(params, jobstructs{e}, 'THERE SEEMS TO BE A CONDITION PREVENTING TO REAP A RANGE');
                            error('This should never happen (job %s)!!!', states{z});
                            goon = false;
                        end
                    else
                        hayprob = false;
                    end
                    if goon
                        if hayprob
                            mprintf(params.ftst, 'WE HAVE BEEN ABLE TO RECOVER FROM THE CONDITION :-)\n');
                        end
                        %fill results
                        fnames = fieldnames(jobstructs{e}.G);
                        fieldsToFill = params.fieldsToFill;
                        for zz=(1):numel(fieldsToFill) 
                          %we shift zz by 1 (output{zz+1}) since the first
                          %value is the global time spent
                          jobstructs{e}.G.(fnames{fieldsToFill(zz)})(indexes(limitsByW(k,1):limitsByW(k,2)),:) = output{zz+1};
                        end
                        terclus = addElapsedTime(terclus, output{1});
                    else
                        if hayprob
                            mprintf(params.ftst, 'WE HAVEN''T BEEN ABLE TO RECOVER FROM THE CONDITION :-(\n');
                            makeFailedMAT(params, jobstructs{e}, 'WE HAVEN''T BEEN ABLE TO RECOVER FROM THE CONDITION :-(');
                            error('This should never happen (job %s)!!!', states{z});
                        end
                    end
                else
                    err = job.Tasks(k).Error;
                    mprintf(params.ftst, 'JOB %g HAS TASK %g FINISHED WITH ERROR!!! %s', jrid, k, showError(err));
                    makeFailedMAT(params, jobstructs{e}, 'TASK WITH UNTREATED ERROR');
                    error('This should never happen (job %s)!!!', states{z});
                end
            end

    end
    if (~strcmp('notjob', states{z})) && (isa(jobstructs{e}.job, 'distcomp.job') || isa(jobstructs{e}.job, 'distcomp.simplejob'))
        destroy(jobstructs{e}.job);
    end
end
rand(params.randmethod, tempseed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeFailedMAT(params, jobstruct, condition)
failed = struct;
failed.job = struct;
jd = {'Configuration', 'Name', 'ID', 'UserName', 'Tag', 'State', 'CreateTime', 'SubmitTime', 'StartTime', 'FinishTime', 'FileDependencies', 'PathDependencies', 'JobData', 'UserData'};
td = {'Configuration', 'Name', 'ID', 'Function', 'NumberOfOutputArguments', 'InputArguments', 'OutputArguments', 'CaptureCommandWindowOutput', 'CommandWindowOutput', 'State', 'Error', 'ErrorMessage', 'ErrorIdentifier', 'CreateTime', 'StartTime', 'FinishTime', 'UserData'};
for k=1:numel(jd);
  failed.job.(jd{k}) = jobstruct.job.(jd{k});
end
failedjob.Tasks = struct;
for k=1:numel(jobstruct.job.Tasks)
  for kk=1:numel(td);
    failed.job.Tasks(k).(td{kk}) = jobstruct.job.Tasks(k).(td{kk});
  end
end
failed.G          = jobstruct.G;
failed.indexNewG  = jobstruct.indexNewG;
failed.rangeid    = jobstruct.rangeid;
failed.generation = jobstruct.generation;
failed.condition  = condition;
save([params.nomdir sprintf('failedjob_%d_%d.mat', jobstruct.generation, jobstruct.rangeid)], 'failed');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function terclus = addElapsedTime(terclus, elapsed)
%if ~any(isnan(elapsed))
terclus.elapsed = terclus.elapsed+sum(elapsed);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mejorEnConserva = findMejor(P, mejorEnConserva, ranges, generationInRange)
[maxi, mejor] = max(P.fitness);               %#ok<ASGLU> % mejor individuo
% mejor = noplanos(mejor);
if (~isempty(mejor)) %&& ( isempty(mejorEnConserva.data.fitness) || (maxi>=mejorEnConserva.data.fitness) )
    mejorEnConserva.rind  = findRange(ranges, mejor);
    mejorEnConserva.gen   = generationInRange(mejorEnConserva.rind);
    mejorEnConserva.num   = mejor;
    mejorEnConserva.data  = getSubStruct(P, mejor, 'extensive');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = findRange(ranges, idx)
solvfun = @(ix) find((ix>=ranges(:,1)) & (ix<=ranges(:,2)));
if numel(idx)==1
    r = solvfun(idx);
else
    r = arrayfun(solvfun, idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%update variables after a job (some jobs) is (are) reaped
function [jobsActive, P, generationInRange, mejorEnConserva, params] = updateGA(params, terclus, toUpdate, jobsActive, P, generationInRange, mejorEnConserva, startTime) %#ok<INUSL>
ranges = terclus.ranges;
for k=1:numel(toUpdate)
    ei    = toUpdate(k);
    rpri  = jobsActive{ei}.rangeid;
    gen   = jobsActive{ei}.generation;
    if params.rewindNans
      %if we are rewinding individuals with nan fitness, localize them
      areNan           = isnan(jobsActive{ei}.G.fitness);
      if isfield(jobsActive{ei}.parentdata, 'notNew')
        %if we also are using a quota, make sure that we don't erroneously
        %try to rewind new individuals, since they have no ancestor.
        areNan         = areNan & jobsActive{ei}.parentdata.notNew;
      end
      %get the indexes of the individuals to be rewound
      rang             = (ranges(rpri,1):ranges(rpri,2))';
      rang             = rang(areNan);
      params.Preverted = rang;
      %get the individuals themselves.
      params.Pold      = getSubStruct(P, jobsActive{ei}.parentdata.idxInP(rang), 'extensive');
    end
    P   = assignSubStruct(P, ranges(rpri,:), 'intensive', jobsActive{ei}.G);
    generationInRange(rpri) = gen;
    
    areError = find(cellfun(@(x)not(isempty(x)) && isfield(x, 'isError'), jobsActive{ei}.G.raster(jobsActive{ei}.indexNewG)));
    if not(isempty(areError))
      errs = jobsActive{ei}.G.raster(jobsActive{ei}.indexNewG(areError)); %#ok<NASGU>
      idxs = jobsActive{ei}.indexNewG(areError) + ranges(rpri,1) - 1; %#ok<NASGU>
      save([params.nomdir sprintf('errors_%04d_%d.mat', gen, rpri)], 'errs', 'idxs');
      error('There has been an error while evaluating at least one individual!!!!');
    end
    
    if(params.plotting)
      params.writeGeneration(params, jobsActive{ei}.G, gen, rpri, jobsActive{ei}.indexNewG, ranges(rpri,1));
      %P.raster(ranges(rpri,:)) = {[]};
    end
        
    %recalculate best individual
    mejorEnConserva = findMejor(P, mejorEnConserva, ranges, generationInRange); %#ok<ASGLU>
    %output stage
    if rpri==1
        params.writeMejor(mejorEnConserva); %write best individual for this generation
    end
    for z=ranges(rpri,1):ranges(rpri,2)
      params.writeIndividual(gen, rpri, z, P);
    end
    fprintf(params.ftst, '   POPULATION GENERATIONS UPDATE: GENERATION=%g, RANGE=%g, ALL=%s    (TIME ELAPSED: %s)\n', gen, rpri, mat2str(generationInRange), mat2str(etime(clock, startTime)));
    
    if (k==numel(toUpdate)) && params.saveRnd
      frnd = fopen([params.nomdir 'rndseed.txt'], 'w');
      fprintf(frnd, 'Generation %04d, range %d, rndseed for ''%s'': uint32(', gen, rpri, params.randmethod);
      fprintf(frnd, mat2str(rand(params.randmethod)));
      fprintf(frnd, ')\n');
      fclose(frnd);
    end
end
%wipe out reaped jobs
jobsActive(toUpdate) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%destroy distributed jobs carefully
function jobsActive = carefulDestroy(jobsActive)
if (~isempty(jobsActive)) && iscell(jobsActive)
    marked = false(size(jobsActive));
    for k=1:numel(jobsActive)
        if isstruct(jobsActive{k}) && isfield(jobsActive{k}, 'job') && ( isa(jobsActive{k}.job, 'distcomp.job') || isa(jobsActive{k}.job, 'distcomp.simplejob')) && (~isempty(jobsActive{k}.job))
            for m=1:numel(jobsActive{k}.job)
                %if ~strcmp(jobsActive{k}.job(m).State, 'destroyed')
                destroy(jobsActive{k}.job(m));
                %end
            end
            marked(k) = true;
        end
    end
    jobsActive(marked) = [];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function str = showError(error)
% if isempty(error)
%     str='<EMPTY ERROR>';
% elseif (~isstruct(error)) || (~all(ismember(lower(fieldnames(error)), {'message', 'identifier', 'stack'})))
%     str=['<NO STRUCT ERROR: ' any2str(error) '>'];
% else
%     lines = cell(size(error.stack));
%     for m=1:numel(error.stack)
%         lines{m} = sprintf('    File: %s\n  Name: %s\n  Line: %g\n', error.stack(m).file, error.stack(m).name, error.stack(m).line);
%     end
%     str = horzcat(...
%         sprintf('ERROR:\n  Message:    <<%s>>\n  Identifier: <<%s>>\n  Stack:\n', error.message, error.identifier), ...
%         lines{:});
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newjob = evaluarRemote(terclus, G, params, generation, rangeid, indexesG)%, drawing)

taskinfo = struct('ngen', generation, 'rid', rangeid, 'ntask', [], 'dir', [], 'g2s', []);
if params.logInEval
  taskinfo.dir = params.absdir;
  taskinfo.g2s = params.genome2str;
end

rndSeeds = rand(size(G))*1e5; %we must multiply by 1e5 since numbers in the range 0-1 will not give rise to diverse random sequences

tempseed = rand(params.randmethod);

nargs = numel(params.fieldsToFill)+1;
if isempty(terclus.jobMgr)
  outputs = cell(1,nargs);
  limitsByW = [1 numel(G)];
  returnIndividuals = repmat(params.plotting && params.getRasters, size(G));
  taskinfo.task = 1;
  [outputs{:}] = doEvaluationGaGeneric(params.evparams,taskinfo,params.multiEvFun, G, indexesG, rndSeeds, returnIndividuals);
  newjob = struct('Tasks', {struct('OutputArguments', {outputs}, 'ErrorMessage', {''})}, ...
                  'State', {'finished'}, ...
                  'UserData', {{limitsByW, returnIndividuals}}, ...
                  'fakeJob', {true});
else
  newjob = sendRemoteJob(terclus, params, G, indexesG, rndSeeds, taskinfo, nargs, generation);
end
rand(params.randmethod, tempseed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newjob = sendRemoteJob(terclus, params, G, indexesG, rndSeeds, taskinfo, nargs, generation)
  if isfield(terclus, 'pauseTimeAfterError')
    pt = terclus.pauseTimeAfterError;
  else
    pt = 10;
  end
  if isfield(terclus, 'submitTimes')
    st = terclus.submitTimes;
  else
    st = 3;
  end
  
  usePicasso = ischar(terclus.jobMgr) && strcmpi(terclus.jobMgr, 'picasso');
  
  for z=1:st
    newjob = [];
    try
      concW       = terclus.concurrentWorkers; %number of concurrent workers

      [limitsByW indsByW] = calculateRanges(numel(G), concW);

      if usePicasso
        newjob = struct('inputFiles', {cell(1, numel(indsByW))}, 'outputFiles', {cell(1, numel(indsByW))});
      else
        newjob  = createJob(terclus.jobMgr,terclus.jobArgs{:});
      end
      for k=1:numel(indsByW)
        returnIndividuals = repmat(params.plotting && params.getRasters, size(G));
        taskinfo.ntask = k;
        taskargs = {@doEvaluationGaGeneric, nargs, {params.evparams, taskinfo,params.multiEvFun, G(limitsByW(k,1):limitsByW(k,2)), indexesG(limitsByW(k,1):limitsByW(k,2)), rndSeeds(limitsByW(k,1):limitsByW(k,2)), returnIndividuals(limitsByW(k,1):limitsByW(k,2))}};
        if usePicasso
          %evparams will be accessed indirectly, to avoid multiple copies of a fat structure
          taskargs{3}{1} = terclus.evparamsFileNameTemplate;
          newjob = createPicassoTask(params, terclus, k, generation, newjob, taskargs{:});
        else
          createTask(newjob,taskargs{:});
        end
      end
      newjob.UserData = {limitsByW, returnIndividuals, params.ftst};
      if usePicasso
        newjob = submitPicassoJob(params, terclus, newjob, numel(indsByW), generation);
      else
        %set(get(newjob, 'Tasks'), 'CaptureCommandWindowOutput', true);
        submit(newjob);
      end
      break;
    catch ME
      mprintf(params.ftst, 'Try %d. There has been an error while trying to submit a job:\n%s\n\nTrying to destroy the spoiled job... ', z, showError(ME));
      if (~isempty(newjob)) && (~usePicasso)
        destroy(newjob);
      end
      mprintf(params.ftst, 'Destroyed!. ');
      if z<st
        mprintf(params.ftst, 'Waiting %s seconds... ', mat2str(pt));
        pause(pt);
        mprintf(params.ftst, 'Now, trying again...\n');
      else
        mprintf(params.ftst, 'After %d tries spaced %s seconds between them, we have not been able to spawn a job. Simulation will abort...\n', z, mat2str(pt));
        rethrow(ME);
      end
    end
  end
  if z>1
    mprintf(params.ftst, 'After %d tries, we have spawned a new job!!!\n', z);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writePNGs(params, G, generation, rangeid, indexNewG, rangeStartIdx)
indexes = rangeStartIdx-1+indexNewG;

res = struct('Generation', {generation}, 'rangeid', {rangeid}, 'indexes', {indexes}, 'indexNewG', {indexNewG}, 'rasters', {G.raster(indexNewG)}); %#ok<NASGU>
save([params.nomdir sprintf('lastFrames%04d_%d.mat', generation, rangeid)], 'res');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P params terclus] = prepareInitialPopulation(P, params, terclus)
if isfield(params, 'initialPopulation')
  if isstruct(params.initialPopulation)
    inds = params.initialPopulation.individual;
    if ~isempty(setxor(inds, 1:numel(inds)))
      error('Some indexes of individuals are not in the range 1..numel(indexes)!!!!');
    end
    fnames = fieldnames(params.initialPopulation);
    for k=1:params.initialPopulation.numFieldsPop
      P.(fnames{k})(inds) = params.initialPopulation.(fnames{k})(1:end);
    end
    sameranges = terclus.numRanges==params.initialPopulation.oldparams.terclus.numRanges;
    if ~sameranges
      rangesExpanded          = zeros(size(P.genome));
      for k=1:size(terclus.ranges,1)
        rangesExpanded(terclus.ranges(k,1):terclus.ranges(k,2)) = k;
      end
      newres = struct('Generation', {repmat(params.initialGeneration, size(P.genome))}, 'rangeid', {rangesExpanded}, 'indexes', {(1:numel(P.genome))'}, 'indexNewG', {(1:numel(P.genome))'}, 'rasters', {[]});
      clear rangesExpanded;
    end
    for k=1:numel(params.initialPopulation.oldparams.terclus.numRanges)
      fname = [params.initialPopulation.basedir filesep sprintf('lastFrames%04d_%d.mat', params.initialPopulation.generation(1), k)];
      if exist(fname, 'file')
        if sameranges
          copyfile(fname, [params.nomdir sprintf('lastFrames%04d_%d.mat', params.initialGeneration, k)]);
        end
        load(fname, 'res');
        if exist('res', 'var')
          if isfield(res.rasters{1}, 'pos') %#ok<NODEF>
            P.raster(res.indexes) = res.rasters(1:end);
          end
          clear res;
        end
      end
    end
    if ~sameranges
      newres.rasters = P.rasters;
      for k=1:numel(terclus.numRanges)
        thisRange = newres.rangeid==k;
        res = struct('Generation', {newres.Generation(thisRange)}, 'rangeid', {newres.rangeid(thisRange)}, 'indexes', {newres.indexes(thisRange)}, 'indexNewG', {newres.indexNewG(thisRange)}, 'rasters', {newres.rasters(thisRange)});
        save([params.nomdir sprintf('lastFrames%04d_%d.mat', res.Generation(1), k)], 'res');
        clear res thisRange;
      end
      clear newres;
    end
    doStart         = false;
  elseif iscell(params.initialPopulation)
    P.genome(1:end) = params.initialPopulation(1:end);
    doStart         = true;
  else
    error('params.initialPopulation must be a struct or a cell, but it isn''t any of these!!!');
  end
  params            = rmfield(params, 'initialPopulation');
else
  P.genome(1:end)   = params.generarGenomas(params.tam_poblacion, params);
  doStart           = true;
end

if doStart
  fprintf('PLEASE WAIT WHILE INITIAL POPULATION IS EVALUATED...\n');
  [P terclus params] = startPopulation(params, terclus, P, params.initialGeneration);
else
  fprintf('INITIAL POPULATION ALREADY EVALUATED...\n');
end
