function subnomdir = GAStepByStepGeneric(params)
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

if ~isfield(params, 'saveRnd')
  params.saveRnd = true;
end

terclus.elapsed = 0;

params.fpob = fopen([params.nomdir,'poblacion.txt'],'a');
params.farb = fopen([params.nomdir,'arbol.txt'],'a');
params.fres = fopen([params.nomdir,'resultado.txt'],'a');
params.ftst = fopen([params.nomdir,'timestats.txt'],'a');
params.allfs = {params.fpob; params.farb; params.fres; params.ftst};


params.writeMejor            = params.writeMejorFun(params.fres, params.genome2str);
params.writeIndividual       = params.writeIndividualFun(params.fpob, params.genome2str);
params.writeGenealogy        = params.writeGenealogyFun(params.farb);
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
params.fieldsToFill = 2:numel(fieldnames(params.indivStruct)); %skip first field 'genome' when reaping results
%it is convenient to pre-fabricate this structure
jobStruct = {'G' 'indexNewG' 'jobs' 'generation' 'isInitial'};
% legend: G==individuals of this range
%         indexNewG== indexes in G of truly new individuals (only they will be evaluated)
%         job== job handle in charge of evaluating this range
%         generation== the generation
jobStruct = reshape([jobStruct; cell(size(jobStruct))], 1, []);
params.jobStruct = struct(jobStruct{:}); clear jobStruct;

if ~isempty(params.makeupParams)
  params = params.makeupParams(params);
end

%here, we will store the best individual of the simulation. It cannot be
%safely stored in the population, since it might be overwritten
mejorEnConserva = struct('gen', [], 'num', [], 'data',  {params.indivStruct});
jobsActive = cell(0,1);
try
    terclus.rangeSizes = diff(terclus.ranges,1,2)+1;

    P = repmatStruct(params.indivStruct, [params.tam_poblacion,1]);
    
    [P params terclus] = prepareInitialPopulation(P, params, terclus);
    
    for i=1:numel(P.genome)
      params.writeGenealogy(params.initialGeneration, i, '=', 0, 0);
      params.writeIndividual(params.initialGeneration, i, P);
    end
    
    %find the best in the population
    mejorEnConserva = findMejor(P, mejorEnConserva, terclus.ranges, params.initialGeneration);
    fprintf('TIME ELAPSED: %g\n', etime(clock, startTime));

    doTermination = false;
    lastTime = clock;
    rind = 1;
    for cont=(params.initialGeneration+1):(params.initialGeneration+params.generaciones) %for each generation
      %for each range, the generation it was last updated
      generationInRange   = cont-1;
      nj                  = params.jobStruct;
      nj.generation       = cont;
      nj.isInitial        = false;
      if params.mutateOutside
        [P nj.indexNewG]  = evolucionGeneric(params, P,cont,rind,generationInRange,terclus,mejorEnConserva); % proxima generacion
        longs             = cellfun('prodofsize', P.genome(nj.indexNewG));
        mprintf(params.ftst, 'SPAWNED GEN=%g, Nº OF STRINGS: %g\n   TIME ELAPSED=%f (CPU: %f), TIME IN CLUSTER=%f\n', cont, numel(nj.indexNewG), etime(clock, startTime), cputime-initcputime, terclus.elapsed);
        mprintf(params.ftst, '   NEW STRING LENGTH   STATS: MIN=%g, MAX=%g, MEAN=%g, MEDIAN=%g\n', min(longs), max(longs), mean(longs), median(longs));
      end
      %nj.job        = evaluarRemote(terclus, nj.G.genome(nj.indexNewG), params, nj.generation, nj.rangeid, nj.indexNewG);
      [params terclus P nj] = params.evaluateGeneration(params, terclus, P, nj);
      mprintf(params.ftst, 'FINISHED GENERATION %d\n   TIME ELAPSED: %f (CPU: %f), TIME IN CLUSTER: %f\n', cont, etime(clock, startTime), cputime-initcputime, terclus.elapsed);
      tooSlowEvaluation = etime(clock, lastTime)>params.tooLongEvTime;
      if tooSlowEvaluation
        mprintf(params.ftst, 'We stopped prematurely!!! tooSlowEvaluation=%s\n', any2str(tooSlowEvaluation));
        doTermination = true;%break;
      end
      lastTime = clock;
      [P, mejorEnConserva] = updateGA(params, terclus, nj, P, mejorEnConserva, startTime);
      if doTermination
        break;
      end
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
nj            = params.jobStruct;
nj.generation = Generation;
nj.indexNewG  = 1:numel(P.genome);
nj.isInitial  = true;

[params terclus P nj] = params.evaluateGeneration(params, terclus, P, nj);

if(params.plotting)
  writePNGs(params, P, Generation, 1:numel(P.genome));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mejorEnConserva = findMejor(P, mejorEnConserva, ranges, generationInRange)
[maxi, mejor] = max(P.fitness);               %#ok<ASGLU> % mejor individuo
% mejor = noplanos(mejor);
if (~isempty(mejor)) %&& ( isempty(mejorEnConserva.data.fitness) || (maxi>=mejorEnConserva.data.fitness) )
    mejorEnConserva.gen   = generationInRange;
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
function [P, mejorEnConserva] = updateGA(params, terclus, nj, P, mejorEnConserva, startTime)
ranges = terclus.ranges;

    gen   = nj.generation;
    
    areError = find(cellfun(@(x)not(isempty(x)) && isfield(x, 'isError'), P.raster(nj.indexNewG)));
    if not(isempty(areError))
      errs = P.raster(nj.indexNewG(areError));  %#ok<NASGU>
      idxs = nj.indexNewG(areError); %#ok<NASGU>
      save([params.nomdir sprintf('errors_%04d.mat', gen)], 'errs', 'idxs');
      error('There has been an error while evaluating at least one individual!!!!');
    end
    
    if(params.plotting)
      writePNGs(params, P, gen, nj.indexNewG);
      %P.raster(ranges(rpri,:)) = {[]};
    end
        
    %recalculate best individual
    mejorEnConserva = findMejor(P, mejorEnConserva, ranges, gen); %#ok<ASGLU>
    %output stage
    params.writeMejor(mejorEnConserva); %write best individual for this generation
    for z=1:numel(P.genome)
      params.writeIndividual(gen, z, P);
    end
    fprintf(params.ftst, '   POPULATION GENERATIONS UPDATE: GENERATION=%g    (TIME ELAPSED: %s)\n', gen, mat2str(etime(clock, startTime)));
    
    if params.saveRnd
      frnd = fopen([params.nomdir 'rndseed.txt'], 'w');
      fprintf(frnd, 'Generation %04d, rndseed for ''%s'': uint32(', gen, params.randmethod);
      fprintf(frnd, mat2str(rand(params.randmethod)));
      fprintf(frnd, ')\n');
      fclose(frnd);
    end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writePNGs(params, G, generation, indexNewG)
indexes = indexNewG;

res = struct('Generation', {generation}, 'indexes', {indexes}, 'indexNewG', {indexNewG}, 'rasters', {G.raster(indexNewG)}); %#ok<NASGU>
save([params.nomdir sprintf('lastFrames%04d.mat', generation)], 'res');


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
