function params = mainGeneric(initparams, postparams, parmode, mode, varargin)

%cosas parametrizadas que el algoritmo espera encontrar
%P.fitness=>fitness de cada individuo
%P.rndSeed=>semilla aleatoria de cada individuo
%P.raster=>informacion adicional de cada individuo, entendida como
%          necesaria para visualizarlo. Si es un struct con el campo
%          'isError', señaliza un error en la evaluación del individuo, que
%          para el algoritmo
% params.indivStruct    => cell array con los campos a almacenar en P,
%                          junto con sus valores por defecto
%
% params.genome2str         => convierte el genoma en un string
% params.writeMejorFun      => escribe el mejor de la poblacion
% params.writeIndividualFun => escribe un individuo de la poblacion
% params.writeGenealogyFun  => escribe la genealogía de un individuo
% params.multiEvFun  => función que evalúa varios individuos. Debe aceptar
%                       una serie de parámetros especificados, y devolver
%                       como salida los elementos correspondientes de
%                       indivStruct (menos .genome)
% params.generarGenomas => función que genera genomas aleatorios
% params.mutateGenome   => función que muta un genoma
% params.makeupParams   => función que permite alterar los parámetros,
%                          quizá copiándolos de un lugar a otro de la
%                          estructura.
%

MAINSUB = 'resultados';

if ~exist('parmode', 'var')
  parmode = 'embedded';
end

if iscell(mode)
    MAINSUB = mode{1};
    mode    = mode{2};
end

MAINSUBD = MAINSUB;
if ~isempty(MAINSUBD)
  MAINSUBD = [MAINSUBD filesep];
end

switch mode
  case 'replay' % additionalArgs = {subdirectory,newsubdirectory}
    subdirectory    = varargin{1};
    newsubdirectory = varargin{2};
    estadofile      = [MAINSUBD,subdirectory,filesep,'estado.mat'];
    load(estadofile,'-mat');      % carga los parámetros del fichero estado
%    params = putDefaultParameters(params); %#ok<NODEF>
    params = modifyParams(params, [], varargin(3:end)); %#ok<NODEF>
    rand(params.randmethod, hex2num(params.seedRandom)); %#ok<NODEF>
    params.nomdir = [MAINSUBD filesep newsubdirectory];
    if exist(params.nomdir, 'dir')
        rmdir(params.nomdir, 's');
    end
    mkdir(params.nomdir);
    save([params.nomdir,filesep,'estado.mat'], 'params','-mat');
  case 'continue' % additionalArgs = {subsubdirectory,newsubdirectory,presion,numgenerations}
    subsubdirectory = varargin{1};
    subdirectory    = [subsubdirectory filesep '..'];
    newsubdirectory = varargin{2};
    presion         = varargin{3};
    numgenerations  = varargin{4};
    estadofile      = [MAINSUBD,subdirectory,filesep,'estado.mat'];
    load(estadofile,'-mat');      % carga los parámetros del fichero estado
    params.nomdir = [MAINSUBD newsubdirectory];
    if exist(params.nomdir, 'dir')
        rmdir(params.nomdir, 's');
    end
    mkdir(params.nomdir);
    poph                          = loadSimulation([MAINSUBD,subsubdirectory], false, false, false);
    lastGen                       = max(poph.generation);
    thisGen                       = poph.generation==lastGen;
    params.initialGeneration   = lastGen;
    params.generaciones        = numgenerations;
    params.PRESIONES           = presion;
    params.num_iteraciones     = numel(presion);
    params.seedRandom          = num2hex(makeRandomSeed);
    params.randmethod          = 'twister';
    save([params.nomdir,filesep,'estado.mat'], 'params','-mat');
    if isfield(params, 'rewindNans') && params.rewindNans
      error('The code for continueing simulations is not prepared for this option!!\nIn order to get it prepared, individuals with NaN fitness should be tracked to previous generations in the poph structure');
    end
    params.initialPopulation   = struct;
    pophnames                  = fieldnames(poph);
    for k=1:poph.numFieldsPop %iterate for fields of the individuals
      fname = pophnames{k};
      params.initialPopulation.(fname) = poph.(fname)(thisGen);
    end
    params.initialPopulation.numFieldsPop = poph.numFieldsPop;
    params.initialPopulation.basedir = poph.basedir;
    params.initialPopulation.oldparams = poph.params;
    params.tam_poblacion       = numel(params.initialPopulation.(fname));
    params.terclus.ranges      = calculateRanges(params.tam_poblacion, params.terclus.numRanges);
    clear thisGen poph;
%    params = putDefaultParameters(params); %#ok<NODEF>
    params = modifyParams(params, postparams, varargin(5:end));
    rand(params.randmethod, hex2num(params.seedRandom)); %#ok<NODEF>
  case {'newsim', 'newsimSeeded'} % additionalArgs = pairs of strings 'option', 'value'
    params                     = initparams; %the default params structure
    clear                        initparams;
    switch mode
      case 'newsim'
        otherArgs = varargin;
      case 'newsimSeeded'
        seed = varargin{2};
        otherArgs = varargin(2:end);
    end
    params.seedRandom          = num2hex(makeRandomSeed);
    params.randmethod          = 'twister';
%---------------DEFINE PARAMETERS-----------------------------
    % PARAMETROS GENERALES
    params.saveRnd             = true;
    params.initialGeneration   = 0;
    params.PRESIONES           = 0.01;       %[1 0.2 0.1 0.0001 0.5];% [0.75 0.5 1/3 0.2 0.1 0.0001];
    params.num_iteraciones     = numel(params.PRESIONES);
    params.ensayos             = 1; %nº de veces que se hace el experimento para cada presion
    params.generaciones        = 500;%1000;%500;         % número de intentos para encontrar un árbol de una relación de aspecto dada
    params.tam_poblacion       = 40;45;100;45;%90;%1000;%300;         % número de individuos en la población
    params.presion             = []; %0.1;        % presion selectiva (0: no ha presión, 1: escalado normal)
    [params.nomdir params.ahora] = directorio(MAINSUB);
    params.tooLongEvTime       = inf;
    
    params.getRasters                    = false;
    params.plotting                      = true; 
    params.evaluateAlways                = true;
    params.logInEval                     = false;
    
    params.elitism = false;
    
%     params.evparams.methodToStop         = 0; %-1 (stop only by cpu), 0 (stop either by cpu or by time spent), 1 (stop only by time spent): 
%     params.evparams.maxCPUTimeEval       = 120;%300; %if empty, use as stop criteria fixedTimeSpent
%     params.evparams.fixedTimeSpent       = 500;
%     params.evparams.surrealDisp          = 1e2;
%     params.evparams.retrieve             = 'ss';
%     params.evparams.imageDimensions      = [500 500];
%     params.evparams.colorMap             = [0 0 0; 1 1 1; 0.5 0.5 0.5; 0 1 0; 0 0 1; 1 0 0];
%     params.evparams.tmpdirs              = {'c:\temp\', '/var/tmp/'};
%     [params.evparams.ss params.symbolS]  = initialSSConfiguration;
%     params.evparams.ss.rad               = params.evparams.ss.rad*0.5;
%     params.evparams.ss.stick.dist        = params.evparams.ss.stick.dist*0.5;
%     params.evparams.ss.allr              = params.evparams.ss.rad+params.evparams.ss.stick.dist;
% 
%     params.symbolS.symbols               = genesSwitchStrNum(params.symbolS.symbols);
%     params.symbolS.symbolOpen            = genesSwitchStrNum(params.symbolS.symbolOpen);
%     params.symbolS.symbolClose           = genesSwitchStrNum(params.symbolS.symbolClose);
% %   ss.organism.geneTypes                = {'mitosis',   'endblock', 'wait', 'fasterClock', 'slowerClock', 'moreDrag', 'lessDrag', 'morePenK', 'lessPenK', 'moreCycFreq', 'lessCycFreq', 'modeNeutral', 'modeCyclic', 'modeExcitable', 'excite'};
% %   ss.organism.geneAcronyms             = {'MT[', 'MT]', 'HLD', 'CK+', 'CK-', 'DG+', 'DG-', 'PK+', 'PK-', 'CF+', 'CF-', 'MDN', 'MDC'};%, 'MDE', 'EXC'};
%     params.genomeGeneration.symbolsFreq  =               [100    100    100    100    100    100    100    100    100    100    100];%    100    100];
%     params.genomeGeneration.mitosisFreq  =  200;
%     params.genomeGeneration.meanLong     = 50;
%     params.genomeGeneration.meanAnid     = 0;
% 
%     % PARAMETROS DE ACI
%     %modes for insertion operator:
%     %     -'one_level': the operator can insert '+', '-' and 'G'
%     %                   everywhere, but '[]' only at root level, to prevent
%     %                   nesting
%     %     -'normal': the operator can insert '+', '-', 'G', and '[]'
%     %                everywhere
%     %     -insert_list: the operator can insert any string in he cell list
%     %                   "insert_list" everywhere 
%     %false, it cannot insert nested brackets
%     params.mode_op_insercion  = 'normal'; %'normal'; 'one_level'; {{'+', '-'}};
%     %if true, only insert symbols inside brackets
%     params.insertOnlyNested   = false; %true; %false;
%     %if true, the deletion operator can delete non-empty brackets. If
%     %false, it can only delete empty brackets
%     params.op_eliminacion_extended  = true;%false;
%     %the way to select a symbol for mutation: by symbol (first select a
%     %symbol, then select a position in the genome that has that symbol) or
%     %by position (first select a position in the genome, then delete it)
%     params.selectSymbolForMutation  = 'byPos'; %'byPos'; %'bySymbol';
%     
%     
%     params.maxLongs.maxInsert  = []; %this is the max long of sub-string for mutations that insert sub-strings
%     params.maxLongs.maxString  = []; %this is the max long of sub-string for mutations that insert sub-strings
%     params.maxLongs.penalizeFitness = true; %flag to appy penalities based on these values to the fitness function
%     params.maxLongs.penalizeOperators = true; %flag to appy penalities based on these values to the results of mutation operators
%     params.elitism = false;%true;
%     params.lambdaMutation = 0.5;
%     params.mutationFreqs            =  [...
%         2        0.05  ;...      %      2   op_alteracion       A simbolo  posicion  
%         3        0.05  ;...      %      3   op_dupaleatoria      R segmento posicion
%         4        0.05  ;...      %      4   op_dupnivel          L segmento posicion
%         5        0.05  ;...      %      5   op_dupsecuencia      C segmento posicion
%         6        0.15  ;...      %      6   op_eliminacion       D posicion
%         7        0.05  ;...      %      7   op_insercion         I simbolo  posicion
%         8        0.00  ;...      %      8   op_transferencia     T segmento posicion
%                 ];   % probabilidad relativa de cada operador
    
    
    % PARAMETROS DE PARALELISMO
    params.terclus            = struct;
    params.terclus.doit       = true;
    params.terclus.tag        = 'springsNEW';
    params.terclus.jmArgs     = {'LookupURL',         'n0000'};
    params.terclus.jobArgs    = {'Tag',               params.terclus.tag, ...
                                 ...%'FileDependencies',  genpath, ...
                                 'PathDependencies',  genpath, ...
                                     };
    params.terclus.pauseTime           = 1;
    params.terclus.pauseTimeAfterError = 10;
    params.terclus.submitTimes         = 3;
    params.terclus.numRanges           = 1;
    params.terclus.ranges              = calculateRanges(params.tam_poblacion, params.terclus.numRanges);
    params.terclus.concJobs            = 1;%2; %number of concurrent jobs in the cluster
    params.terclus.concurrentWorkers   = 40;45;%48;%2; %number of concurrent workers in a job
    
    params.terclus.inputFileNameTemplate      = 'TaskInput%03d.mat';
    params.terclus.outputFileNameTemplate     = 'TaskOutput%03d.mat';
    params.terclus.evparamsFileNameTemplate   = 'evparams.mat';
    params.terclus.maxPicassoMem              = '400mb';
    params.terclus.maxPicassoTime             = '6:00:00';
    params.terclus.scriptFileName             = 'launchScript';
    params.terclus.jobIDFileName              = 'jobIDFile.txt';
    params.terclus.remoteSSH                  = false;
    params.terclus.ssh.classpath              = ['GaGeneric' filesep 'picasso' filesep 'ssh' filesep 'ganymed-ssh2-build250.jar'];
    params.terclus.ssh.nameurl                = 'picasso-cluster.scbi.uma.es';
    params.terclus.ssh.username               = 'josedavid';
    params.terclus.ssh.password               = '';
    params.terclus.ssh.nomdir                 = '';
    params.terclus.ssh.useMatlab              = false;
    params.terclus.showOutput                 = false;
    params.terclus.ssh.rmCommand              = ['rm Task*.mat; rm ' params.terclus.evparamsFileNameTemplate '; rm ' params.terclus.scriptFileName];
    
    
    params = modifyParams(params, postparams, otherArgs);
    
    if strcmp(mode, 'newsimSeeded')
      if ~iscell(seed)
        seed = {seed};
      end
      switch numel(seed)
        case 1
          seed = repmat(seed, params.tam_poblacion, 1);
        case params.tam_poblacion
          seed = reshape(seed, [], 1);
        otherwise
          error('If there are more than one seed, the population size and the number of seed must match!!!');
      end
      params.initialPopulation = seed;
      clear seed otherArgs;
    end

    rand(params.randmethod, hex2num(params.seedRandom));
    fprintf('seedRandom = %s\n', params.seedRandom);
    
    save([params.nomdir,filesep,'estado.mat'], 'params','-mat');
%---------------END DEFINE PARAMETERS-------------------------
    otherwise
        error('mode not understood <%s>', any2str(mode));
end



printReadableEstado([params.nomdir,filesep,'readableEstado.txt'], params, 'params', 'w', 'terclus.ssh.password');

%place it here to not save it to estado.mat
switch parmode
    case 'local'
      params.terclus.jobMgr = findResource('scheduler', 'type', 'local');
      params.terclus.concurrentWorkers = params.terclus.jobMgr.ClusterSize;
    case 'remote'
      params.terclus.jobMgr = findResource('scheduler', 'type', 'torque', params.terclus.jmArgs{:});
      %params.terclus.jobMgr = findResource('scheduler', 'type', 'jobmanager', params.terclus.jmArgs{:});
    case 'embedded'
      params.terclus.jobMgr = [];
    case 'picasso'
      params.terclus.jobMgr = 'picasso';
      if params.terclus.remoteSSH
        if params.terclus.pauseTime<60
           params.terclus.pauseTime=60;
        end
        javaaddpath(params.terclus.ssh.classpath);
%   javaaddpath([params.terclus.SSHjavadir, 'jsch.jar']);
%   javaaddpath([params.terclus.SSHjavadir, 'sshutils.jar']);
      end
    otherwise
      error('First argument must be either ''local'' or ''remote'', but it is ''%s''!!!\n', num2str(parmode));
end


if isfield(params, 'GAFun')
  gafun = params.GAFun;
else
  gafun = @GAGeneric;
end

for iter=1:params.num_iteraciones
    params.num_iteracion= iter;
    params.presion = params.PRESIONES(iter);
    fprintf('Iteration nº %g, pressure %.10f\n', iter, params.presion);
    
    for i = 1:params.ensayos
        params.ensayo       = i;
        fprintf('ENSAYO %g\n', i);
        subnomdir = gafun(params); %#ok<NASGU>
%         makeAllGraphics(subnomdir, true);
%         makeLongGenomeGraphics(subnomdir, inf, true);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seed = makeRandomSeed
seed = sum(100*clock);
%clk  = clock; seed = mod(floor(sum(clk([1 6 5 4 3 2]).*[1e10 1e8 1e6 1e4 1e2 1])), 2^32);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nomdir, ahora] = directorio(nombre)

if ~isempty(nombre)
  nombre = [nombre filesep];
end

ahora = datestr(now, '_yyyy_mm_dd_HH_MM_SS');
ahora(ahora=='_') = 'YMDhms';
nomdir = [nombre,ahora];
mkdir(nomdir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = modifyParams(params, postparams, args)
if isstruct(postparams)
  params = replaceParams(params, postparams);
end
k=1;
while k<numel(args)
  if ~ischar(args{k})
    error('I expect ''option'', ''value'' pairs!!!');
  end
  try
    eval(['params.' args{k} ';']);
    if ischar(args{k+1})
      eval(['params.' args{k} ' = ' args{k+1} ';']);
    elseif iscell(args{k+1}) && (numel(args{k+1})==1)
      points = find(args{k}=='.');
      if isempty(points)
        params.(args{k})                                   = args{k+1}{1};
      else
        starts = [1 (points+1)];
        ends   = [(points-1) numel(args{k})];
        strs   = arrayfun(@(s,e) args{k}(s:e), starts, ends, 'uniformoutput', false);
        switch numel(starts)
          case 2
            params.(strs{1}).(strs{2})                               = args{k+1}{1};
          case 3
            params.(strs{1}).(strs{2}).(strs{3})                     = args{k+1}{1};
          case 4
            params.(strs{1}).(strs{2}).(strs{3}).(strs{4})           = args{k+1}{1};
          case 5
            params.(strs{1}).(strs{2}).(strs{3}).(strs{4}).(strs{5}) = args{k+1}{1};
          otherwise
            error('error');
        end
      end
    else
      error('error');
    end
    if any(strcmp(args{k}, {'tam_poblacion', 'terclus.numRanges'}))
      params.terclus.ranges = calculateRanges(params.tam_poblacion, params.terclus.numRanges);
    end
    if strcmp(args{k}, 'statNamesFull') && isfield(params, 'buildFieldStatsFun')
      params = params.buildFieldStatsFun(params);
    end
  catch ME
    error('Sorry, but the pair <%s>, <%s> canot be evaluated!!!', args{k}, any2str(args{k+1}));
  end
  k=k+2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = replaceParams(params, postparams, flag)
fn = fieldnames(postparams);
for k=1:numel(fn)
  fk = fn{k};
  if isstruct(postparams.(fk))
    params.(fk) = replaceParams(params.(fk), postparams.(fk), strcmp(fk, 'terclus'));
  else
    params.(fk) = postparams.(fk);
    if strcmp(fk, 'tam_poblacion') || ...
       ( (nargin>2) && flag && strcmp(fk, 'numRanges') )
      params.terclus.ranges = calculateRanges(params.tam_poblacion, params.terclus.numRanges);
    end
  end
end
