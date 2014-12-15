function conf = molMotConf(seedVersion, pophVersion, fitnessVersion)

if (~exist('seedVersion', 'var')) || isempty(seedVersion)
  seedVersion = 1;
end
if (~exist('pophVersion', 'var')) || isempty(pophVersion)
  pophVersion = 3;
end
if (~exist('fitnessVersion', 'var')) || isempty(fitnessVersion)
  fitnessVersion = 1;
end

geneTypes              = {...
  'MSS', 'mitosisSpringShift'; ...
  'MST', 'mitosisStatic';...
  'SWF', 'switchForcement';...%¡¡¡!!!switchForcement(spring:1:6 ¿relative to springtoshift?, both:logical, #generations, alsoSwitchType:logical¿¿¿???)
  ...%¡¡¡!!!switchForcementType(spring:1:6 ¿relative to springtoshift?, both:logical, #generations)
  ...%'permuteGeneDistribution', ... %arg 1: 3 options
  ...%'permuteChangeToDevelopment', ... %arg 1: 6 options, arg 2: 5 options
  ...%'increaseR', ...
  ...%'decreaseR', ...
  ...%'increaseK', ...
  ...%'decreaseK', ...
  ...%V1 & V2: a column is swtiched. This preserves the same springs to shift, but the springs are glued by the opposite ends  
  ...%D1 & D2: a diagonal of the 2x2 matrix is switched. This preserves the oppositeness, and so is better than H1 & H2 (rows are switched)
  ...%V: the same as concatenating V1V2: both columns are switched (thus, the rows are switched). this is to enhance expressiveness of the code
  ...%H: the same as concatenating D1D2: both diagonals are switched (thus, the columns are switched). This is a very important operator: it handles handedness
  'PSM', 'permuteShiftMatrix'; ... %arg 1: V1, V2, D1, D2, H, V
  'MSP', 'mirrorSprings'; ...
  'CGD', 'changeGeneDistribution'; ... %arg 1: -1, 0, 1, to mean 0, 0.5, 1
  'CSD', 'changeAttachedSpringsDistribution'; ...
  'CMT', 'changeMitosisTime'; ...
  'WST', 'waitSomeTime'; ...
  'CSR', 'changeSpringR'; ...
  'CSK', 'changeSpringK'; ...
  'CET', 'changeEdgeTime' ...
  };
name2Abrev   = reshape([geneTypes(:,2)'; geneTypes(:,1)'], 1, []);
abrev2Name   = reshape([geneTypes(:,1)'; geneTypes(:,2)'], 1, []);
geneTypesStr = reshape([geneTypes(:,1)'; arrayfun(@(x)x, 1:size(geneTypes,1), 'UniformOutput', false)], 1, []);
geneTypesStr = struct(geneTypesStr{:});
geneTypes    = reshape([geneTypes(:,2)'; arrayfun(@(x)x, 1:size(geneTypes,1), 'UniformOutput', false)], 1, []);
geneTypes    = struct(geneTypes{:});
geneTables   = struct('name2Abrev', struct(name2Abrev{:}), 'abrev2Name', struct(abrev2Name{:}), 'name2num', geneTypes, 'abrev2num', geneTypesStr);

conf.params.genome2str          = convertGenomeStrSimple(geneTables);
conf.params.evaluateGenome      = @evaluarMolMot;
conf.params.multiEvFun          = @evaluarMolMots;
conf.params.writeMejorFun       = @DAwriteMejorFun;
%rows:
%       1: variable/field name
%       2: default value
%       3: write format in poblacion.txt (empty string means omission)
%       4: order in poblacion.txt
%       5: transform function for writing in poblacion.txt
%       6: read format from poblacion.txt
conf.params.statNamesFull           = {'fitness', 'rndSeedDev', 'rndSeedEv', 'devTimeSpent', 'evTimeSpent', 'cpuTime', 'ncells', 'ncellsD', 'ncellsTeo', 'npoints', 'nsprings1', 'nsprings2', 'developed', 'overload1', 'overload2', 'atpChg1', 'atpChg2', 'offsetAbs', 'offsetRel'  ; ...
                                       nan,       nan,          nan,         nan,            nan,           nan,       nan,      nan,       nan,         nan,       nan,         nan,         false,       nan,         nan,         nan,       nan,       nan,         nan          ; ...
                                       '%17s',    '%s',         '%s',        '%-06.3f',      '%-06.3f',     '%-06.3f', '%03d',   '%03d',    '%03d',      '%03d',    '%03d',      '%03d',      '%01d',      '%13d',      '%13d',      '%02d',    '%02d',    '%13d',      '%13d'       ;...
                                       01,        06,           07,          03,             04,            05,        08,       09,        10,          11,        12,          13,          14,          15,          16,          17,        18,        19,          20,          ; ...
                                       @mat2str,  @num2hex,     @num2hex,    '',             '',            '',        '',       '',        '',          '',        '',          '',          '',          '',          '',          '',        '',        '',          ''           ; ...
                                       '%f',      '%s',         '%s',        '%f',           '%f',          '%f',      '%f',     '%f',      '%f',        '%f',      '%f',        '%f',        '%f',        '%f',        '%f',        '%f',      '%f',      '%f',        '%f'        };                                                                                      
if pophVersion >= 2
    conf.params.statNamesFull = [conf.params.statNamesFull, {'timeAttached'; nan; '%7d'; 21; ''; '%f'}];
  if pophVersion >= 3
    conf.params.statNamesFull = [conf.params.statNamesFull, {'spectralGap';  nan; '%7d'; 22; ''; '%f'}];
  end
end
conf.params.indivStructFull         = [{'genome', 'rasterDev', 'rasterEv'  ; ...
                                       {[]},     {[]},        {[]}         ; ...
                                       '%s',     '',          ''           ; ...
                                       inf,      nan,         nan          ; ...
      convertGenomeStrWrapped(geneTables),       '',          ''           ; ...
                                       '%[^\n]', '',          ''          }, ...
                                       conf.params.statNamesFull];
conf.params.evparams.statNames   = conf.params.statNamesFull(1:2,:);
conf.params.indivStruct          = conf.params.indivStructFull(1:2,:);
conf.params.writeIndividualFun   = makeWriteFun(conf.params);
conf.params.writeGenealogyFun    = @DAwriteGenealogyFun;
[conf.params.genome,            ...
 conf.params.mutateGenome,      ...
 conf.params.makeupParams,      ...
 conf.params.generarGenomas]    = DAmakeMutationConfiguration(geneTypes);

conf.params.selection.mode         = 'FUSS'; %'FUSS'; 'PACO';
conf.params.selection.FUSS.epsMode = 'rel'; %'abs'; 'rel';
conf.params.selection.FUSS.epsDist = 0.05;
conf.params.selection.FUSS.mode    = 'levelled1'; %'levelled1'; 'plain';
conf.params.selection.FUSS.ntimes  = 2; 1; 3;
conf.params.rewindNans             = true;%true; false;
conf.params.newsQuotaMode          = 'plain'; %'none';  'plain'; 'rand';      'custom';
conf.params.newsQuota              = 0.15; %NOTHING; 0.3;     [0.25 0.35]; @(generacion, rangeid, generationInRange, range, P, params)ceil(rand*10);

conf.fun.system                 = AtomSystemSimple(geneTypes);

conf.ssParams.geneTables        = geneTables;
conf.ssParams.geneTypes         = geneTypes;
conf.ssParams.tickWidth         = 0.5;
conf.ssParams.stepsByTick       = 50;
conf.ssParams.ndims             = 2;
conf.ssParams.U                 = 1;5;
conf.ssParams.rad               = 2;
conf.ssParams.allr              = 4;
conf.ssParams.penk              = 100;
conf.ssParams.K                 = 20;
conf.ssParams.stretchedFactor   = 0.75;
conf.ssParams.compressedFactor  = 1.25;
if iscell(seedVersion)
  if numel(seedVersion)>1
    atomLength  = seedVersion{2};
  else
    atomLength  = [];
  end
  seedVersion   = seedVersion{1};
else
  atomLength    = [];
end
if isstruct(seedVersion)
  conf.samples.ss = seedVersion;
  switch numel(atomLength)
    case 0
      warning('molMotConf:atomLength', 'The value of atomLength is estimated to be %s, but, since the seed has been imposed, it is only a guess', any2str(atomLength));
      atomLength = conf.samples.ss.pos(2,2)-conf.samples.ss.pos(4,2);
    case 1
      warning('molMotConf:atomLength', 'The value of atomLength is assigned the value %s, but, since the seed has been imposed, it is the user''s responsability to ensure it', any2str(atomLength));
    case 2
      alpoints = atomLength;
      atomLength = realsqrt(sum(realpow(conf.samples.ss.pos(alpoints(1),:)-conf.samples.ss.pos(alpoints(2),:), 2), 2));
      warning('molMotConf:atomLength', 'The value of atomLength is assigned the distance %s from the points %s, but, since the seed has been imposed, it is the user''s responsability to ensure it', any2str(atomLength), any2str(alpoints));
    otherwise
      error('jarllll!');
  end
else
  y = 1.3;
  switch seedVersion
    case {1, 2}
      energyDens = sign(0.5-mod(seedVersion, 2));
      conf.samples.ss = generateSampleSeed1(conf, [0 1.25*y; 1 1.5*y; 0 0.25*y; 1 0*y], energyDens);
      atomLength = conf.samples.ss.pos(2,2)-conf.samples.ss.pos(4,2);
    case {3, 4}
      energyDens = sign(0.5-mod(seedVersion, 2));
      conf.samples.ss = generateSampleSeed1(conf, [0 1.5*y; 1 1.5*y; 0 0*y; 1 0*y], energyDens);
      atomLength = conf.samples.ss.pos(2,2)-conf.samples.ss.pos(4,2);
    case {5, 6}
      energyDens = sign(0.5-mod(seedVersion, 2));
      x = 1.2;
      y = x*1.75;
      conf.samples.ss = generateSampleSeed1(conf, [0*x 0*y; 1*x 0*y; 0.5*x sin(pi/3)*y; 0.5*x tan(pi/6)/2*y], energyDens);
      atomLength = realsqrt(sum(realpow(conf.samples.ss.pos(3,:)-conf.samples.ss.pos(1,:), 2), 2));
    otherwise
      error('seed version not understood!!!');
  end
end

conf.params.GAFun                                   = @GAGenericDecoupledDev;
conf.params.writeRangeID                            = true;

conf.params.evparams.verbose                        = false;
conf.params.evparams.numTicksForStabilization       = 2;
conf.params.evparams.methodToStop                   = 0; %0, -1, +1
conf.params.evparams.dev.fixedTimeSpent             = 500;
conf.params.evparams.dev.maxCPUTimeEval             = 100;
conf.params.evparams.dev.devEThreshold              = 0.01;
conf.params.evparams.walker.fixedTimeSpent          = 1000;
conf.params.evparams.walker.maxCPUTimeEval          = 500;
conf.params.evparams.walker.devEThreshold           = -inf;
conf.params.evparams.relax.fixedTimeSpent           = 100;
conf.params.evparams.relax.maxCPUTimeEval           = 100;
conf.params.evparams.relax.devEThreshold            = 0.01;

conf.params.jiggleSelection                         = eps(1);
conf.params.evparams.record.positions               = true;
conf.params.evparams.record.useWalks                = 'longest'; %'aggregate';
conf.params.evparams.fitnessCalc                    = 'mode2';%'plain';%'mode1';
conf.params.evparams.spectralGap.rotCutoff          = 1e-12;
conf.params.evparams.spectralGap.gapToRecord        = 4;
conf.params.evparams.imageDimensions                = 500;
conf.params.evparams.relax.doRelax                  = true;
conf.params.evparams.relax.holdOnRow                = true;
conf.params.evparams.walker.fixUnixWindows.doFix    = true;
conf.params.evparams.walker.fixUnixWindows.epsilon  = eps(1e1);
conf.params.evparams.walker.fixUnixWindows.maxilon  = eps(1e10);
conf.params.evparams.walker.minAmountOfToes         = [];%minAmountOfToesFun;%@(conns, vs) min(8, ceil(numel(vs)/2)); %[];
conf.params.evparams.walker.cutoff                  = 1.5;
conf.params.evparams.walker.brownianFactor          = 0;1; %-0.1; %for the moment, do not use brownian motion, as we are not able of stripping small structures from large fitness from brownian excursions instead of good processivity
conf.params.evparams.walker.rewiringPlace           = 'after'; %'none'; 'before'; 'after';
conf.params.evparams.walker.rewiringMode            = {'geotopological', 'georecurtopo'}; %'geotopological'; %'topological'
conf.params.evparams.walker.rewireLength            = [atomLength*2      atomLength/2];   %atomLength*2;     %nan;
conf.params.evparams.walker.ABCMode                 = 'heurconcave'; %'meanpos'; 'heurconcave';
conf.params.evparams.walker.tooSmallNeighbour       = 2/3;
conf.params.evparams.walker.springK                 = 100;
conf.params.evparams.walker.pointRad                = 3;
conf.params.evparams.walker.pointAllr               = 4;
conf.params.evparams.walker.row.ballRad             = 3;
conf.params.evparams.walker.row.ballAllr            = 4;
conf.params.evparams.walker.row.ballSep             = 8;
conf.params.evparams.walker.initialDepth            = 1;
conf.params.evparams.walker.row.ballPenk            = 500;
conf.params.evparams.walker.row.typeSpec            = [2 2];
conf.params.evparams.walker.row.active              = true;
conf.params.evparams.walker.row.allPoints           = false;
conf.params.evparams.walker.row.lockFreeze          = true;
%conf.params.evparams.walker.row.freezeK             = 1/conf.params.evparams.walker.row.ballPenk;
conf.params.evparams.walker.leg.dynLengths          = false;
conf.params.evparams.walker.leg.dynLengthFactor     = 0.5;
conf.params.evparams.walker.leg.atpK                = 200;
conf.params.evparams.walker.leg.timeAtp             = 50;
conf.params.evparams.walker.leg.timeRelax           = 20;
conf.params.evparams.walker.leg.stateInit           = {'sticky', 'relaxation'};
conf.params.evparams.walker.leg.stateInitTime       = [nan,      0];
conf.params.evparams.walker.correctMistakeRotation  = true;
conf.params.evparams.walker.correctMistakeSign      = true;
conf.params.evparams.walker.correctMistakeRewiring  = true;
conf.params.evparams.walker.eigCalc                 = 'full'; %'full'; 'sparse';

conf.fun.developGenome          = @DADevelopGenome;
conf.fun.gnmEquals              = @DAgenomesAreEqual;

conf.fun.loadSim                = makeLoadFun(@loadFun1, conf.params);
conf.fun.phylotree              = phyloFun(@getBasicInf1);

conf.samples.params             = conf.params;
conf.samples.params.evparams.ss = conf.samples.ss;

conf.samples.postparams.saveRnd                   = true;
conf.samples.postparams.PRESIONES                 = 0.01;       %[1 0.2 0.1 0.0001 0.5];% [0.75 0.5 1/3 0.2 0.1 0.0001];
conf.samples.postparams.ensayos                   = 1; %nº de veces que se hace el experimento para cada presion
conf.samples.postparams.generaciones              = 100;%1000;%500;         % número de intentos para encontrar un árbol de una relación de aspecto dada
conf.samples.postparams.tam_poblacion             = 100;45;100;45;%90;%1000;%300;         % número de individuos en la población
conf.samples.postparams.tooLongEvTime             = inf;
conf.samples.postparams.plotting                  = true; 
conf.samples.postparams.evaluateAlways            = true;
conf.samples.postparams.logInEval                 = false;
conf.samples.postparams.logInMutation             = false;
conf.samples.postparams.terclus.concurrentWorkers = conf.samples.postparams.tam_poblacion;40;45;%48;%2; %number of concurrent workers in a job
conf.samples.postparams.elitism                   = false;

conf.samples.mainfun                              = DAmainFun(conf.samples.params, conf.samples.postparams);
end

function fun = minAmountOfToesFun
  fun = @(conns, vs) min(8, ceil(numel(vs)/2));
end

function fun = DAmainFun(params, postparams)
fun = @(varargin)mainGeneric(params, postparams, varargin{:});
end

function cgs = convertGenomeStrSimple(geneTables)
cgs = @(genome) DAConvertGenome(genome, geneTables);
end

function cgs = convertGenomeStrWrapped(geneTables)
cgs = @(genome) DAConvertGenome(genome{1}, geneTables);
end

function as = AtomSystemSimple(geneTypes)
as = @(varargin) DynamicAtomSystem(geneTypes, varargin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = DAwriteMejorFun(fres, switchGenomeStr)
  fun = @(mejorEnConserva) ...
    fprintf(fres,'%g %g %g %s %s %s %s\n',mejorEnConserva.gen, mejorEnConserva.rind, mejorEnConserva.num, num2hex(mejorEnConserva.data.rndSeedDev), num2hex(mejorEnConserva.data.rndSeedEv), mat2str(mejorEnConserva.data.fitness), switchGenomeStr(mejorEnConserva.data.genome{1})); % actualizamos el fichero de resultados
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = DAwriteGenealogyFun(farb)
  fun = @(generationNew, rangeindNew, indexNew, changeStr, generationOld, rangeindOld, indexOld) ...
    fprintf(farb,'%d %d %d %s %d %d %d\n',generationNew,rangeindNew,indexNew,changeStr,generationOld, rangeindOld, indexOld);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ssSeed = generateSampleSeed1(conf, initpos, fixedForceDensity)
K                 = conf.ssParams.K;%20;
stretchedFactor   = conf.ssParams.stretchedFactor;%0.75;
compressedFactor  = conf.ssParams.compressedFactor;%1.25;
% initpos           = [0 0; 0 1; 1 0; 1 1];
%initpos      = [0 1.5; 1 1.5; 0 0; 1 0]*10;%[0 1.25; 1 1.5; 0 0.25; 1 0]*10;
%initpos      = [0 1.25; 1 1.5; 0 0.25; 1 0]*10;
initpos = initpos*10;
%initpos(:,2) = initpos(:,2)*1.3;
% %pos                   = [1 0; -1 0; 0 sqrt(3); 0 1/sqrt(3)]*10; %triangled
% %pos                   = [0 1; 1 1; 0 -0.5; 1 -0.5];
% %pos                   = [0 1; 1 1; 0 0; 1 0];
% %pos                   = [0 1; 1 1; 0 -10; 1 -10];
springToFix       = 1;
%fixedForceDensity = -1;%+1;%-1;

ss = conf.fun.system(conf.ssParams.ndims, 0, conf.ssParams.U, conf.ssParams.tickWidth, conf.ssParams.stepsByTick);

ss.pointDefaultVals(4:6) = {conf.ssParams.rad, conf.ssParams.allr, conf.ssParams.penk};

sspec                 = springSpec('K', stretchedFactor, compressedFactor, K);
[springsR, springsK]  = createAtomSprings4P1S(initpos, springToFix, fixedForceDensity, sspec);

[ss newAtom]  = addAtomComplex(ss, zeros(1,4), initpos, springsK, springsR, {[]});
ss.lockedUntil(1) = ss.t;

ssSeed = ss;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = makeLoadFun(loadFun, params)
fun = @(varargin) loadFun(params, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function poph = loadFun1(params, varargin)
indivStruct   = params.indivStructFull; %get full info on fields
names         = indivStruct(1,:); %extract info
order         = cell2arrayMX(indivStruct(4,:));
reads         = indivStruct(6,:);
readed        = cellfun('prodofsize', reads)>0; %find the fields to be read
[order order] = sort(order); %order the field according to the order in the file
names         =   names(order);
reads         =   reads(order);
readed        =  readed(order);
names         =   names(readed); %filter out non-written fields
reads         =   reads(readed);
reads         = [reads; array2cell(repmat(' ', size(reads)))]; %compose scan string
reads{end}    = '';
reads         = horzcat('%f %f %f ', reads{:});
names         = [{'generation', 'rangeid', 'individual'} names];
alsoRange     = true;
poph          = loadGeneric(reads, names, alsoRange, varargin{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basicinf = getBasicInf1(poph, idx)
basicinf = sprintf('Index in struct: %04d, Generation: %03d, individual: %03d, ncellsD: %d, np: %d, ns: %d,\nFITNESS: %08.6f, chg1: %d, chg2: %d, overload1: %d, overload2: %d,\noffsetA: %d, offsetR: %d devTimeSpent: %04d, evTimeSpent: %04d, cpuTimeSpent: %04d', ...
                   idx, ...
                   poph.generation(idx), ...
                   poph.individual(idx), ...
                   poph.ncellsD(idx), ...
                   poph.npoints(idx), ...
                   poph.nsprings1(idx), ...
                   poph.fitness(idx), ...
                   poph.atpChg1(idx), ...
                   poph.atpChg2(idx), ...
                   poph.overload1(idx), ...
                   poph.overload2(idx), ...
                   poph.offsetAbs(idx), ...
                   poph.offsetRel(idx), ...
                   poph.devTimeSpent(idx), ...
                   poph.evTimeSpent(idx), ...
                   poph.cpuTime(idx));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fun = phyloFun(getBasicInfFun)
fun = @(varargin)phyloGeneric(getBasicInfFun, @printLastFrame, varargin{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printLastFrame(poph, getBasicInfFun, indDisplayed, indRequested)
gid = poph.generation(indDisplayed);
iid = poph.individual(indDisplayed);
gir = poph.generation(indRequested);
iir = poph.individual(indRequested);
if isfield(poph, 'rangeid')
  rid = poph.rangeid(indDisplayed);
  fname = [poph.basedir filesep sprintf('lastFrames%04d_%d.mat', gid, rid)];
else
  fname = [poph.basedir filesep sprintf('lastFrames%04d_1.mat', gid)];
end
drawn = false;
fuseSprings = false;
margen = 1.1;
pltargs = {'axisEqual', false, 'axisSquare', true, 'circleFaces', 50};
if exist(fname, 'file')
  load(fname, 'res');
  if exist('res', 'var')
    stDev = res.rasterDevs{res.indexes==iid}; 
    stEv  = res.rasterEvs{res.indexes==iid}; 
    if all(isfield(stDev, {'pos', 'springEnds'}))  && all(isfield(stEv, {'pos', 'springEnds'}))
      stEv.row = poph.params.evparams.walker.row;
      asp = BallPlotter('axisWindow', getAxisWindow(stDev.pos, margen), pltargs{:});
      if fuseSprings
        t = stDev.t; stDev = fuseRepeatedSprings(combineBasicSystems(DynamicAtomSystem([], {}), stDev), false); stDev.t = t;
      end
      drawSSAndText(asp, stDev, getBasicInfFun(poph, indRequested));
      assignin('base', 'dss', stDev);
      %maximize(gcf);
      asp2 = BallPlotter('axisWindow', getAxisWindow(stEv.pos, margen), pltargs{:});
      if fuseSprings
        t = stEv.t; stEv = fuseRepeatedSprings(combineBasicSystems(DynamicAtomSystem([], {}), stEv), false); stEv.t = t;
      end
      drawSSAndText(asp2, stEv, getBasicInfFun(poph, indRequested));
      %maximize(gcf);
      drawn = true; %#ok<NASGU>
      assignin('base', 'zss', stEv);
    else
      fprintf('UNABLE TO SHOW INDIVIDUAL AS THE DATA HAS NOT BEEN UNDERSTOOD!!!!\n');
    end
  else
    fprintf('UNABLE TO SHOW INDIVIDUAL AS THE VARIABLE res HAS NOT BEEN FOUND!!!!\n');
  end
else
  fprintf('UNABLE TO SHOW INDIVIDUAL AS THE FILE %s HAS NOT BEEN FOUND!!!\n', fname);
end
end
