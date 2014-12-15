function conf = molMot3DConf(useFirstOrder, seedVersion, pophVersion, fitnessVersion)

if (~exist('useFirstOrder', 'var')) || isempty(useFirstOrder)
  useFirstOrder = false;
end
if (~exist('seedVersion', 'var')) || isempty(seedVersion)
  seedVersion = 1;
end
if (~exist('pophVersion', 'var')) || isempty(pophVersion)
  pophVersion = 1;
end
if (~exist('fitnessVersion', 'var')) || isempty(fitnessVersion)
  fitnessVersion = 1;
end

conf.params.genome2str          = @convertGenomeStr;
conf.params.evaluateGenome      = @evaluarMolMot;
conf.params.multiEvFun          = @evaluarMolMots3D;
conf.params.writeMejorFun       = @DAwriteMejorFun;
%rows:
%       1: variable/field name
%       2: default value
%       3: write format in poblacion.txt (empty string means omission)
%       4: order in poblacion.txt
%       5: transform function for writing in poblacion.txt
%       6: read format from poblacion.txt
conf.params.statNamesFull           = {'genome', 'fitness', 'rndSeed', 'devTimeSpent', 'evTimeSpent', 'cpuTime', 'npoints', 'nsprings1', 'atpChg1', 'atpChg2', 'offsetAbs', 'offsetRel', 'timeAttached', 'spectralGap'  ;
                                       {[]},     nan,       nan,       nan,            nan,           nan,       nan,       nan,         nan,       nan,       nan,         nan,         nan,            nan            ; ...
                                       '',       '%17s',    '%s',      '%-06.3f',      '%-06.3f',     '%-06.3f', '%03d',    '%03d',      '%02d',    '%02d',    '%13d',      '%13d',      '%7d',          '%7d'          ; ...
                                       nan,      01,        06,        03,             04,            05,        11,        12,          17,        18,        19,          20,          21,             22             ; ...
                                  @genoma2Char3, @mat2str,  @num2hex,  '',             '',            '',        '',        '',          '',        '',        '',          '',          '',             ''             ; ...
                                       '%s',     '%f',      '%s',      '%f',           '%f',          '%f',      '%f',      '%f',        '%f',      '%f',      '%f',        '%f',        '%f',           '%f'          };                                                                                      
conf.params                                         = buildFieldStats(conf.params);
conf.params.buildFieldStatsFun                      = @buildFieldStats;
conf.params.writeGenealogyFun                       = @DAwriteGenealogyFun;

conf.params.genome.fold.numPoints                   = 50;
conf.params.genome.fold.lmin                        = 3.8;
conf.params.genome.fold.l0                          = 10;
conf.params.genome.fold.maxRecursion                = 20;
conf.params.genome.fold.maxReps                     = 200;
conf.params.genome.fold.nlast                       = 10;
conf.params.genome.fold.nconns                      = 4;
conf.params.genome.fold.considerAll                 = false;
conf.params.genome.fold.returnEmptyInGA             = true;
conf.params.genome.poisson                          = {0.5, 0:10, poisspdf(0:10, 0.5)};
conf.params.genome.mutationFactorRange              = [0.5 1.5];

conf.params.mutateGenome                            = @mutateFold;
conf.params.makeupParams                            = [];
conf.params.generarGenomas                          = @createRandomFolds;

conf.params.selection.mode                          = 'FUSS'; %'FUSS'; 'PACO';
conf.params.selection.FUSS.epsMode                  = 'rel'; %'abs'; 'rel';
conf.params.selection.FUSS.epsDist                  = 0.05;
conf.params.selection.FUSS.mode                     = 'levelled1'; %'levelled1'; 'plain';
conf.params.selection.FUSS.ntimes                   = 2; 1; 3;
conf.params.rewindNans                              = true;%true; false;
conf.params.nansReplaced                            = true;
conf.params.newsQuotaMode                           = 'plain'; %'none';  'plain'; 'rand';      'custom';
conf.params.newsQuota                               = 0.05; %NOTHING; 0.3;     [0.25 0.35]; @(generacion, rangeid, generationInRange, range, P, params)ceil(rand*10);

conf.params.GAFun                                   = @GAGeneric;%@GAGenericDecoupledDev;
conf.params.writeRangeID                            = true;
conf.params.writeGeneration                         = @writeGen;

conf.params.evparams.verbose                        = false;
conf.params.evparams.numTicksForStabilization       = 2;
conf.params.evparams.methodToStop                   = 0; %0, -1, +1
conf.params.evparams.useFirstOrder                  = useFirstOrder;
conf.params.evparams.doCheckVecs                    = false; %too time consuming unless absolutely necessary
conf.params.evparams.dev.fixedTimeSpent             = 500;
conf.params.evparams.dev.maxCPUTimeEval             = 100;
conf.params.evparams.dev.devEThreshold              = 0.01;
conf.params.evparams.walker.fixedTimeSpent          = 1000;
conf.params.evparams.walker.maxCPUTimeEval          = 1000;
conf.params.evparams.walker.devEThreshold           = -inf;
conf.params.evparams.relax.fixedTimeSpent           = 100;
conf.params.evparams.relax.maxCPUTimeEval           = 100;
conf.params.evparams.relax.devEThreshold            = 0.01;

conf.params.jiggleSelection                         = eps(1);
conf.params.evparams.record.positions               = false;
conf.params.evparams.record.useWalks                = 'longest'; %'aggregate';
conf.params.evparams.record.toes                    = 'light1';
conf.params.evparams.record.nrecToes                = 20;
conf.params.evparams.fitnessCalc                    = 'hoh';%'mode3'; 'mode2';%'plain';%'mode1'; %'hoh';
conf.params.evparams.fitness.HOH.mode               = 'abs'; %'absabs'; %'abs'; %'norm'; %'rel1';
conf.params.evparams.fitness.HOH.zeroIsNan          = true; %true; false;
conf.params.evparams.spectralGap.rotCutoff          = 1e-12;
conf.params.evparams.spectralGap.gapToRecord        = 3;
conf.params.evparams.relax.doRelax                  = true;
conf.params.evparams.relax.holdOnRow                = true;
conf.params.evparams.walker.fixUnixWindows.doFix    = true;
conf.params.evparams.walker.fixUnixWindows.epsilon  = eps(1e1);
conf.params.evparams.walker.fixUnixWindows.maxilon  = eps(1e10);
conf.params.evparams.walker.minAmountOfToes         = [];%minAmountOfToesFun;%@(conns, vs) min(8, ceil(numel(vs)/2)); %[];
conf.params.evparams.walker.cutoff                  = 1.5;
conf.params.evparams.walker.brownianFactor          = 0;1; %-0.1; %for the moment, do not use brownian motion, as we are not able of stripping small structures from large fitness from brownian excursions instead of good processivity
conf.params.evparams.walker.tooChangedScale         = 0.2;
conf.params.evparams.walker.rewiringPlace           = 'none'; %'none'; 'before'; 'after';
% conf.params.evparams.walker.rewiringMode            = {'geotopological', 'georecurtopo'}; %'geotopological'; %'topological'
% conf.params.evparams.walker.rewireLength            = [atomLength*2      atomLength/2];   %atomLength*2;     %nan;
conf.params.evparams.walker.ABCMode                 = 'heurconcave'; %'meanpos'; 'heurconcave';
conf.params.evparams.walker.tooSmallNeighbour       = 2/3;
conf.params.evparams.walker.springK                 = 100;
conf.params.evparams.walker.pointRad                = (conf.params.genome.fold.lmin-0.6)/2;
conf.params.evparams.walker.pointAllr               = conf.params.genome.fold.lmin/2;
conf.params.evparams.walker.row.ballRad             = (conf.params.genome.fold.lmin-0.6)/2;
conf.params.evparams.walker.row.ballAllr            = conf.params.genome.fold.lmin/2;
conf.params.evparams.walker.row.ballSep             = conf.params.genome.fold.lmin;
conf.params.evparams.walker.initialDepth            = 1;
conf.params.evparams.walker.row.ballPenk            = conf.params.evparams.walker.springK*2;
conf.params.evparams.walker.row.typeSpec            = [2 2];
conf.params.evparams.walker.row.active              = true;
conf.params.evparams.walker.row.allPoints           = true;
conf.params.evparams.walker.row.lockFreeze          = true;
conf.params.evparams.walker.row.attachMode          = 0; %0; 1; ¿2?; ¿3?;
conf.params.evparams.walker.row.spchg.doit          = true;
conf.params.evparams.walker.row.spchg.constT        = 1;
%conf.params.evparams.walker.row.freezeK             = 1/conf.params.evparams.walker.row.ballPenk;
conf.params.evparams.walker.leg.dynLengths          = false;
conf.params.evparams.walker.leg.dynLengthFactor     = 0.5;
conf.params.evparams.walker.leg.atpK                = conf.params.evparams.walker.springK;
conf.params.evparams.walker.leg.timeAtp             = 50;
conf.params.evparams.walker.leg.timeRelax           = 20;
conf.params.evparams.walker.leg.stateInit           = {'sticky', 'relaxation'};
conf.params.evparams.walker.leg.stateInitTime       = [nan,      0];
conf.params.evparams.walker.correct.rotation        = true;
conf.params.evparams.walker.correct.sign            = true;
conf.params.evparams.walker.correct.pairIndex       = true;
conf.params.evparams.walker.correct.singleRow       = true;
conf.params.evparams.walker.correct.rewiring        = true;
conf.params.evparams.walker.correct.attach          = true;
conf.params.evparams.walker.correct.rowRad          = true;
conf.params.evparams.walker.useHull                 = false;
conf.params.evparams.walker.extendInterface         = true;
conf.params.evparams.walker.eigCalc                 = 'full'; %'full'; 'sparse';
conf.params.evparams.genome.quant.nbits             = 16;
conf.params.evparams.genome.quant.nfrac             = 9;
conf.params.evparams.genome.quant.q                 = [];%quantizer('ufixed', [conf.params.evparams.genome.quant.nbits conf.params.evparams.genome.quant.nfrac]);
conf.params.evparams.genome.L0                      = conf.params.genome.fold.l0;
conf.params.evparams.genome.Kmut                    = conf.params.evparams.walker.springK;
conf.params.evparams.genome.Kmut_chain              = conf.params.evparams.walker.row.ballPenk;
conf.params.evparams.walker.d3.borderRating         = 'handoverhand1'; %'handoverhand1'; 'method1'
conf.params.evparams.walker.d3.distRatio            = 1/conf.params.genome.fold.lmin;
conf.params.evparams.walker.d3.borderMode           = 'dist2Near'; %'notconnected'; 'anyone'; 'onlyconnected'; 'dist2'; 'dist2Near';
conf.params.evparams.walker.d3.borderDist           = conf.params.genome.fold.lmin*3;
conf.params.evparams.walker.d3.filterATPTooClose    = conf.params.genome.fold.lmin*0.9;
conf.params.evparams.walker.d3.ATPConnMode          = 'simplest'; %'neighsPassive'; %'simplest';
conf.params.evparams.walker.d3.rotationMode         = 'handoverhand';%'inchworm'; %'handoverhand';
conf.params.evparams.walker.d3.hoh.pivot            = @pivotFunMolMot3DConf;%conf.params.genome.fold.numPoints; %select pivoting point (default should be the last one)
conf.params.evparams.walker.d3.hoh.orientPoint      = @orientFunMolMot3DConf; %select orientating point (default should be the second one)
conf.params.evparams.walker.d3.hoh.pinpoint         = 1; %select pinpoint point (default should be the first one)
conf.params.evparams.walker.d3.hoh.swingAngle       = pi/12; %15º==pi/12, 10º==pi/18, 5º==pi/36
conf.params.evparams.walker.d3.hoh.fusionMode       = 'justOnePoint';
conf.params.evparams.walker.d3.hoh.duplicateBy      = 'rotation'; %'mirroring'; %'rotation';
conf.params.evparams.genome.minConns                = 3;
conf.params.evparams.genome.maxChainDist            = conf.params.genome.fold.lmin*1.1;
conf.params.evparams.genome.minChainDist            = conf.params.genome.fold.lmin*0.9;
conf.params.evparams.genome.gap.useIt               = false;
conf.params.evparams.genome.gap.toRecord            = conf.params.evparams.spectralGap.gapToRecord;
conf.params.evparams.genome.gap.interval            = [];
if conf.params.genome.fold.returnEmptyInGA
  conf.params.evparams.genome.fold                  = conf.params.genome.fold;
  conf.params.evparams.genome.fold.generarGenomas   = conf.params.generarGenomas;
end

if strcmp(conf.params.evparams.fitnessCalc, 'hoh')
  conf.params.statNamesFull           = [conf.params.statNamesFull, ...
                                       {'num_adv', 'sum_adv', ; ...
                                        nan,       nan,       ; ...
                                        '%7d',     '%7d',     ; ...
                                        50,        51,        ; ...
                                        '',        '',        ; ...
                                        '%f',      '%f'      }];
  conf.params                         = buildFieldStats(conf.params);
end


if conf.params.evparams.useFirstOrder
  conf.params.evparams.doCheckVecs = true; %it is necessary
  k1 = conf.params.evparams.walker.springK/5;
  k2 = conf.params.evparams.walker.row.ballPenk/5;
  conf.params.evparams.walker.springK = k1;
  conf.params.evparams.walker.row.ballPenk = k2;
  conf.params.evparams.walker.leg.atpK = k1;
  conf.params.evparams.ss.u = 0;
  conf.params.evparams.genome.Kmut = k1;
  conf.params.evparams.genome.Kmut_chain = k2;
end

conf.fun.developGenome          = @DADevelopGenome;
conf.fun.gnmEquals              = @genomeEquals;

conf.fun.loadSim                = makeLoadFun(@loadFun1, conf.params);
conf.fun.phylotree              = phyloFun(@getBasicInf1);
conf.fun.getInfo                = @getBasicInf1;
conf.fun.printSS                = @printLastFrame;
conf.fun.replayGenome           = @replayGenome;
conf.fun.g2s                    = @genoma2Char3;
conf.fun.s2g                    = @char3ToGenoma;

conf.ssParams.tickWidth         = 0.5;
conf.ssParams.stepsByTick       = 50;
conf.ssParams.ndims             = 3;
conf.ssParams.U                 = 1;5;
if conf.params.evparams.useFirstOrder
conf.ssParams.U                 = 0;
end
conf.ssParams.rad               = conf.params.evparams.walker.row.ballRad;
conf.ssParams.allr              = conf.params.evparams.walker.row.ballAllr;
conf.ssParams.penk              = conf.params.evparams.walker.row.ballPenk;
conf.ssParams.K                 = conf.params.evparams.walker.springK;
conf.ssParams.stretchedFactor   = 0.75;
conf.ssParams.compressedFactor  = 1.25;
conf.samples.ss                 = BasicSystem(conf.ssParams.ndims, 0, conf.ssParams.U, conf.ssParams.tickWidth, conf.ssParams.stepsByTick);
conf.samples.ss.pointDefaultVals(4:6) = {conf.ssParams.rad, conf.ssParams.allr, conf.ssParams.penk};

conf.samples.params             = conf.params;
conf.samples.params.evparams.ss = conf.samples.ss;

conf.samples.postparams.saveRnd                   = true;
conf.samples.postparams.PRESIONES                 = 0.01;       %[1 0.2 0.1 0.0001 0.5];% [0.75 0.5 1/3 0.2 0.1 0.0001];
conf.samples.postparams.ensayos                   = 1; %nº de veces que se hace el experimento para cada presion
conf.samples.postparams.generaciones              = 200;%1000;%500;         % número de intentos para encontrar un árbol de una relación de aspecto dada
conf.samples.postparams.tam_poblacion             = 200;45;100;45;%90;%1000;%300;         % número de individuos en la población
conf.samples.postparams.tooLongEvTime             = inf;
conf.samples.postparams.plotting                  = true; 
conf.samples.postparams.getRasters                = false; 
conf.samples.postparams.evaluateAlways            = false;true;
conf.samples.postparams.logInEval                 = false;
conf.samples.postparams.logInMutation             = false;
conf.samples.postparams.terclus.concurrentWorkers = conf.samples.postparams.tam_poblacion;40;45;%48;%2; %number of concurrent workers in a job
conf.samples.postparams.elitism                   = false;

conf.samples.mainfun                              = DAmainFun(conf.samples.params, conf.samples.postparams);
end

function pivot  =  pivotFun(ss); pivot  = size(ss.pos,1);   end
function orient = orientFun(ss); orient = size(ss.pos,1)-1; end
function fun      = nanQuota075
  fun = @(gen, rid, genInRange, range, P, params) ((sum(isnan(P.fitness))==numel(P.fitness)) * numel(P.fitness)) + ((sum(isnan(P.fitness))~=numel(P.fitness)) * round(sum(isnan(P.fitness))*0.75));
end

function equ = genomeEquals(g1, g2)
if iscell(g1)
  p1 = g1{1};
  m1 = g1{2};
else
  p1 = g1;
  m1 = [];
end
if iscell(g2)
  p2 = g2{1};
  m2 = g2{2};
else
  p2 = g2;
  m2 = [];
end
mn1 = min(p1);
mn2 = min(p2);
p1(:,1) = p1(:,1)-mn1(1);
p1(:,2) = p1(:,2)-mn1(2);
p1(:,3) = p1(:,3)-mn1(3);
p2(:,1) = p2(:,1)-mn2(1);
p2(:,2) = p2(:,2)-mn2(2);
p2(:,3) = p2(:,3)-mn2(3);
equ = isequalwithequalnans(p1, p2) && isequalwithequalnans(m1, m2);
end

function output = genoma2Char3(genoma)
if iscell(genoma)
  g = genoma{1};
else
  g = genoma;
end
tabla = ['A':'Z', 'a':'z', '0':'9', '!?'];
g = g(:);
gb = zeros(numel(g), 16, 'uint16');
for k=1:16
  gb(:,k) = bitget(g, k);
end
vals = uint16(realpow(2, 0:5));
n1 = sum(bsxfun(@times, gb(:,1:6),   vals), 2);
n2 = sum(bsxfun(@times, gb(:,7:12),  vals), 2);
n3 = sum(bsxfun(@times, gb(:,13:16), vals(1:4)), 2);
g  = char(reshape(tabla(1+[n1 n2 n3]'), 1, []));
if iscell(genoma)
  genoma{1} = g;
else
  genoma    = g;
end
output = any2str(genoma, false, true);
end

function output = char3ToGenoma(str)
output = eval(str);
if iscell(output)
  g = output{1};
else
  g = output;
end
if ischar(g)
  tabla = ['A':'Z', 'a':'z', '0':'9', '!?'];
  tablainv = repmat(' ', 1, 255);
  tablainv(tabla) = 1:64;
  vg = uint16(reshape(tablainv(g), 3, []))'-1;
  nums = vg(:,1)+vg(:,2)*64+vg(:,3)*4096;
  g = reshape(nums, [], 3);
end
if iscell(output)
  output{1} = g;
else
  output    = g;
end
end

function writeGen(params, G, generation, rangeid, indexNewG, rangeStartIdx)
%indexes = (1:numel(G.genome))';
indexes = rangeStartIdx-1+indexNewG;
%genomes = G.genome;
genomes = G.genome(indexNewG);
% %rasters = cellfun(@(x) convertCompressedPos(x.pos, params.evparams), G.raster, 'uniformoutput', false);
% %rasters = cellfun(@(x) convertCompressedPos(x.pos, params.evparams), G.raster(indexNewG), 'uniformoutput', false);
% rasters = G.raster(indexNewG);
% evparams = params.evparams;
% for k=1:numel(rasters)
%   [g snapped err] = convertCompressedPos(rasters{k}.pos, evparams, false);
%   if err
%     rasters{k}    = rasters{k}.pos;
%   else
%     rasters{k}    = g;
%   end
% end
%res = struct('Generation', {generation}, 'rangeid', {rangeid}, 'indexes', {indexes}, 'indexNewG', {indexNewG}, 'genomes', {genomes}, 'rasters', {rasters}); %#ok<NASGU>
res = struct('Generation', {generation}, 'rangeid', {rangeid}, 'indexes', {indexes}, 'indexNewG', {indexNewG}, 'genomes', {genomes}); %#ok<NASGU>
save([params.nomdir sprintf('lastFrames%04d_%d.mat', generation, rangeid)], 'res');
end

function fun = minAmountOfToesFun
  fun = @(conns, vs) min(8, ceil(numel(vs)/2));
end

function fun = DAmainFun(params, postparams)
fun = @(varargin)mainGeneric(params, postparams, varargin{:});
end

function output = convertGenomeStr(input)
if ischar(input)
  output = eval(input);
else
  output = any2str(input, false, true);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params                                = buildFieldStats(params)
params.indivStructFull                         = [{'raster', ; ...
                                                  {[]},      ; ...
                                                  '',        ; ...
                                                  nan,       ; ...
                                                  '',        ; ...
                                                  '',        }, ...
                                                  params.statNamesFull];
params.evparams.statNames                      = params.statNamesFull(1:2,:);
params.fieldsToFill                            = 1:size(params.indivStructFull,2); %DO NOT skip first field 'genome' when reaping results
params.indivStruct                             = params.indivStructFull(1:2,:);
params.writeIndividualFun                      = makeWriteFun(params);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = DAwriteMejorFun(fres, switchGenomeStr)
  fun = @(mejorEnConserva) ...
    fprintf(fres,'%g %g %g %s %s %s\n',mejorEnConserva.gen, mejorEnConserva.rind, mejorEnConserva.num, num2hex(mejorEnConserva.data.rndSeed), mat2str(mejorEnConserva.data.fitness), switchGenomeStr(mejorEnConserva.data.genome{1})); % actualizamos el fichero de resultados
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = DAwriteGenealogyFun(farb)
  fun = @(generationNew, rangeindNew, indexNew, changeStr, generationOld, rangeindOld, indexOld) ...
    fprintf(farb,'%d %d %d %s %d %d %d\n',generationNew,rangeindNew,indexNew,changeStr,generationOld, rangeindOld, indexOld);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = makeLoadFun(loadFun, params)
fun = @(varargin) loadFun(params, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function poph = loadFun1(confparams, varargin)
basedir       = varargin{1};
fname         = [basedir    filesep '..' filesep 'estado.mat'];
if exist(fname, 'file')
  load(fname, 'params');
  if isfield(params, 'indivStructFull')
    indivStruct = params.indivStructFull; %get full info on fields
  end
end
if not(exist('indivStruct', 'var'))
  indivStruct = confparams.indivStructFull; %get full info on fields
end
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
basicinf = sprintf('Index in struct: %04d, Generation: %03d, individual: %03d, np: %d, ns: %d,\nFITNESS: %08.6f, chg1: %d, chg2: %d, spectralGap: %s,\noffsetA: %d, offsetR: %d devTimeSpent: %04d, evTimeSpent: %04d, cpuTimeSpent: %04d', ...
                   idx, ...
                   poph.generation(idx), ...
                   poph.individual(idx), ...
                   poph.npoints(idx), ...
                   poph.nsprings1(idx), ...
                   poph.fitness(idx), ...
                   poph.atpChg1(idx), ...
                   poph.atpChg2(idx), ...
                   mat2str(poph.spectralGap(idx)), ...
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

function structs = replayGenome(poph, mode, index, varargin)
rid = 1;
switch mode
  case 'poph'
    idp = index;
    gen = poph.generation(index);
    idx = poph.individual(index);
  case 'gen'
    gen = index(1);
    idx = index(2);
    idp = poph.triplet2idx(gen, rid, idx);
end
if ~poph.params.evaluateAlways
  idp = getFirstUnchangedAncestor(poph, idp);
  gen = poph.generation(idp);
  idx = poph.individual(idp);
end
fname = [poph.basedir filesep sprintf('lastFrames%04d_%d.mat', gen, rid)];
if exist(fname, 'file')
  load(fname, 'res');
  genome = res.genomes{res.indexes==idx};
  structs = evaluarMolMot(genome, poph.params.evparams, poph.rndSeed{idp}, [], varargin{:});
end
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
margen = 1.1;
pltargs = {'axisEqual', false, 'axisSquare', true, 'circleFaces', 10, 'showBalls', 'allFiltered', 'showRowContacts', true};
%pltargs = {'axisEqual', false, 'axisSquare', true, 'circleFaces', 10, 'showBalls', 'none', 'showRowContacts', false};
if exist(fname, 'file')
  load(fname, 'res');
  if exist('res', 'var')
    initial = res.genomes{res.indexes==iid};
    stDev   = showGenome(initial, poph.params);
    stDev.genome = initial;
%     if ~isfloat(initial)
%       initial = convertCompressedPos(initial, poph.params.evparams);
%     end
%     [sp1 sp2] = find(triu(distanceMatrixSquared(initial)<=(realpow(poph.params.evparams.genome.L0, 2)), 1));
%     r         = realsqrt(sum(realpow(initial(sp1,:)-initial(sp2,:), 2), 2));
%     stDev     = addSprings(addPoints(poph.params.evparams.ss, initial), [sp1 sp2], poph.params.evparams.walker.springK, r);
%     final = res.rasters{res.indexes==iid};
%     if ~isfloat(final)
%       final = convertCompressedPos(final, poph.params.evparams);
%     end
%     stEv      = addSprings(addPoints(poph.params.evparams.ss, final), [sp1 sp2], poph.params.evparams.walker.springK, r);
    if all(isfield(stDev, {'pos', 'springEnds'}))  %&& all(isfield(stEv, {'pos', 'springEnds'}))
      stEv.row = poph.params.evparams.walker.row;
      asp = BallPlotter('axisWindow', getAxisWindow(stDev.pos, margen), pltargs{:});
      if isempty(getBasicInfFun)
        drawSS(asp, stDev);
      else
        drawSSAndText(asp, stDev, getBasicInfFun(poph, indRequested));
      end
      set(gcf, 'name', sprintf('GENOME (requested g:%d i:%d; displayed g:%d i:%d', gir, iir, gid, iid));
      assignin('base', 'dss', stDev);
      %maximize(gcf);
%       asp2 = BallPlotter('axisWindow', getAxisWindow(stEv.pos, margen), pltargs{:});
%       drawSSAndText(asp2, stEv, getBasicInfFun(poph, indRequested));
%       set(gcf, 'name', sprintf('AFTER'));
%       %maximize(gcf);
%       assignin('base', 'zss', stEv);
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
