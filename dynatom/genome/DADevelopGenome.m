%%%%%%%%%%%%%%%%%%%%%%%%
function [ss varargout] = DADevelopGenome(genome, evparams, rndSeed, interaction, rec)

ss = evparams.ss;
% if ischar(genome)
%   genome = evparams.convertStr(genome);
% end

if ~isempty(genome)
  ss.atomGenome{1} = genome;
end

if (nargin>=3) && (~isempty(rndSeed))
  if iscell(rndSeed) && (numel(rndSeed)==1)
    rndSeed = rndSeed{1};
  end
  if ischar(rndSeed)
    rndSeed = hex2num(rndSeed);
  end
  rand('twister', rndSeed);
end

outputrec = nargout>1;

if outputrec
  if ~exist('rec', 'var')
    rec = makeNewRecorder(ss, 'MemRecorder', 10000, 10000, 100, 100);
  end
else
  rec = [];
end

numTicksForStabilization      = evparams.numTicksForStabilization;
recordInitialStateOut         = true;
if nargin<4
  interaction                 = [false false];
elseif numel(interaction)<2
  interaction                 = [interaction interaction];
end
useKeyboard                   = interaction(1);
beVerbose                     = interaction(2);
evparams.verbose              = evparams.verbose || beVerbose;
clearOldModifications         = false; %we manage this ourselves
timeBetween                   = ss.tick*numTicksForStabilization;

interfun                      = getIntercalateFun(cputime, evparams);
ss.newcpustart                = cputime;
if isfield(evparams, 'recFreq')
  recFreq = evparams.recFreq;
else
  recFreq = [];
end

%useKeyboard = false;

ss          = prepareSimulationStructure(ss);
ss.eOld     = calculateEnergy(ss);
[ss rec recordInitialState numIteration goonKeyboard] = simulateSystemIntercalatedCoarse(ss, rec, ss.tick:ss.tick:timeBetween, interfun, recordInitialStateOut, useKeyboard, clearOldModifications, recFreq); %#ok<NASGU> %varargin = {clearOldModifications, odefun}

if (~any(isfield(ss, {'simError', 'checkFailed'}))) && ... 
   ( any(isnan(ss.pos(:))) || any(isnan(ss.vel(:))) || any(isinf(ss.pos(:))) || any(isinf(ss.vel(:))) || ...
     any(abs(ss.pos(:))>1e5) || any(abs(ss.vel(:))>1e5)  )
  ss.simError = 'Detected nonsensical values in position and/or velocities!!!!';
end

if outputrec
  varargout{1} = rec;
end

end

function fun = getIntercalateFun(cpustart, evparams)
fun = @(ss, rec, numIteration) intercalateFun(cpustart, evparams, ss, rec, numIteration);
end

function [ss rec recordInitialState abortSimulation] = intercalateFun(cpustart, evparams, ss, rec, numIteration) %#ok<INUSD>
eNew = calculateEnergy(ss);
if evparams.verbose && (mod(ss.t,10)==0)
  fprintf('Simulated time: %05.01f, elapsed CPU time: %06.01f\n', ss.t, cputime-cpustart);
end
recordInitialState = false;
if any(isfield(ss, {'simError', 'checkFailed'}))
  abortSimulation    = true;
else
  
  
  abortSimulation = ( (~isfield(ss, 'pShifting_M')) || isempty(ss.pShifting_M) && all(isnan(ss.lockedUntil)) ) && ...
                    (~any(ss.rdyn.radapt))  && (abs(ss.eOld-eNew)<evparams.devEThreshold);
  
  if ~abortSimulation
    if evparams.methodToStop>(-1)
      abortSimulation = abortSimulation || (ss.t>=evparams.fixedTimeSpent);
    end
    if evparams.methodToStop<1 
      abortSimulation = abortSimulation || ( (cputime-cpustart)>evparams.maxCPUTimeEval );
    end
  end
  
  ss.eOld = eNew;
end
if ~abortSimulation
  %[rdon] = clearOldPerturbations(ss, ss.rdyn, ss.r, ss.t);
  [ss.rdyn ss.r] = clearOldPerturbations(ss, ss.rdyn, ss.r, ss.t);
end
end
