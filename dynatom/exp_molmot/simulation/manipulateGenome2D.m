%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ss rec goon] = manipulateGenome2D(genome, rndSeedDev, evparams, interaction, beVerbose, rec)

if iscell(genome)
  sowSeed(rndSeedDev);
  evparams.fixedTimeSpent                 = evparams.dev.fixedTimeSpent;
  evparams.maxCPUTimeEval                 = evparams.dev.maxCPUTimeEval;
  evparams.devEThreshold                  = evparams.dev.devEThreshold;

  if beVerbose; fprintf('Develop the structure\n'); end
  [ss rec] = DADevelopGenome(genome, evparams, [], interaction, rec);
  ss.ncellsTeo = DAnumCellsInGenome(genome, ss.atomDefaultVals{4}{1}.geneDistribution, ss.geneTypes);
  ss.developed  = all(isnan(ss.lockedUntil)) && isempty(ss.pShifting_M) && (~any(ss.rdyn.radapt));
  if any(isfield(ss, {'simError', 'checkFailed'}))
    ss.overload1 = nan;
  else
    ss.overload1 = calculateOverload(evparams, ss);
  end
elseif isstruct(genome)
  ss = genome;
  ss = combineBasicSystems(removePoints(removeSprings(removeAtoms(evparams.ss, (1:size(evparams.ss.atomPoints, 1))', false, false), (1:numel(evparams.ss.r))'), (1:numel(evparams.ss.m))'), ss, [evparams.ss.pointVars, evparams.ss.springVars]);
  fields = {'ncellsTeo', 'overload1', 'developed'; nan, nan, true};
  for k = 1:size(fields, 1)
    fs = fields{1,k};
    if isfield(genome, fs)
      ss.(fs) = genome.(fs);
    else
      ss.(fs) = fields{2,k};
    end
  end
  ss.ncellsTeo = genome.ncellsTeo;
  ss.overload1 = genome.overload1;
  %ss.ncellsD   = genome.ncellsD;
  ss.developed = genome.developed;
  clear genome;
else
  error('genome not understood!!!!');
end

goon = true;
