%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ss rec goon genome spectrum] = manipulateGenome3D(genome, rndSeedDev, evparams, interaction, beVerbose, rec)

goon        = true;
spectrum    = [];

%extract the genome and the mutation to be applied
if isstruct(genome)
  mutation   = genome.mutation;
  genome     = genome.genome;
elseif iscell(genome)
  mutation   = genome{2};
  genome     = genome{1};
elseif ischar(genome) && strcmp(genome, 'LAZY_GENERATION')
  %we must generate the genome by ourselves
  sowSeed(rndSeedDev);
  rndSeedDev = []; %this way, the seed is not sown again
  genome     = evparams.genome.fold.generarGenomas(1, evparams, true);
  genome     = genome{1};
  mutation   = [];
else
  mutation   = [];
end

dst         = evparams.genome.L0;

if isstruct(genome)
  if all(isfield(genome, {'pos', 'vel', 'r', 'k', 'm'}))
    ss      = combineBasicSystems(removePoints(removeSprings(evparams.ss, (1:numel(evparams.ss.r))'), (1:numel(evparams.ss.m))'), genome, [evparams.ss.pointVars, evparams.ss.springVars]);
  else
    error('This struct cannot be read!!!!');
  end
elseif isnumeric(genome)
  if isfloat(genome)
    pos     = genome;
  else
    pos     = convertCompressedPos(genome, evparams);
  end
  conns     = triu(distanceMatrixSquared(pos)<=(dst*dst), 1);
  [sp1 sp2] = find(conns);
  r         = realsqrt(sum(realpow(pos(sp1,:)-pos(sp2,:), 2), 2));
  ss        = addSprings(addPoints(evparams.ss, pos), [sp1 sp2], evparams.genome.Kmut, r);
end

if ~isempty(mutation)
  spEnds       = mutation(:,[1 2]);
  factors      = mutation(:,3);
  %the mutation consists of modifying the length of some springs
  [sps sps]    = ismember(spEnds, ss.springEnds, 'rows');
  for z=1:numel(sps)
    ss.r(sps(z))    = ss.r(sps(z)).*factors(z);
  end
  evparams.fixedTimeSpent                 = evparams.dev.fixedTimeSpent;
  evparams.maxCPUTimeEval                 = evparams.dev.maxCPUTimeEval;
  evparams.devEThreshold                  = evparams.dev.devEThreshold;
  ss.useSAP                               = true;
  %the sapFilter is a matrix which can be used to prevent some connections
  %between nearby points from being made. In this case, we want to prevent
  %the points already joined by a spring from experiencing the interactions
  sapFilter                               = zeros(size(ss.pos,1), size(ss.pos,1), 'int8');
  sapFilter(sub2ind(size(sapFilter), ss.springEnds(:), [ss.springEnds(:,2); ss.springEnds(:,1)])) = -10;
  ss.stick.rad(1:end)                     = evparams.walker.pointAllr;
  ss.stick.allr(1:end)                    = evparams.walker.pointAllr;
  ss.stick.penk(1:end)                    = evparams.genome.Kmut_chain;
  ss.sap                                  = createSweepAndPrune(ss.pos, ss.stick.allr, [], false, sapFilter);
  ss.k((ss.springEnds(:,2)-ss.springEnds(:,1))==1) = evparams.genome.Kmut_chain;
  evparams.ss                             = ss;
  if beVerbose; fprintf('Calculate mutation\n'); end;
  sowSeed(rndSeedDev);
  [ss rec]     = DADevelopGenome([], evparams, [], interaction, rec);
  if any(isfield(ss, {'simError', 'checkFailed'}))
    goon         = false;
  else
    %now, let's snap it to grid, and evaluate
    [g pos err]    = convertCompressedPos(ss.pos, evparams);
    goon           = not(err);
    if goon
      genome       = g;
      conns        = triu(distanceMatrixSquared(pos)<=(dst*dst), 1);
      [sp1 sp2]    = find(conns);
      ss.pos       = pos;
      r            = realsqrt(sum(realpow(pos(sp1,:)-pos(sp2,:), 2), 2));
      ss           = addSprings(removeSprings(ss, (1:numel(ss.r))'), [sp1, sp2], zeros(size(r)), r);
      ss.useSAP    = false;
      %make sure that:
      %    -the structure is still nominally rigid, at least
      %    -the vertices in the chain are still in a certain range
      [goon spectrum]   = validChain(pos, conns, sp1, sp2, r, evparams);
    end
  end
else
  if isfield(evparams, 'onlyValidChains')
    onlyValidChains = evparams.onlyValidChains;
  else
    onlyValidChains = true;
  end
  if onlyValidChains
    %this is necessary for otherwise perfectly fit individuals, but failing
    %to be a valid chain
    [goon spectrum] = validChain(pos, conns, sp1, sp2, r, evparams);
  else
    goon            = true;
  end
end

ss.k(1:end) = evparams.walker.springK;
