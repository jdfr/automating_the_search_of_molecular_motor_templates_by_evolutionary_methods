function structs = evaluarMolMot(genome, evparams, rndSeedDev, rndSeedEv, interaction, outputrec, recFreq)

initcpu = cputime;

if ~exist('recFreq', 'var')
  recFreq = [];
end
if numel(recFreq)==1
  recFreq = [recFreq, recFreq];
end

if (numel(outputrec)==1)
  outputrec = [outputrec outputrec];
end

if numel(interaction)<2
  interaction = [interaction interaction];
end
beVerbose     = interaction(2);

if outputrec(1)
  rec = makeNewRecorder(evparams.ss, 'MemRecorder', 10000, 10000, 100, 100);
else
  rec = [];
end

evparams = transformCorrectParams(evparams);

if      isfield(evparams, 'doCheckVecs')
  evparams.ss.checkVecs = evparams.doCheckVecs;
elseif ~isfield(evparams.ss, 'checkVecs')
  evparams.ss.checkVecs = false;
end
if      isfield(evparams, 'useFirstOrder')
  evparams.ss.firstOrder = evparams.useFirstOrder;
elseif ~isfield(evparams.ss, 'firstOrder')
  evparams.ss.firstOrder = false;
end
%first order integration requires a careful choice of spring stiffnesses.
%This command is a useful shortcut to consistently change all stiffnesses 
%for 3D simulations:
%   k1=20; k2=40; poph.params.evparams.walker.springK = k1; poph.params.evparams.walker.row.ballPenk = k2; poph.params.evparams.walker.leg.atpK = k1; poph.params.evparams.ss.u = 0; poph.params.evparams.genome.Kmut = k1; poph.params.evparams.genome.Kmut_chain = k2;

ssCanonical = evparams.ss;
goon = true;
spectrum = [];
if ~isempty(recFreq)
  evparams.recFreq = recFreq(1);
end
if size(ssCanonical.pos,2)<3
  [ss rec] = manipulateGenome2D(genome, rndSeedDev, evparams, interaction, beVerbose, rec);
else
  uniqueSeed = rndSeedDev;
  if ~isempty(uniqueSeed)
    sowSeed(uniqueSeed);
    rndSeedDev = rand*1e5;
    rndSeedEv = rand*1e5;
  end
  [ss rec goon genome spectrum] = manipulateGenome3D(genome, rndSeedDev, evparams, interaction, beVerbose, rec);
  rndSeedDev = uniqueSeed;
end

devTimeSpent = ss.t-evparams.ss.t;

structs = struct('ssDev', [], 'ssLegs', [], 'ssAfter', [], 'ssRelaxed', [], 'stats', []);
structs.ssDev = ss;
if outputrec(1)
  structs.ssDev.rec = rec; clear rec;
end

goon = goon && (~any(isfield(ss, {'simError', 'checkFailed'})));
%goon = false; structs.ssAfter = ss; structs.ssRelaxed = ss;
if goon
  sowSeed(rndSeedEv);
  if beVerbose; fprintf('Calculate legs for the structure\n'); end
  [evparams legerror] = prepareLegsSimulation(ss, evparams);
  goon = isempty(legerror);
  if (~goon) && beVerbose
    fprintf('Error calculating legs: <%s>\n', legerror);
  end
  structs.ssLegs = evparams.ss;
else
  structs.ssLegs    = ss;
end

goon = goon && not(isempty(structs.ssLegs)) && isfield(structs.ssLegs, 'legs');% && (~ischar(structs.ssLegs.legs));

if goon

  if outputrec(2)
    rec = makeNewRecorder(evparams.ss, 'MemRecorder', 10000, 10000, 100, 100);
  else
    rec = [];
  end
  
  isrecord = isfield(evparams, 'record');
  evparams.ss.record.cm = isrecord && evparams.record.positions;
  if evparams.ss.record.cm
    nrecords          = round(evparams.fixedTimeSpent/(evparams.ss.tick/evparams.ss.stepsByTick))+5;
    %evparams.ss.rc.P = zeros(nrecords,2);
    evparams.ss.rc.P  = zeros(nrecords,1);
    evparams.ss.rc.F  = false(nrecords,1);
    evparams.ss.rc.i  = 1;
  end
  evparams.ss.record.toes1 = isrecord && isfield(evparams.record, 'toes') && strcmp(evparams.record.toes, 'light1');
  if evparams.ss.record.toes1
    nrecords          = evparams.record.nrecToes;
    evparams.ss.rc.Pt = zeros(nrecords,6);
    evparams.ss.rc.Ft = zeros(nrecords,1);
    evparams.ss.rc.Tt = zeros(nrecords,1);
    evparams.ss.rc.it = 1;
  end

  if beVerbose; fprintf('Do the (hopefully) walking\n'); end
  if ~isempty(recFreq)
    evparams.recFreq = recFreq(2);
  end
  [ss rec] = DADevelopGenome([], evparams, [], interaction, rec);
  
  structs.ssAfter = ss;
  if outputrec(2)
    structs.ssAfter.rec = rec; clear rec;
  end
  
else
  
  structs.ssAfter   = ss;
  structs.ssRelaxed = ss;
  
end

evTimeSpent = ss.t-evparams.ss.t;

goon = goon && (~any(isfield(ss, {'simError', 'checkFailed'})));

if goon
  if evparams.relax.doRelax
    ss                  = relaxWalker(ss, evparams, interaction, beVerbose);
    structs.ssRelaxed   = ss;
  end
end

goon = goon && (~any(isfield(ss, {'simError', 'checkFailed'})));

if isfield(ss, 'simError')
  structs.ssError       = ss;
elseif isfield(ss, 'checkFailed')
  structs.ssFailed      = ss;
end

if isfield(structs.ssDev, 'atomPoints')
  atomPoints = structs.ssDev.atomPoints;
else
  atomPoints = [];
end

if beVerbose; fprintf('Evaluate the structure\n'); end
structs.stats           = evaluateFitness(genome, spectrum, ssCanonical, structs.ssDev, structs.ssLegs, ss, evparams, ...
                                          goon, rndSeedDev, rndSeedEv, ...
                                          initcpu, devTimeSpent, evTimeSpent, atomPoints);

%structs.stats.fitness = rand;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = relaxWalker(ss, evparams, interaction, beVerbose)
ss.brownianFractor           = 0;
if ~evparams.relax.holdOnRow
  ss.chgStep                 = [];
end
ss.postProcess               = [];
for k=1:numel(ss.legs)
  k1                         = mod(k,numel(ss.legs))+1;
  if ~isempty(ss.legs(k).atpb)
    [ss indexesAfterRemove]  = removeSprings(ss, ss.legs(k).atpb);
    if ~isempty(ss.legs(k1).atpb)
      ss.legs(k1).atpb       = indexesAfterRemove(ss.legs(k1).atpb);
    end
    ss.legs(k).atpb = [];
  end
  if ~isempty(ss.legs(k).atp)
    [ss indexesAfterRemove]  = removePoints(ss, ss.legs(k).atp);
    if ~isempty(ss.legs(k1).atp)
      ss.legs(k1).atp        = indexesAfterRemove(ss.legs(k1).atp);
    end
    ss.legs(k).atp = [];
  end
end
ss                         = prepareSimulationStructure(ss);
evparams.fixedTimeSpent    = evparams.relax.fixedTimeSpent;
evparams.maxCPUTimeEval    = evparams.relax.maxCPUTimeEval;
evparams.devEThreshold     = evparams.relax.devEThreshold;

evparams.ss                = ss;

if beVerbose; fprintf('Allow the structure to relax after (hopefully) walking\n'); end

ss = DADevelopGenome([], evparams, [], interaction);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = evaluateFitness(genome, spectrum, ssCanonical, ssDev, ssLegs, ssFin, evparams, goon, rndSeedDev, rndSeedEv, initcpu, devTimeSpent, evTimeSpent, atomPoints)

aux   = [evparams.statNames(1,:); array2cell(evparams.statNames(2,:))];
stats = struct(aux{:});

useCells = isfield(stats, 'ncells');

if ~useCells
  stats.genome = {genome};
end

if goon
  if useCells
    ssFin.atomPoints        = atomPoints;
    stats.overload2         = calculateOverload(evparams, ssFin);
  end

  stats.atpChg1             = ssFin.legs(1).atpCount;
  stats.atpChg2             = ssFin.legs(2).atpCount;

  modelScales               = [max(reshape(distanceMatrix(ssLegs.pos), [], 1)); ...
                               max(reshape(distanceMatrix(ssFin.pos),  [], 1))];
  tooChanged                = isfield(evparams.walker, 'tooChangedScale') && ( (min(modelScales)/max(modelScales)) < evparams.walker.tooChangedScale );

  w1                        = ssLegs.m/sum(ssLegs.m);
  cm1                       = sum(bsxfun(@times, ssLegs.pos, w1));
  w2                        = ssFin.m(1:numel(ssLegs.m))/sum(ssFin.m(1:numel(ssLegs.m)));
  cm2                       = sum(bsxfun(@times, ssFin.pos(1:numel(ssLegs.m),:), w2));

  displacement              = cm2-cm1;

  if isfield(evparams, 'record') && evparams.record.positions
    contact          = ssFin.rc.F(1:ssFin.rc.i-1);
    cmposx           = ssFin.rc.P(1:ssFin.rc.i-1,1);
    events           = diff(contact);
    evpos            = find(events);

    %the states are the events, that is to say: 
    %    - if transition from 0 to 1, states(k)==+1
    %    - if transition from 1 to 0, states(k)==-1
    %    - the first state is prepended, and determined as (+1) if 
    %      the history series starts in 1, or (-1) if it starts in 0.
    states          = [contact(1)*2-1; events(evpos)];
    %positions of the center of mass of the structure at each state.
    %Note that we add 1 to evpos, since we measure each transition
    %[k,k+1] at the point k+1. Also note that we add the last position.
    %This has the following effect:
    %  if evpos(end)+1 == numel(cmposx), the difference will be 0 and 
    %  it will not add up
    %  if evpos(end)+1  < numel(cmposx), the difference will add up
    %  only if states(end)==1
    cmposst         = [cmposx(1); cmposx(evpos+1); cmposx(end)];
    stepSim         = ssCanonical.tick/ssCanonical.stepsByTick;
    cmtst           = [0; evpos; ssFin.rc.i-1];
    %differences between adyacent states' positions
    xdiff           = diff(cmposst);
    tdiff           = diff(cmtst)*stepSim;
    switch evparams.record.useWalks
      case 'aggregate'
        %add the differences in pairs of states [1,-1], discard the [-1, 1]
        %ones. Also, if the last state is 1 (meaning transition 
        %from 0 to 1), add it up.
        stats.offsetAbs    = abs( sum(xdiff(states>0)) );
        stats.timeAttached = abs( sum(tdiff(states>0)) );
      case 'longest'
        stats.offsetAbs    = max(abs( xdiff(states>0) ) );
        stats.timeAttached = max(abs( sum(tdiff(states>0)) ) );
        if isempty(stats.offsetAbs)
          stats.offsetAbs = 0;
          stats.timeAttached = 0;
        end
      otherwise
       error('evparams.record.useWalks %s not understood!!!!', any2str(evparams.record.useWalks));
    end
  else
    stats.offsetAbs = abs(displacement(1));
  end
  stats.offsetRel   = stats.offsetAbs/max(modelScales);
  
  if isfield(ssFin.record, 'toes1') && ssFin.record.toes1
    strides.ts          = ssFin.rc.Tt(1:(ssFin.rc.it-1));
    strides.legs        = ssFin.rc.Ft(1:(ssFin.rc.it-1));
    strides.oLegs       = 3-strides.legs; % if 1=>2, if 2=>1
    strides.pos{1}      = ssFin.rc.Pt(1:(ssFin.rc.it-1),1:3);
    strides.pos{2}      = ssFin.rc.Pt(1:(ssFin.rc.it-1),4:6);
    strides.lengths     = strides.pos{1}(:,1)-strides.pos{2}(:,1);
    %this means that the legs properly alternate in adhesion to the
    %filament
    strides.alternated  = strides.legs(1:end-1)~=strides.legs(2:end);
    %strides in a hand-over-hand fashion will have alterning signs
    strides.isStride    = sign(strides.lengths(1:end-1))~=sign(strides.lengths(2:end));
    stats.strides       = strides;
  end

  if ~tooChanged
    switch evparams.fitnessCalc
      case 'mode1'
        stats.fitness = realpow(stats.offsetRel, 2)*(stats.atpChg1+stats.atpChg2);
      case 'mode2'
        stats.fitness = stats.offsetRel;
      case 'mode3'
        stats.fitness = stats.offsetAbs + (stats.atpChg1+stats.atpChg2);
      case 'plain'
        stats.fitness = stats.offsetAbs;
      case 'hoh'
        areok           = find(strides.alternated & strides.isStride)+1;
        advances        = abs(strides.lengths(areok));
        stats.num_adv   = numel(advances);
        stats.sum_adv   = sum(advances);
        if evparams.fitness.HOH.zeroIsNan && (stats.num_adv==0)
          stats.fitness = nan;
        else
          switch evparams.fitness.HOH.mode
            case 'abs'
             %the walked distance in correct strides
              stats.fitness = stats.sum_adv;
            case 'absabs'
              stats.fitness = stats.sum_adv + stats.offsetAbs;
            case {'norm', 'rel1'}
              %the walked distance in correct strides, normalized by the max
              %stride
              advances = advances/max(advances);
              stats.fitness = sum(advances);
              if strcmp(evparams.fitnessHOHMode, 'rel1')
                %the walked distance in correct strides, normalized by the max
                %stride, averaged over the number of performed strides
                stats.fitness = stats.fitness/stats.num_adv;
              end
            otherwise
              error('evparams.fitnessHOHMode=%s, not understood!!!\n', any2str(evparams.fitnessHOHMode));
          end
        end
      otherwise
        error('evparams.fitnessCalc %s not understood!!!!', any2str(evparams.fitnessCalc));
    end
  end

  if isfield(evparams, 'spectralGap')
    rotCutoff   = evparams.spectralGap.rotCutoff;
    gapToRecord = evparams.spectralGap.gapToRecord;
  else
    rotCutoff   = 1e-12;
    gapToRecord = 4;
  end
  if isfield(evparams.walker, 'd3') && strcmp(evparams.walker.d3.rotationMode, 'handoverhand')
    ssANM = ssDev;
  else
    ssANM = ssLegs;
  end
%   switch evparams.walker.d3.rotationMode
%     case 'inchworm'
%       ssANM = ssLegs;
%     case 'handoverhand'
%       ssANM = ssDev;
%   end
  if not(any(isnan(ssANM.pos(:)))) && not(any(isinf(ssANM.pos(:))))
    if isempty(spectrum)
      [spectrum spectrum spectrum] = ANMdecomposition(ssANM);
      spectrum          = diag(spectrum);
    end
    if not(any(isnan(spectrum)))
      bigEnough         = spectrum>rotCutoff;
      bigEnough(1:3)    = false; %these ALWAYS ARE DISPLACEMENTS/ROTATIONS
      spectrum          = spectrum(bigEnough);
      spectrum          = log10(spectrum/spectrum(1));
      stats.spectralGap = spectrum(gapToRecord+1);
    end
  end
  
end
stats.npoints      = numel(ssLegs.m);
stats.nsprings1    = numel(ssLegs.r);
if useCells
  stats.rndSeedDev   = rndSeedDev;
  stats.rndSeedEv    = rndSeedEv;
  stats.nsprings2  = numel(ssFin.r);
  stats.ncells     = size(atomPoints, 1);
  stats.ncellsD    = size(unique(sort(atomPoints, 2), 'rows'), 1);
  stats.ncellsTeo  = ssLegs.ncellsTeo;
  stats.overload1  = ssLegs.overload1;
  stats.developed  = ssLegs.developed;
else
  stats.rndSeed    = rndSeedDev;
end
stats.devTimeSpent = devTimeSpent;
stats.evTimeSpent  = evTimeSpent;
stats.cpuTime      = cputime-initcpu;



