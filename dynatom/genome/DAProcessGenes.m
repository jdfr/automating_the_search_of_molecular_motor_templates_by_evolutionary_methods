%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ss r rdyn veryChanged] = DAProcessGenes(ss, r, rdyn, time)
%do point shifting makeup
[ss modifiedPointsOrSprings] = finishPointShiftings(ss);
veryChanged = modifiedPointsOrSprings;
%find atoms whose genome has finished to express a gene
atomsToUnlock = find(ss.lockedUntil<=ss.t);
if ~isempty(atomsToUnlock)
  %find atoms whose genome is empty
  silentAtoms = cellfun(@isempty, ss.atomGenome(atomsToUnlock));
  %earmark those atoms
  ss.lockedUntil(atomsToUnlock(silentAtoms)) = nan;
  %indexes of atoms with non-empty genomes
  atomsToExpressGene = atomsToUnlock(~silentAtoms);
    if ~isempty(atomsToExpressGene);
    %the atoms with non-empty genomes must express the next gene
    [ss.executingGene(atomsToExpressGene) ss.atomGenome(atomsToExpressGene)] ...
      = cellfun(@(genome) deal(genome{1}, genome(2:end)), ...
                ss.atomGenome(atomsToExpressGene), 'UniformOutput', false);
    %for each atom expressing a new gene, do it
    k=1;
    geneTypes = ss.geneTypes;
    newAtoms = [];
    mitosisHappended = false;
            while k<=numel(atomsToExpressGene)
              atomIndex = atomsToExpressGene(k);
              endGeneTime = ss.t;
              continuar = true;
              while continuar
                gene = ss.executingGene{atomIndex};
                switch gene.type
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.switchForcement
                    %error('This gene is disabled!!!');
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.mitosisSpringShift
                    if isempty(ss.pShifting_Current)
                      ss.r = r;
                      ss.rdyn = rdyn;
                      [ss endGeneTime newAtom] = atomMitosisSpringShift(ss, atomIndex, gene);
                      mitosisHappended = true;
                      r = ss.r;
                      rdyn = ss.rdyn;
                      modifiedPointsOrSprings = true;
                      ss.pShifting_Current = numel(ss.pShifting_Tf);
                    else
                      ss.lockedUntil(atomIndex)   = ss.pShifting_Tf(ss.pShifting_Current);
                      ss.atomGenome{atomIndex}    = [ss.executingGene(atomIndex) ss.atomGenome{atomIndex}];
                      ss.executingGene{atomIndex} = [];
                      endGeneTime                 = ss.lockedUntil(atomIndex);
                    end
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.mitosisStatic
                    ss.r = r;
                    ss.rdyn = rdyn;
                    [ss newAtom] = atomMitosisStatic(ss, atomIndex, gene);
                    mitosisHappended = true;
                    r = ss.r;
                    rdyn = ss.rdyn;
                    %update the new atom's state
                    if isempty(ss.atomGenome{newAtom})
                      ss.lockedUntil(newAtom) = nan;
                      ss.executingGene{newAtom} = [];
                    else
                      ss.lockedUntil(newAtom)   = ss.t;
                      ss.executingGene{newAtom} = ss.atomGenome{newAtom}{1};
                      ss.atomGenome{newAtom}    = ss.atomGenome{newAtom}(2:end);
                      atomsToExpressGene(end+1) = newAtom;
                    end
                    modifiedPointsOrSprings = true;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.mirrorSprings
                    dualSpring    = [6;   5;   4;   3;   2;   1];
                    springIndexes = ss.atomSprings(atomIndex,:)';
                    Rs            = rdyn.r0(springIndexes);
                    %Ks            = ss.k(springIndexes);
                    deltaTs       = ss.atomState{atomIndex}.changeEdgeTime;
                    rdyn          = perturbateR(rdyn, r, ss.t, springIndexes, Rs(dualSpring));
                    %ss            = perturbToChange(ss, true, 'r', springIndexes, ss.t, deltaTs, Rs(dualSpring));
                    %ss            = perturbToChange(ss, true, 'k', springIndexes, ss.t, deltaTs, Ks(dualSpring));
                    modifiedPointsOrSprings   = true;
                    endGeneTime               = ss.t+deltaTs;
                    ss.lockedUntil(atomIndex) = endGeneTime;

                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.permuteShiftMatrix
                    switch gene.args{1}
                      case 'V1'
                        ss.atomState{atomIndex}.shiftMatrix(:,1) = ss.atomState{atomIndex}.shiftMatrix([2,1],1);
                      case 'V2'
                        ss.atomState{atomIndex}.shiftMatrix(:,2) = ss.atomState{atomIndex}.shiftMatrix([2,1],2);
                      case 'D1'
                        ss.atomState{atomIndex}.shiftMatrix([1,4]) = ss.atomState{atomIndex}.shiftMatrix([4,1]);
                      case 'D2'
                        ss.atomState{atomIndex}.shiftMatrix([2,3]) = ss.atomState{atomIndex}.shiftMatrix([3,2]);
                      case 'H1'
                        ss.atomState{atomIndex}.shiftMatrix(1,:) = ss.atomState{atomIndex}.shiftMatrix(1,[2,1]);
                      case 'H2'
                        ss.atomState{atomIndex}.shiftMatrix(2,:) = ss.atomState{atomIndex}.shiftMatrix(2,[2,1]);
                      case 'H'
                        ss.atomState{atomIndex}.shiftMatrix      = ss.atomState{atomIndex}.shiftMatrix(:,[2,1]);
                      case 'V'
                        ss.atomState{atomIndex}.shiftMatrix      = ss.atomState{atomIndex}.shiftMatrix([2,1],:);
                      otherwise
                        error('gene permuteShiftMatrix cannot have argument %s\n', mat2str(gene.args{1}));
                    end
            %       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       case geneTypes.changeSpringToShift
            %         if ~ismember(gene.args{1}, 1:6)
            %           error('gene changeSpringToShift cannot have argument %s\n', mat2str(gene.args{1}));
            %         end
            %         ss.atomState{atomIndex}.springToShift = gene.args{1};
            %       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       case geneTypes.changeSpringToShiftSign
            %         switch gene.args{1}
            %           case -1
            %             ss.atomState{atomIndex}.springToShiftSign = -1;
            %           case +1
            %             ss.atomState{atomIndex}.springToShiftSign = +1;
            %           case 0
            %             ss.atomState{atomIndex}.springToShiftSign = -ss.atomState{atomIndex}.springToShiftSign;
            %           otherwise
            %             error('gene changeSpringToShiftSign cannot have argument %s\n', mat2str(gene.args{1}));
            %         end
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.changeGeneDistribution
                    ss.atomState{atomIndex}.geneDistribution=sign(gene.args{1});
            %         newval = ss.atomState{atomIndex}.geneDistribution+gene.args{1};
            %         newval = max(min(newval,1),0); %ratio: max 1, min 0
            %        ss.atomState{atomIndex}.geneDistribution = newval;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.changeAttachedSpringsDistribution
                    ss.atomState{atomIndex}.attachedSpringsDistribution=sign(gene.args{1});
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.changeMitosisTime
                    newval = ss.atomState{atomIndex}.mitosisTime+gene.args{1};
                    newval = max(min(newval,10),1); %TODO: these limits are hard-coded: THEY MUST BE PARAMETERIZED!!!!
                    ss.atomState{atomIndex}.mitosisTime = newval;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.waitSomeTime
                    timeToWait = gene.args{1};
                    timeToWait = max(min(timeToWait,10),0); %TODO: these limits are hard-coded: THEY MUST BE PARAMETERIZED!!!!
                    endGeneTime=ss.t+timeToWait;
                    ss.lockedUntil(atomIndex) = endGeneTime;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.changeSpringR
                    relativeSpringsToSwitch = findSpringsRelativeToShift(ss.atomState{atomIndex}.shiftMatrix, gene.args{1}, false);
                    spring = ss.atomSprings(atomIndex,relativeSpringsToSwitch);
                    %maxR = min(max(ss.r(ss.atomSprings(atomIndex,:)))*5, 20);
                    %minR = 0; 
                    minR = 0;
                    thisr = rdyn.r0(ss.atomSprings(atomIndex,:));
                    maxR = max(thisr)*5; %TODO: these limits are hard-coded: THEY MUST BE PARAMETERIZED!!!!
                    factor = min(max(gene.args{2}, 0.5), 2); %TODO: these limits are hard-coded: THEY MUST BE PARAMETERIZED!!!!
                    newr = max(min(factor*thisr(relativeSpringsToSwitch), maxR), minR);

                    rdyn          = perturbateR(rdyn, r, ss.t, spring, newr);
                    modifiedPointsOrSprings = true;
                    endGeneTime=ss.t+ss.atomState{atomIndex}.changeEdgeTime;
                    ss.lockedUntil(atomIndex) = endGeneTime;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.changeSpringK
                    %error('This gene is disabled!!!');
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  case geneTypes.changeEdgeTime
                    chgT = gene.args{1}+ss.atomState{atomIndex}.changeEdgeTime;
                    chgT = max(min(chgT,10),0.1); %TODO: these limits are hard-coded: THEY MUST BE PARAMETERIZED!!!!
                    ss.atomState{atomIndex}.changeEdgeTime = chgT;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  otherwise 
                    error('gene for atom %g cannot be recognized: %s\n', atomIndex, gene.type);
                end
                %gene expression goes on if the expressed gene is instantaneous and
                %there remain genes to be expressed
                isNowSilent = isempty(ss.atomGenome{atomIndex});
                continuar = (endGeneTime==ss.t) && (~isNowSilent);
                if continuar
                  ss.executingGene{atomIndex} = ss.atomGenome{atomIndex}{1};
                  ss.atomGenome{atomIndex}    = ss.atomGenome{atomIndex}(2:end);
                end
                %if the atom has become silent, update its state accordingly
                if isNowSilent
                  ss.lockedUntil(atomIndex) = nan;
                  ss.executingGene{atomIndex} = [];
                end
              end
              k=k+1;
            end
            veryChanged = veryChanged || (k>1);
            
            if ss.useSAP && mitosisHappended
              ss = makeSAPSets(ss);
            end
  end;
end;
  if modifiedPointsOrSprings
    ss = prepareSimulationStructure(ss);
    [ss.shiftedPoints nevermind onlyShiftedIndexes] = unique(ss.pShifting_M);
    %this sparse matrix has 1 in the (i,j)-th position, meaning that the
    %i-th shifted point is in the j-th shifting
    % ss.shiftingMatrix = sparse(ss.pShifting_M, 1:numel(ss.pShifting_M), 1, size(ss.pos,1), numel(ss.pShifting_M));
    ss.shiftingMatrix = sparse(onlyShiftedIndexes, 1:numel(ss.pShifting_M), 1, numel(ss.shiftedPoints), numel(ss.pShifting_M));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = makeSAPSets(ss)
sets = zeros(size(ss.m), 'uint8');
forcedPoints = unique(ss.pShifting_M);
if numel(forcedPoints)>255
  error('Too many cells are dividing simultaneosly!!!!');
end
sets(forcedPoints) = uint8(1):numel(forcedPoints);
ss.sap = createSweepAndPrune(ss.pos, ss.stick.allr, sets, ss.sap.mode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [ss pointsModified] = finishPointShiftings(ss)
%do shifting makeup: points _M and _2 are to be fused
finishedShifts = ss.pShifting_Tf<=(ss.t+eps(ss.t));
pointsModified = any(finishedShifts);
if pointsModified
  shiftsToFuse = finishedShifts & ss.pShifting_fuse;
  %points _M are to be erased
  pointsToFuse_M = ss.pShifting_M(shiftsToFuse);
  %points _2 are to be added _M's contexts
  pointsToFuse_2 = ss.pShifting_2(shiftsToFuse);
  %remove the finished shifting vars
  for k=1:numel(ss.pShiftingVars)
    name = ss.pShiftingVars{k};
    ss.(name)(finishedShifts,:) = [];
  end

  %points to delete for each fusing
  toDelete = zeros(size(pointsToFuse_M));
  pointsToFuse_M_working = pointsToFuse_M;
  pointsToFuse_2_working = pointsToFuse_2;
  newIndexes = 1:size(ss.pos,1);
  %do it the heavy, slow but safe way
  for k=1:numel(pointsToFuse_M)
    current_delete   = pointsToFuse_M_working(k);
    current_receiver = pointsToFuse_2_working(k);
    if current_delete==current_receiver
      %toDelete{k} = [];
      continue
    end
    toDelete(k)      = current_delete;
    %if this _M is repeated, replace its ocuurences by _2, which will also
    %be erased
    pointsToFuse_M_working(pointsToFuse_M_working==current_delete) = current_receiver;
    %if this _M is also an _2, update the other shifting's _2
    pointsToFuse_2_working(pointsToFuse_2_working==current_delete) = current_receiver;
    %update re-index var
    newIndexes(newIndexes==current_delete) = current_receiver;
    newIndexes(current_delete)             = current_receiver;
    %add _M's atom list to _2's atom list
    ss.pointAtoms{current_receiver} = [ss.pointAtoms{current_receiver}; ss.pointAtoms{current_delete}];
  end
  toDelete = unique(toDelete);
  if (toDelete(1)==0)
    toDelete = toDelete(2:end);
  end
  %update the indexes in all point index vars
  for k=1:numel(ss.pointIndexVars)
    name = ss.pointIndexVars{k};
    ss.(name) = auxReindex(newIndexes, ss.(name));
  end
  %remove the points
  ss = removePoints(ss, toDelete);
    
  if isempty(ss.pShifting_Tf)
    ss.pShifting_Current = [];
  else
    ss.pShifting_Current = max(ss.pShifting_Tf);
  end
  
  if ss.useSAP
    ss = makeSAPSets(ss);
  end
  
%   CORREGIR: ASI AUTOMATIZADO ES MU BONITO, PERO NO TIENE EN CUENTA QUE UN MISMO PUNTO _2 PUEDE RECIBIR CONTRIBUCIONES DE VARIOS _M 
%   CORREGIR: LOS PUNTOS _M PUEDEN IR A VARIOS _2 (EN UN PAR DISTINTO EN CADA CASO): HAY QUE FUSIONAR CON CUIDADO. ES MAS, PUEDE HABER VARIAS fUSIONES eNcADENADAS _m->_2_m->_2 
%   %contexts of _M and _2 are suppossed to be identical, except for the atoms
%   %they are attached to. So, all of _M's context is lost but its field
%   %"pointAtoms"
%   ss.pointAtoms(pointsToFuse_2) = cellfun(@(a,b) [a;b], ss.pointAtoms(pointsToFuse_2), ss.pointAtoms(pointsToFuse_M), 'UniformOutput', false);
%   %prepare the new indexes
%   newIndexes = 1:size(ss.pos,1);
%   newIndexes(pointsToFuse_M) = pointsToFuse_2;
%   %update the indexes in all point index vars
%   for k=1:numel(ss.pointIndexVars)
%     name = ss.pointIndexVars{k};
%     ss.(name) = auxReindex(newIndexes, ss.(name));
%   end
%   %remove the old points
%   ss = removePoints(ss, pointsToFuse_M);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function returns the index(es) of one (two) springs, relative to the
%shift matrix. If both==true, two springs are returned: the specified and
%its opposite.
function springsToSwitch = findSpringsRelativeToShift(shiftMatrix, spring, both)
switch spring
  case 'V1'
    ends = shiftMatrix(:,1); %this is the same as the spring to shift
  case 'V2'
    ends = shiftMatrix(:,2); %this is the same as the spring to receive the shift 
  case 'H1'
    ends = shiftMatrix(1,:);
  case 'H2'
    ends = shiftMatrix(2,:);
  case 'D1'
    ends = shiftMatrix([1 4]);
  case 'D2'
    ends = shiftMatrix([2 3]);
end
ends = sort(ends);
springsToSwitch = twoCombinationToIndex(ends(1),ends(2),4);
if both
  springsToSwitch = [springsToSwitch; 6-springsToSwitch+1];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ss endTime newAtom] = atomMitosisSpringShift(ss, atomToDivide, gene)

ss = handleForcementAfterMitosis(ss, atomToDivide);
state = ss.atomState{atomToDivide};

%split the genome
[ss.atomGenome{atomToDivide} genomeB] = DAsplitGenome(ss.atomGenome{atomToDivide}, state.geneDistribution);

%find the parent's points
parentPoints  = ss.atomPoints(atomToDivide, :);
parentSprings = ss.atomSprings(atomToDivide,:);

springsK = ss.k(parentSprings);
springsR = ss.r(parentSprings);

%create child
[ss newAtom] = addAtomComplex(ss, zeros(1,4), ss.pos(parentPoints,:), springsK, springsR, {genomeB}, {state});
%make sure that, initially, new points are exactly equal to parent points
newPoints = ss.atomPoints(newAtom,:);
for k=1:numel(ss.pointVars)
  name = ss.pointVars{k};
  if ~strcmp(name, 'pointAtoms')
    ss = reassignField(ss, name, 3, parentPoints, newPoints);
  end
end

localSprings      = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
%In an atom, a spring is between the two points. This gives the index for
%the dual spring: the spring running along the other two points
dualSpring        = [6;   5;   4;   3;   2;   1];

%pointsA1 go from points B1 to pointsB2. 
%pointsB2 go from points A2 to pointsA1.
%
pointsA1 = ss.atomPoints(atomToDivide, state.shiftMatrix(:,1));
pointsA2 = ss.atomPoints(atomToDivide, state.shiftMatrix(:,2));
pointsB1 = ss.atomPoints(newAtom,      state.shiftMatrix(:,1));
pointsB2 = ss.atomPoints(newAtom,      state.shiftMatrix(:,2));
% %points in the shifting spring from the original atom
% pointsA1 = ss.springEnds(ss.atomSprings(atomToDivide,           state.springToShift), :);
% %points in the shifting spring from the new atom
% pointsB1 = ss.springEnds(ss.atomSprings(newAtom,                state.springToShift), :);
% %points in the shifting spring's dual one, from the original atom
% pointsA2 = ss.springEnds(ss.atomSprings(atomToDivide,dualSpring(state.springToShift)),:);
% %points in the shifting spring's dual one, from the new atom
% pointsB2 = ss.springEnds(ss.atomSprings(newAtom,     dualSpring(state.springToShift)),:);

endTime = ss.t+state.mitosisTime;

%DISTRIBUTION OF OTHER ATOMS' VERTICES
%CASE -1: all vertices go to cell A. As they are already there, there 
%         is no need to change anything.
%CASE  0: the distribution makes sure that the other atoms remain attached
%         to the vertices that are not shifted to be joined. That is to
%         say, other atoms are to be attached to pointsB1 and pointsA2.
%         they are already attached to pointsA2, so we only need to care
%         about transfering other atoms from pointsA1 to pointsB1.
%CASE +1: all vertices go to cell B.
switch state.attachedSpringsDistribution
  case -1; fromPoints = [];                   toPoints = [];
  case 0;  fromPoints =  pointsA1;            toPoints =  pointsB1;
  case 1;  fromPoints = [pointsA1; pointsA2]; toPoints = [pointsB1; pointsB2];
  otherwise; error('the distribution of attached springs is ordered by three options, and %s is not one of them\n', mat2str(state.attachedSpringsDistribution));
end
for k = 1:numel(fromPoints)
  fromPoint   = fromPoints(k);
  toPoint     = toPoints(k);
  otherAtoms    = ss.pointAtoms{fromPoint}(ss.pointAtoms{fromPoint}~=atomToDivide);
  %rearrange springs' indexes
  springsAttachedToPoint = ss.springEnds==fromPoint;
  springsAttachedToPoint(parentSprings,:) = false;
  ss.springEnds(springsAttachedToPoint)   = toPoint;
  %the toPoint acquires all atoms but the parent
  ss.pointAtoms{toPoint} = [newAtom; otherAtoms];
  %the fromPoint only keeps one atom: the parent
  ss.pointAtoms{fromPoint} = atomToDivide;
  %the other atoms are switched from the parent's to the daughter's point
  for m=1:numel(otherAtoms)
    atom = otherAtoms(m);
    ap = ss.atomPoints(atom,:);
    ss.atomPoints(atom,ap==fromPoint) = toPoint;
  end
end


% %DISTRIBUTION OF OTHER ATOMS' VERTICES
% %easy to become unstable: what if an atom becomes stretched to only two
% %points?
% for k=1:4
%   %for each parent's point, we check whether all other springs (and other
%   %atoms) must be switched to the daughter's point
%   if state.attachedSpringsDistribution(k)=='2'
%     parentPoint   = parentPoints(k);
%     daughterPoint = ss.atomPoints(newAtom,k);
%     otherAtoms    = ss.pointAtoms{parentPoint}(ss.pointAtoms{parentPoint}~=atomToDivide);
%     %rearrange springs' indexes
%     springsAttachedToPoint = ss.springEnds==parentPoint;
%     springsAttachedToPoint(parentSprings,:) = false;
%     ss.springEnds(springsAttachedToPoint)   = daughterPoint;
%     %the daughter's point acquires all atoms but the parent
%     ss.pointAtoms{daughterPoint} = [newAtom; otherAtoms];
%     %the parent point only keeps one atom: the parent
%     ss.pointAtoms{parentPoint} = atomToDivide;
%     %the other atoms are switched from the parent's to the daughter's point
%     for m=1:numel(otherAtoms)
%       atom = otherAtoms(m);
%       ap = ss.atomPoints(atom,:);
%       ss.atomPoints(atom,ap==parentPoint) = daughterPoint;
%     end
%   end
% end

% if state.springToShiftSign<=0
%   pointOrder = [2 1];
% else
%   pointOrder = [1 2];
% end
%ss.pShifting_M  = [ss.pShifting_M;  pointsA1(pointOrder)'; pointsB2(pointOrder)'];

%two shiftings per mitosis, to be symmetric
ss.pShifting_M  = [ss.pShifting_M;  pointsA1'; pointsB2'];
ss.pShifting_1  = [ss.pShifting_1;  pointsB1'; pointsA2'];
ss.pShifting_2  = [ss.pShifting_2;  pointsB2'; pointsA1'];
ss.pShifting_fuse = [ss.pShifting_fuse; true(2,1); false(2,1)];
ss.pShifting_T0 = [ss.pShifting_T0; repmat(ss.t,             4,1)];
ss.pShifting_Tf = [ss.pShifting_Tf; repmat(endTime,          4,1)];
ss.pShifting_TD = [ss.pShifting_TD; repmat(-state.mitosisTime,4,1)];
% 
%_V are speeds designed to make points _M and _2 to meet at time _Tf
% ss.pShifting_V  = [ss.pShifting_V;  (ss.pos(pointsA2,:)-ss.pos(pointsB1,:))/state.mitosisTime/2];

ss.lockedUntil([atomToDivide newAtom]) = endTime;
ss.executingGene{atomToDivide} = gene;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ss newAtom] = atomMitosisStatic(ss, atomToDivide, gene)

ss = handleForcementAfterMitosis(ss, atomToDivide);
state = ss.atomState{atomToDivide};

%split the genome
[ss.atomGenome{atomToDivide} genomeB] = DAsplitGenome(ss.atomGenome{atomToDivide}, state.geneDistribution);

%find the parent's points
parentPoints  = ss.atomPoints(atomToDivide, :);
parentSprings = ss.atomSprings(atomToDivide,:);

springsK = ss.k(parentSprings);
springsR = ss.r(parentSprings);

%create child
[ss newAtom] = addAtomComplex(ss, parentPoints, zeros(0,4), springsK, springsR, {genomeB}, {state});

ss.lockedUntil([atomToDivide newAtom]) = ss.t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%after a mitosis event, delayed enforcement must be handled
function ss = handleForcementAfterMitosis(ss, atomIndex)

stateAS = ss.atomState{atomIndex}.afterStabilized;

%TAKE CARE OF STRESS ENFORCEMENTS
concluded = find(stateAS.delayedSwitchesDoForcement(:,2)<=0);
%this is handled in a loop (non-vectorized) to be safe: if several delayed
%enforcements over the same spring end by the same time, we must make sure
%that negations are performed the right amount of times. 
for k=1:numel(concluded)
  springToSwap = stateAS.delayedSwitchesDoForcement(concluded(k),1);
  stateAS.doForcement(springToSwap) = ~stateAS.doForcement(springToSwap);
end
stateAS.delayedSwitchesDoForcement(concluded,:) = [];
stateAS.delayedSwitchesDoForcement(:,2) = stateAS.delayedSwitchesDoForcement(:,2)-1;

%TAKE CARE OF TYPES OF STRESS ENFORCEMENTS
concluded = find(stateAS.delayedSwitchesForcementType(:,2)<=0);
%this is handled in a loop (non-vectorized) to be safe: if several delayed
%enforcements over the same spring end by the same time, we must make sure
%that negations are performed the right amount of times. 
for k=1:numel(concluded)
  springToSwap = stateAS.delayedSwitchesForcementType(concluded(k),1);
  stateAS.forcementType(springToSwap) = ~stateAS.forcementType(springToSwap);
end
stateAS.delayedSwitchesForcementType(concluded,:) = [];

ss.atomState{atomIndex}.afterStabilized = stateAS;      
stateAS.delayedSwitchesForcementType(:,2) = stateAS.delayedSwitchesForcementType(:,2)-1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
