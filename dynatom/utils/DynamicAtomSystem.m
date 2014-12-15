%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = DynamicAtomSystem(geneTypes, otherArgs)

ss = BasicSystem(otherArgs{:});

ss.postProcess     = @magicallyMovedPoints;
ss.chgStep         = @DAProcessGenes;
ss.brownianFractor = false;
ss.useAnisodamp    = false;
ss.useSAP          = false;

ss.atomPoints             = zeros(0,4); %atom's four points
ss.atomSprings            = zeros(0,6); %atom's six springs, ordered: 12 13 14 23 24 34
ss.atomGenome             = cell(0,1);
ss.executingGene          = cell(0,1);
ss.lockedUntil            = zeros(0,1);
ss.atomState              = cell(0,1);
afterStabilized           = struct(... %substruct of atomDefaultState
  'doForcement',                  false(6,1), ... %force each spring if the flag is true
  'forcementType',                false(6,1), ... %if doForcement is true, it means to (false:increase stress)/(true:reverse stress)
  'delayedSwitchesDoForcement',   zeros(0,2), ... %this will be a list of pairs (spring, amounts of generations). After each mitosis, every amount will be decreased by one. When an amount reachs 0, it is deleted and doForcement(spring)   is switched
  'delayedSwitchesForcementType', zeros(0,2)  ... %this will be a list of pairs (spring, amounts of generations). After each mitosis, every amount will be decreased by one. When an amount reachs 0, it is deleted and forcementType(spring) is switched
  );
atomDefaultState          = struct(... %atom specific vars
  'geneDistribution',            -1,    ... % [-1,0,1] == genes [all to 1st daughter, distributed among two daughters, all to 2nd daughter]
  'attachedSpringsDistribution', 0,      ... % [-1,0,1] == attached go to [all to 1,distributed in 1 & 2 in non-shifted vertices,all to 2]
  'shiftMatrix',                 [1 3; 2 4], ... [a b;c d] the springs to shift are "ac" in the 1st daughter and "bd" in the 2nd daughter
  ...%'compressionToDevelopment',    repmat(0.7,6,1), ...
  ...%'stretchingToDevelopment',     repmat(1.3,6,1), ...
  ...%'changeToDevelopment',         [1 1 1 1 1 1]', ...
  'mitosisTime',                 3,      ... % time units
  'mitosisPreTime',              5,      ... % time units
  'changeEdgeTime',              3,      ... % length units
  'afterStabilized',             afterStabilized ...
  );


ss.atomConstants = struct;

ss.geneTypes              = geneTypes;
ss.atomVars               = {'atomPoints', 'atomSprings', 'atomGenome', 'atomState',        'executingGene', 'lockedUntil'};
ss.atomDefaultVals        = {zeros(1,4),   zeros(1,6),    {[]},         {atomDefaultState}, {[]},             NaN};
ss.pointIndexVars{end+1}  = ss.atomVars{1};
ss.springIndexVars{end+1} = ss.atomVars{2};

ss.pointAtoms                = cell(0,1);  %for each point, the atoms it belongs to
ss.pointVars{end+1}          = 'pointAtoms';
ss.pointDefaultVals{end+1}   = {[]};

ss.springAtom                = zeros(0,1); %for each spring, the atom it belongs to
ss.springVars{end+1}         = 'springAtom';
ss.springDefaultVals{end+1}  = 0;

ss.atomIndexVars             = {ss.pointVars{end}, ss.springVars{end}};

ss.devTimes = struct('T_after_growth', {[]}, 'T_after_stabilization', {[]});

ss.pShifting_Current = [];
ss.pShifting_M  = zeros(0,1);
ss.pShifting_1  = zeros(0,1);
ss.pShifting_2  = zeros(0,1);
ss.pShifting_fuse = false(0,1);
% ss.pShifting_V  = zeros(0,size(ss.pos,2));
ss.pShifting_T0 = zeros(0,1);
ss.pShifting_Tf = zeros(0,1);
ss.pShifting_TD = zeros(0,1);
% ss.pShiftingVars = {'pShifting_M', 'pShifting_1', 'pShifting_2', 'pShifting_V', 'pShifting_T0', 'pShifting_Tf'};
ss.pShiftingVars = {'pShifting_M', 'pShifting_1', 'pShifting_2', 'pShifting_fuse', 'pShifting_T0', 'pShifting_Tf', 'pShifting_TD'};
ss.shiftingMatrix = []; %working var: sparse matrix containing 1 in the (i,j)-th position, meaning that the i-th point is in the j-th shifting
ss.shiftedPoints  = []; %working var: UNIQUE-ED version of ss.pShifting_M
ss.pointIndexVars = [ss.pointIndexVars, ss.pShiftingVars(1:3)];

newSystemVars   = {'tick'};
%ss.systemVars   = [ss.systemVars, newSystemVars];

ss.allStateVars = [ss.allStateVars, ss.atomVars, 'devTimes', 'atomDefaultVals', 'geneTypes', 'atomConstants', ss.pShiftingVars{:}, 'pShifting_Current', ss.pointVars{end}, ss.springVars{end}, newSystemVars];

ss.fieldSets    = [ss.fieldSets, 'atomVars', 'atomIndexVars', 'pShiftingVars'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dpos, dvel] = magicallyMovedPoints(ss, dpos, dvel, r, t)
%Given several triplets of points (M,1,2), control the points M's kinematics:  
%   -before the interval, their movements are stuck to points 1's
%   -after  the interval, their movements are stuck to points 2's
%   -during the interval, they smoothly move along the lines from points 1
%    to points 2 
%
%dpos==vel, dvel==acel

if ~isempty(ss.pShifting_M)
  %get data
  pos = ss.pos;
  pM = ss.pShifting_M;
  if isfield(ss, 'shiftedPoints')
    shifted = ss.shiftedPoints;
    shiftingMatrix = ss.shiftingMatrix;
  else
    [shifted nevermind onlyShiftedIndexes] = unique(ss.pShifting_M);
    %this sparse matrix has 1 in the (i,j)-th position, meaning that the
    %i-th shifted point is in the j-th shifting
    % ss.shiftingMatrix = sparse(ss.pShifting_M, 1:numel(ss.pShifting_M), 1, size(ss.pos,1), numel(ss.pShifting_M));
    shiftingMatrix = sparse(onlyShiftedIndexes, 1:numel(ss.pShifting_M), 1, numel(shifted), numel(ss.pShifting_M));
  end
  p2 = ss.pShifting_2;
  Tf = ss.pShifting_Tf;
  TD = ss.pShifting_TD;
% PREVIOUSLY, THE CODE WAS PREPARED TO HANDLE TRIPLETS OUTSIDE ITS LIFE
% CYCLE. HOWEVER, THE CODE WAS UNABLE TO HANDLE POINTS ATTACHED TO SEVERAL
% TRIPLETS. SINCE IT IS QUITE IMPORTANT TO ACHIEVE THIS, WE DROP THE
% ABILITY TO HANDLE TRIPLETS OUTSIDE ITS LIFE CYCLE (TOO EXPENSIVE TO DO
% IT PROPERLY). SO, WE ***RELY*** ON THE SYSTEM TO HALT THE SIMULATION
% WHENEVER A TRIPLET HAS FINISHED ITS LIFE CYCLE, TO WIPE IT OUT. THIS ALSO
% IS BETTER, SINCE THE MAKING OF THE INCIDENCE (SPARSE) MATRIX CAN BE
% OUTSOURCED
    %triplet_dpos = (dpos(p2,:)+dpos(pM,:))/2 + bsxfun(@rdivide, pos(p2,:)-pos(pM,:), (Tf-t+1e-3)./(-TD)); %2.02  
    triplet_dpos = pos(p2,:)-pos(pM,:);
    dv = (Tf-t+1e-3)./(-TD);
    triplet_dpos(:,1) = triplet_dpos(:,1)./dv;
    triplet_dpos(:,2) = triplet_dpos(:,2)./dv;
    triplet_dpos = triplet_dpos + (dpos(p2,:)+dpos(pM,:))/2;
    dpos(shifted,:) = shiftingMatrix*triplet_dpos;
    %cancel accelerations
    dvel(shifted,:) = shiftingMatrix*((dvel(p2,:)+dvel(pM,:))/2);
end
end