%for each time in Ts, apply the given function
function playSimulation(rec, Ts, fun)
neededFieldSets = {}; %all fields %{'dynVars', 'flagVars', 'numFlagVars', 'cumsumFlagVars', 'doubledFlagVars', 'doubledNumFlagVars'};
%locate Ts members which are stored as 'dyn' timesteps
[nevermind, dynloc] = ismember(Ts, rec.dynTs(1:(rec.dynIndex-1)));
% edges = (rec.dynTs(1:(rec.dynIndex-2))+rec.dynTs(2:(rec.dynIndex-1)))/2;
% edges = [edges(1)-(edges(2)-edges(1)); edges; edges(end)+(edges(end)-edges(end-1))];
% [dynloc dynloc] = histc(Ts, edges);
% dynloc(dynloc==numel(edges))=0;
%define boundaries for 'notdyn' states
limits = [rec.notdynTs(1:(rec.notdynIndex-1)), [rec.notdynTs(2:(rec.notdynIndex-1)); Inf]];
%define initially selected 'notdyn'
structureK = 1;
structure = addFieldSets(rec, rec.notdyn{structureK}, neededFieldSets{:});
%get customized selector for this 'notdyn' state
structure = makeSimulationAmounts(structure); %makeDynamicVars(makeSimulationAmounts(structure));
%for each timestep to play
for k=1:numel(Ts)
  t = Ts(k);
  %chech whether it is within the time interval for the current
  %'notdyn'
  if (t<limits(structureK,1)) || (t>=limits(structureK,2))
    %if not, find the desired 'notdyn'
    structureK = find((t>=limits(:,1)) & (t<limits(:,2)));
    if numel(structureK)==0;
      error('TIME %g OUT OF RANGE (min:%g, max:%g)', t, limits(1,1), limits(end,2));
    end
    if numel(structureK)>1;
      error('MULTIPLE PERIODS MATCHED');
    end
    structure = addFieldSets(rec, rec.notdyn{structureK}, neededFieldSets{:});
    %get customized selector for this 'notdyn' state
    structure = makeSimulationAmounts(structure); %makeDynamicVars(makeSimulationAmounts(structure));
  end
  %let's rock!
  if abs(rec.notdynTs(structureK)-t)<(128*eps(t))
    %if the timestep is stores as 'notdyn', use it directly. We
    %re-load 'structure', since it may have been modified previously
    structure = addFieldSets(rec, rec.notdyn{structureK}, neededFieldSets{:});
    %get customized selector for this 'notdyn' state
    structure = makeSimulationAmounts(structure);
    %structure = makeDynamicVars(structure);
    fun(t, structure);
  elseif dynloc(k)~=0
    %if the timestep is recorded in 'dyn', retrieve dynamical vars
    %and dump them in the 'notdyn' to customize it. The key is that,
    %for each structure, all the 'dyn' states comprise the same set
    %of dynamical vars, so the same structure can be used over and
    %over
    structure = dumpEvolvedState(structure, rec.dyn{dynloc(k)});
    fun(t, structure);
  else
    %this seems to be an error: t is not a timestep in 'rec'
    fun(t, []);
  end
%         if mod(k,500)==0
%           fprintf('joer, que lento %d\n', k);
%         end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%adds field sets (metadata) to structure. It can add only the sets
%whose names are provided. If no names are provided, all sets are added
function st = addFieldSets(rec, st, varargin)
  fs = rec.fieldSets;
  if isempty(varargin)
    transferred = fs.fieldSets;
  else
    transferred = varargin;
  end
  for k=1:numel(transferred)
    st.(transferred{k}) = fs.(transferred{k});
  end

 