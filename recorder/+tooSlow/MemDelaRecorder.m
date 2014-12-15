%this recorder saves the simulation in memory
classdef MemDelaRecorder < MemRecorder
  
  properties
    dyndela;    %delaunay triangulation
  end
  
  methods
    %constructor
    function rec = MemDelaRecorder(ss, dynspace, notdynspace, incdyn, incnotdyn)
      rec = rec@MemRecorder(ss, dynspace, notdynspace, incdyn, incnotdyn);
      rec.dyndela = cell(dynspace,1);
    end
    
    function rec = recordAllState(rec, ss)
      rec = recordAllState@MemRecorder(rec, ss);
      rec.notdyn{rec.notdynIndex-1}.delaunay = ss.delaunay(:,[1 2]);
    end
    
    function rec = recordDynState(rec, ss, Ts, states)
      rec = recordDynState@MemRecorder(rec, ss, Ts, states);
      if numel(Ts)>1
        error('This must be called yomestep bny timestep!!!');
      end
      %make sure that room is available
      if numel(rec.dyndela)<numel(rec.dyn)
        rec.dyndela{numel(rec.dyn)} = [];
      end
      rec.dyndela{rec.dynIndex-1} = ss.delaunay(:,[1 2]);
    end
    
    function playSimulation(rec, Ts, fun)
      neededFieldSets = {}; %all fields %{'dynVars', 'flagVars', 'numFlagVars', 'cumsumFlagVars', 'doubledFlagVars', 'doubledNumFlagVars'};
      %locate Ts members which are stored as 'dyn' timesteps
      [nevermind, dynloc] = ismember(Ts, rec.dynTs(1:(rec.dynIndex-1)));
      %define boundaries for 'notdyn' states
      limits = [rec.notdynTs(1:(rec.notdynIndex-1)), [rec.notdynTs(2:(rec.notdynIndex-1)); Inf]];
      %define initially selected 'notdyn'
      structureK = 1;
      structure = addFieldSets(rec, rec.notdyn{structureK}, neededFieldSets{:});
      %get customized selector for this 'notdyn' state
      structure = makeDynamicVars(makeSimulationAmounts(structure));
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
          structure = makeDynamicVars(makeSimulationAmounts(structure));
        end
        %let's rock!
        if rec.notdynTs(structureK)==t
          %if the timestep is stores as 'notdyn', use it directly. We
          %re-load 'structure', since it may have been modified previously
          structure = addFieldSets(rec, rec.notdyn{structureK}, neededFieldSets{:});
          %get customized selector for this 'notdyn' state
          structure = makeSimulationAmounts(structure);
          structure = makeDynamicVars(structure);
          fun(t, structure);
        elseif dynloc(k)~=0
          %if the timestep is recorded in 'dyn', retrieve dynamical vars
          %and dump them in the 'notdyn' to customize it. The key is that,
          %for each structure, all the 'dyn' states comprise the same set
          %of dynamical vars, so the same structure can be used over and
          %over
          structure = dumpEvolvedState(structure, rec.dyn{dynloc(k)});
          structure.delaunay = rec.dyndela{dynloc(k)};
          fun(t, structure);
        else
          %this seems to be an error: t is not a timestep in 'rec'
          fun(t, []);
        end
%         if mod(k,500)==0
%           fprintf('joer, que lento %d\n', k);
%         end
      end
    end
    
  end
    
  
end