%this recorder saves the simulation in memory
classdef MemRecorder < Recorder
  
  properties
    dyn;         %dynamic state
    dynTs;       %dynamic state's time
    notdyn;      %non-dynamic state
    notdynTs;    %dynamics state's time
    incdyn;      %increment step
    incnotdyn;   %increment step
    dynIndex;    %first empty slot in 'dyn'
    notdynIndex; %first empty slot in 'notdyn'
  end
  
  methods
    %constructor
    function rec = MemRecorder(ss, dynspace, notdynspace, incdyn, incnotdyn)
      rec = rec@Recorder(ss);
      rec.TYPE = 'MemRecorder';
      rec.incdyn = incdyn;
      rec.incnotdyn = incnotdyn;
      rec.dyn = cell(dynspace,1);
      rec.dynTs = zeros(dynspace,1);
      rec.notdyn = cell(notdynspace,1);
      rec.notdynTs = zeros(notdynspace,1);
      rec.dynIndex = 1;
      rec.notdynIndex = 1;
    end
    
    function rec = recordAllState(rec, ss)
      %if ss.t is the same as the last rec.notdynTs, overprint it
      if (rec.notdynIndex>1) && (rec.notdynTs(rec.notdynIndex-1)==ss.t)
        rec.notdyn{rec.notdynIndex-1} = extractAllSystemState(rec, ss);
        return
      end
      %make sure that room is available
      if rec.notdynIndex>numel(rec.notdyn)
        rec.notdyn{end+rec.incnotdyn,:} = [];
        rec.notdynTs(end+rec.incnotdyn,:) = 0;
      end
      %extract state and record it
      rec.notdynTs(rec.notdynIndex) = ss.t;
      rec.notdyn{rec.notdynIndex} = extractAllSystemState(rec, ss);
      rec.notdynIndex = rec.notdynIndex+1;
      recordingCallback(rec, ss); %do not forget to call recordingCallback!!!
    end
    
    function rec = recordDynState(rec, ss, Ts, states)
      %make sure that room is available
      while (rec.dynIndex+numel(Ts))>numel(rec.dyn)
        rec.dyn{end+rec.incdyn,:} = [];
        rec.dynTs(end+rec.incdyn,:) = 0;
      end
      %extract states and record them
      indexes = rec.dynIndex+(0:(numel(Ts)-1));
      rec.dynTs(indexes) = Ts;
      rec.dyn(indexes) = mat2cell(states, size(states,1), ones(1, size(states,2)));
      rec.dynIndex = indexes(end) + 1;
      recordingCallback(rec, ss, Ts, states); %do not forget to call recordingCallback!!!
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
    
    %find all recorded timestep within the given time interval
    function Ts = findTimestepsBetween(rec, t1, t2)
      allTs = allTimeSteps(rec);
      Ts = allTs((allTs>=t1) & (allTs<=t2));
    end
    
    %find all recorded timesteps in this MemRecorder
    function Ts = allTimeSteps(rec)
      %Ts = sort(unique([rec.dynTs(1:(rec.dynIndex-1)); rec.notdynTs(1:(rec.notdynIndex-1))]));
      %elements are presumed to be pre-sorted
      Ts = mergesorted(rec.dynTs(1:(rec.dynIndex-1)), rec.notdynTs(1:(rec.notdynIndex-1)));
      %remove repeated elements
      Ts(diff(Ts)==0) = [];
    end
      
    %Recorder interface, but useless for MemRecorder
    function rec = flushRecording(rec); end

  end
    
  
end