%This abstract class is the interface which BasicSystem uses to record the
%dynamics
classdef Recorder
  
  properties
     TYPE;
     callback;
  end
  
  properties (GetAccess = protected, SetAccess = protected)
    emptyState;
    fieldSets;
  end
  
  methods (Abstract)
    %record all system state
    rec = recordAllState(rec, ss)
    %record dynamic states for several times
    rec = recordDynState(rec, ss, Ts, states)
    %for each time in Ts, apply the given function
    playSimulation(rec, Ts, fun);
    %find all recorded timsteps between two times
    Ts = findTimestepsBetween(rec, t1, t2);
    %find all recorded timesteps in this Recorder
    Ts = allTimeSteps(rec);
    %signaler for disk and cache based recorders
    rec = flushRecording(rec);
  end

  methods
    function rec = Recorder(ss)
      rec.fieldSets = getBasicSystemFieldSets(ss, true);
      rec.TYPE = 'Recorder';
      rec.callback = [];
    end
    %for each timestep in the interval [t1,t2], apply the given function
    function playSimulationBewteen(rec, t1, t2, fun)
      Ts = findTimestepsBetween(rec, t1, t2);
      playSimulation(rec, Ts, fun);
    end
    %return the system for a given T
    function ss = getTimeStep(rec, T)
      ss = [];
      playSimulation(rec, T, @assignSS);
      function assignSS(t, system); ss = system; end
    end
  end
  
  methods (Access = protected)
    %do callbacks, if required. This **MUST** be called by recordAllState
    %and recordDynState if necessary
    function recordingCallback(rec, varargin) %, ss, Ts, states)
      if ~isempty(rec.callback) % && ishandle(rec.callback)
        handle = rec.callback;
        handle(varargin{:});
      end
    end
    
    %extracts info from BasicSystem into a structure
    function st = extractAllSystemState(rec, ss)
      if isempty(rec.emptyState)
      %   %THIS IS VALID ONLY IF NO SYSTEM STATE IS ALLOWED TO BE NESTED IN
      %   %SUBSTRUCTURES
      %   aux = reshape([ss.allStateVars;cell(1,numel(ss.allStateVars))],[],1);
      %   rec.emptyState = struct(aux{:});
        rec.emptyState = struct;
        for k=1:numel(ss.allStateVars)
          rec.emptyState = reassignField(rec.emptyState, ss.allStateVars{k}, 0, []);
        end
      end
      % st = rec.emptyState;
      % for k=1:numel(ss.allStateVars)
      %   %THIS IS VALID ONLY IF NO SYSTEM STATE IS ALLOWED TO BE NESTED IN
      %   %SUBSTRUCTURES
      %   campo = ss.allStateVars{k};
      %   st.(campo) = ss.(campo);
      % end
      st = copyFields(ss, rec.emptyState,  ss.allStateVars);
    end
    
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
    end
  end

end