%this recorder only saves the energy of the simulation
classdef EnergyRecorder < Recorder
  
  properties
    e;
    ts;
    inc;      %increment step
    index;    %first empty slot in 'e'
  end
  properties (GetAccess = protected, SetAccess = protected)
    straightforward; %flag to tell if harvesting energy information is straightforward
  end
  
  methods
    %constructor
    function rec = EnergyRecorder(ss, space, inc)
      rec = rec@Recorder(ss);
      rec.TYPE = 'EnergyRecorder';
      rec.inc = inc;
      rec.e = zeros(space,1);
      rec.ts = zeros(space,1);
      rec.index = 1;
      rec.straightforward = false;
    end
    
    function rec = recordAllState(rec, ss)
      %if ss.t is the same as the last rec.ts, overprint it
      if (rec.index>1) && (ss.t==rec.ts(rec.index-1))
        rec.e(rec.index-1) = calculateEnergy(ss);
        return
      end
      %make sure that room is available
      if rec.index>numel(rec.e)
        rec.e(end+rec.inc,:)  = 0;
        rec.ts(end+rec.inc,:) = 0;
      end
      %extract state and record it
      rec.ts(rec.index) = ss.t;
      rec.e(rec.index)  = calculateEnergy(ss);
      rec.index = rec.index+1;
      %to harvest energy information is straightforward if all of these
      %hold:
      %   -nº of dimensions is 2
      %   -all points are dynamical_p
      %   -no point is dynamical_m, no spring is neither dynamical_r nor dynamical_k  
      rec.straightforward = (size(ss.pos,2)==2) && all(ss.dynamical_p) && (~(any(ss.dynamical_r) || any(ss.dynamical_k) || any(ss.dynamical_m) ));
      recordingCallback(rec, ss); %do not forget to call recordingCallback!!!
    end
    
    function rec = recordDynState(rec, ss, Ts, states)
      %make sure that room is available
      while (rec.index+numel(Ts))>numel(rec.e)
        rec.e(end+rec.inc,:)  = 0;
        rec.ts(end+rec.inc,:) = 0;
      end
      indexes1 = rec.index;
      indexes2 = rec.index+numel(Ts)-1;
      %if the first Ts is the same as the last rec.ts, overprint it
      if (rec.index>1) && (Ts(1)==rec.ts(rec.index-1))
        indexes1 = indexes1-1;
        indexes2 = indexes2-1;
      end
      rec.ts(indexes1:indexes2) = Ts;
      if rec.straightforward
        %fast way
        npoints = size(ss.pos,1);
        svel1   = ss.selector.vel(1,1);
        svel2   = ss.selector.vel(end,2);
        %now, extract posX and posY
        selPX = ss.selector.pos(1,1):ss.selector.pos(1,2);
        selPX = reshape(selPX, numel(selPX)/2, 2);
        selPY = selPX(:,2);
        selPX = selPX(:,1);
        %make them relative to spring ends
        selPX1 = selPX(ss.springEnds(:,1));
        selPX2 = selPX(ss.springEnds(:,2));
        selPY1 = selPY(ss.springEnds(:,1));
        selPY2 = selPY(ss.springEnds(:,2));
        %now, calculate energy in a vectorized way
        vp = realpow(states(svel1:svel2, :), 2);
        e_kinetic = sum(bsxfun(@times, vp(1:npoints,:)+vp((npoints+1):end,:), ss.m));
        springLengths = realsqrt(realpow(states(selPX2,:)-states(selPX1,:), 2) + realpow(states(selPY2,:)-states(selPY1,:), 2));
        springDisplacements = bsxfun(@minus, ss.r, springLengths);
        e_potential = sum(bsxfun(@times, ss.k, realpow(springDisplacements,2)));
        rec.e(indexes1:indexes2) = (e_kinetic+e_potential)/2;
        rec.index = indexes2+1;
      else
        %do it the veeeeeery slow way
        for k=1:numel(Ts)
          ssOtro = dumpEvolvedState(ss, states(:,k));
          rec.e(rec.index) = calculateEnergy(ssOtro);
          rec.index = rec.index+1;
        end
      end
      recordingCallback(rec, ss, Ts, states); %do not forget to call recordingCallback!!!
    end
    
    %specially tailored "playSimulation" since we only have energy
    %information, only energy is passed down to the function
    function playSimulation(rec, Ts, fun)
      es = rec.e(1:(rec.index-1));
      TS = rec.ts(1:(rec.index-1));
      [nevermind, loc] = ismember(TS, Ts);
      clear nevermind;
      loc = loc(loc~=0);
      arrayfun(fun, TS(loc), es(loc));
    end
    
    %find all recorded timestep within the given time interval
    function Ts = findTimestepsBetween(rec, t1, t2)
      allTs = allTimeSteps(rec);
      Ts = allTs((allTs>=t1) & (allTs<=t2));
    end
    
    %find all recorded timesteps in this EnergyRecorder
    function Ts = allTimeSteps(rec)
      Ts = rec.ts(1:(rec.index-1));
    end
      
    %Recorder interface, but useless for MemRecorder
    function rec = flushRecording(rec); end

  end
    
  
 
end