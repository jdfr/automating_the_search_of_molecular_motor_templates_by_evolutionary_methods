%this recorder does nothing when recording but calling the callbacks
classdef DummyRecorder < Recorder
  
  methods
    %constructor
    function rec = DummyRecorder(ss)
      rec = rec@Recorder(ss);
      rec.TYPE = 'DymmyRecorder';
    end
    
    function rec = recordAllState(rec, ss)
      recordingCallback(rec, ss); %do not forget to call recordingCallback!!!
    end
    
    function rec = recordDynState(rec, ss, Ts, states)
      recordingCallback(rec, ss, Ts, states); %do not forget to call recordingCallback!!!
    end
    
    function playSimulation(rec, Ts, fun)
      error('playSimulation is not supported by DummyRecorder');
    end
    
    function Ts = findTimestepsBetween(rec, t1, t2)
      Ts = [];
      error('findTimestepsBetween is not supported by DummyRecorder');
    end
    
    function Ts = allTimeSteps(rec)
      Ts = [];
      error('allTimeSteps is not supported by DummyRecorder');
    end
      
    function rec = flushRecording(rec); end

  end
    
  
end