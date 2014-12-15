%find all recorded timesteps in this MemRecorder
function Ts = allTimeSteps(rec)
%Ts = sort(unique([rec.dynTs(1:(rec.dynIndex-1)); rec.notdynTs(1:(rec.notdynIndex-1))]));
%elements are presumed to be pre-sorted
Ts = mergesorted(rec.dynTs(1:(rec.dynIndex-1)), rec.notdynTs(1:(rec.notdynIndex-1)));
%remove repeated elements
Ts(diff(Ts)==0) = [];
