%find all recorded timesteps in this EnergyRecorder
function Ts = allTimeSteps(rec)
  Ts = rec.ts(1:(rec.index-1));
end
