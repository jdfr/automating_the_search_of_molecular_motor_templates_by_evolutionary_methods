%return the system for a given T
function ss = getTimeStep(rec, T)
  ss = [];
  playSimulation(rec, T, @assignSS);
  function assignSS(t, system); ss = system; end
end
