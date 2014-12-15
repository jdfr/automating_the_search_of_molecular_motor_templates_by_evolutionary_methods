%for each timestep in the interval [t1,t2], apply the given function
function playSimulationBewteen(rec, t1, t2, fun)
Ts = findTimestepsBetween(rec, t1, t2);
playSimulation(rec, Ts, fun);
