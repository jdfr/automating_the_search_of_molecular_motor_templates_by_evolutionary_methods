%find all recorded timestep within the given time interval
function Ts = findTimestepsBetween(rec, t1, t2)
allTs = allTimeSteps(rec);
Ts = allTs((allTs>=t1) & (allTs<=t2));
