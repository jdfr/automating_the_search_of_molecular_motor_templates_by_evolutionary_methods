%helper function for using the state without changing it
function [funout] = makeSimpleIntercalatedFun(fun)

funout = @(ss, rec, numIteration) simpleIntercalated(ss, rec, numIteration, fun);

%this function never aborts nor modifies the state
function [ss rec recordInitialState abortSimulation] = simpleIntercalated(ss, rec, numIteration, fun)
  fun(ss, rec, numIteration);
  recordInitialState = false;
  abortSimulation = false;
