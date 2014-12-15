%simulate the system in batchs, each one having n consecutive timesteps
%from tSteps (n=numStepsIntercalate). The last batch may have a lesser
%amount of timesteps, if the amount of timestep is not divisible by the
%batch size. After each batch, a function callback is executed. This
%function must have the folowing interface:
%        [ss rec recordNextInitialState abortSimulation] = funIntercalate(ss, rec, numIteration) 
%Where:
%   -recordNextInitialState: is passed over to the next simulation batch
%                            the recordInitialState flag.
%   -abortSimulation: if true, the remaining batches are aborted
%Parameters:
%   -ss: system structure. It is critical that the derived fields are
%        consistent. To guarantee it:
%          *If the number of dimensions change:
%              +run makeSimulationAmounts
%          *If points are added or removed:
%              +run makeSimulationAmounts
%              +run makeDynamicVars
%          *If springs are added or removed:
%              +run makeSimulationAmounts
%              +run makeDynamicVars
%              +run makeSpringEndsMatrix
%          *If any flag in dinamical_X vectors flag is switched:
%              +run makeSimulationAmounts
%              +run makeDynamicVars
%          *If springs' end points are changed:
%              +run makeSimulationAmounts
%              +run makeSpringEndsMatrix
%        The function prepareSimulationStructure is provided as a
%        convenience for running these three functions together
%   -rec: recording system
%   -tSteps: the time steps to be recorded. BEWARE: they must be
%            consecutive and greater than the initial ss.t
%   -numStepsIntercalate: number of timesteps to be simulated in a single
%                         batch
%   -funIntercalate: the callback to be executed at the end of each batch
%   -recordInitialState: if the flag is true, the whole initial state is
%                        recorded. This **MUST** be true if the ss
%                        structure has changed since the last simulation
%                        period. This holds if any of the following fields
%                        is changed (modified, shrinked or expanded):
%                          *t, u, pos, vel, m, springEnds, k, r, c
%                          *dynamical_p, dynamical_m, dynamical_r, dynamical_c 
%                          *perturbations, severalPerts, enforcements
%                          *odeOptions
%                          *ANY ADDITIONAL FIELDS ADDED TO THE STRUCTURE
%                           WHICH ARE USED IN ANY modifyXXX CALLBACK
%   -useKeyboard: if true, the user can abort by pressing ctrl-a just after
%                 an iteration
%   -clearOldModifications: if true, this flag tells the simulator to wipe
%                           out dynamic's modification which are older than
%                           the simulated system ss.t, after the simulation
%   -odefun: ODE solver to be used (the odeOptions is a field of ss)
%
%Returns:
%   -ss:  the updated system structure
%   -rec: the updated recording system 
%   -numIteration: the number of performed iterations 
%   -abnormalEnding: true iff the loop is terminated before stability is
%                    achieved (if the termination condition is not that the
%                    system's energy has changed below energyThreshold)
function [ss rec recordInitialState numIteration goonKeyboard] = simulateSystemIntercalatedCoarse(ss, rec, deltaTs, funIntercalate, recordInitialState, useKeyboard, varargin); %varargin = {clearOldModifications, odefun}
  if numel(deltaTs)==0; return; end
  if ~exist('useKeyboard', 'var'); useKeyboard = true; end;
  %batch's timesteps' indexes
  continuar = true;
  numIteration = 0;

  if useKeyboard
    FS = stoploop('use me to stop this!!!') ;
  end
  while continuar
    numIteration = numIteration+1;
    %solving!
    tSteps   = ss.t+deltaTs;
    [ss rec] = simulateSystem(ss, rec, tSteps, recordInitialState, varargin{:});
    %intercalated callback
    [ss rec recordInitialState abortSimulation] = funIntercalate(ss, rec, numIteration);
    %keyboard termination is TRUE if it IS enabled AND the user IS NOT pressing ctrl-a
    goonKeyboard = (~useKeyboard) || (~FS.Stop()); %any(keyinfo~=[17; 65; 162]);
    continuar = (~abortSimulation) && goonKeyboard;
  end
  if useKeyboard
    FS.Clear() ; % Clear up the box
    clear FS ; % this structure has no use anymore    
  end
end

