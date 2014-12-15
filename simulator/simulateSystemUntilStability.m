%simulate system in batchs of numel(deltaTs) intermediate timesteps, until
%the difference in energy drops below energyThreshold.
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
%   -deltaTs: Each batch is composed of the timesteps s.t+deltaTs
%   -energyThreshold: after each batch, the absolute difference in energy
%                     between first and last timesteps of the batch is
%                     checked. If it is below energyThreshold, we stop. If
%                     not, we simulate a new batch
%   -verbose: flag to display a line of information after each batch is
%             done
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
%   -customCondition: if not empty, this should be a callback to a boolean function: 
%                          customGoOn = customCondition(ss, numIteration, energyBeforeStep, energyAfterStep, absThresholdExceeded);
%                     If customCondition is not empty, it is called after
%                     each simulation batch, and the resulting customGoOn
%                     flag overrides the normal loop ending test.
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

function [ss rec numIteration abnormalEnding] = simulateSystemUntilStability(ss, rec, deltaTs, energyThreshold, verbose, recordInitialState, customCondition, useKeyboard, varargin); %varargin={clearOldModifications, odefun}
  if ~exist('customCondition', 'var'); customCondition = [];   end;
  if ~exist('useKeyboard',     'var'); useKeyboard     = true; end;
  time = ss.t;
  continuar = true;
  numIteration = 0;
  %initial energy
  eOld = calculateEnergy(ss);
  if useKeyboard
    FS = stoploop('use me to stop this!!!') ;
  end
  while continuar
    numIteration = numIteration+1;
    tSteps = time+deltaTs;
    %solving!
    [ss rec] = simulateSystem(ss, rec, tSteps, recordInitialState, varargin{:});
    recordInitialState = false;
    %the state must not be recorded after
    %calculate difference in energy
    eNew = calculateEnergy(ss);
    %test for energy variation
    exceeded = (abs(eNew-eOld) > energyThreshold);
    %get custom condition, if any
    customGoOn = isempty(customCondition) || customCondition(ss, numIteration, eOld, eNew, exceeded);
    %keyboard termination is TRUE if it IS enabled AND the user IS NOT pressing ctrl-a
    goonKeyboard = (~useKeyboard) || (~FS.Stop()); %any(keyinfo~=[17; 65; 162]);
    %test for going on
    continuar = goonKeyboard && ((isempty(customCondition) && exceeded) || ((~isempty(customCondition)) && customGoOn));
    %display information
    if verbose
        disp(sprintf('Nºiters: %g, t: %g, temp diff: %g, nºpasos: %d; E diff: %g, treshold: %g', numIteration, ss.t, ss.t-time, numel(deltaTs), abs(eNew-eOld), energyThreshold));
    end
    time = ss.t;
    eOld = eNew;
  end
  if useKeyboard
    FS.Clear() ; % Clear up the box
    clear FS ; % this structure has no use anymore    
  end
  
  abnormalEnding = (~goonKeyboard) || ((~isempty(customCondition)) && (~customGoOn));
end
