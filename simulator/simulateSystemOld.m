%simulate a sistem until the time t=tSteps(end), recording the system's
%state at each time in tSteps. The recording is done through the 'rec'
%closure
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
%   -clearOldModifications: if true, this flag tells the simulator to wipe
%                           out dynamic's modification which are older than
%                           the simulated system ss.t, after the simulation
%   -odefun: ODE solver to be used (the odeOptions is a field of ss)
%
%Returns:
%   -ss:  the updated system structure
%   -rec: the updated recording system 
function [ss rec] = simulateSystem(ss, rec, tSteps, recordInitialState, clearOldModifications, odefun)
  %defalt values for optional arguments
  if nargin<5; clearOldModifications = true; end
  if nargin<6; odefun = @ode23; end
  
  %make the vector state for passing it over to the ODE solver.
  evolvedState = makeEvolvedState(ss);

  %if needed, record the given state
  if recordInitialState && ~isempty(rec)
    rec = recordAllState(rec, ss);
  end

  %set up the tspan argument for the ODE solver
  if numel(tSteps)==1
    %a dummy middle step for not being given all intermediate integration steps
    tspan = [ss.t, (ss.t+tSteps)/2, tSteps];
  else
    tspan = [ss.t, tSteps];
  end

  %solving!
  [Ts, results] = odefun(@systemDynamics, tspan, evolvedState, ss.odeOptions, ss);
  
  %get the state back
  ss.t = Ts(end);
  ss = dumpEvolvedState(ss, results(end,:)');
  
  %record the timesteps
  if ~isempty(rec)
    results = results';
    if numel(tSteps)==1
      rec = recordDynState(rec, ss, Ts(end), results(:,end));
    else
      rec = recordDynState(rec, ss, Ts(2:end), results(:,2:end));
    end
  end
  
  %prune old modifications
  if clearOldModifications
    ss = pruneOldModificationsToDerivatives(ss, false, false);
  end
end
