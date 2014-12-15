%now, modify basic behaviour.
function [dpos dvel dm dk dr dc du] = nonNaturalDynamics(ss, t, dpos, dvel, dm, dk, dr, dc, du)
  %now, modify basic behaviour. For each dynamic var:
  %  -first, we perform the perturbations 
  %  -second, we let customizable modifications to the system's dynamics
  %  -finally, we perform the enforcements to dynamics
  %of course, this code can be condensed into a loop by addressing
  %properties dynamically. However, it severely slows down the
  %function, so we must unfold the loop. Also, we may customize
  %getModification into several versions, each one for the specific
  %field, and save time. However, we would need 12 customized versions
  %of a fairly complex code, quite unmanageable! So, it pays to keep
  %getModification as generic as possible. Just an example: to save time,
  %the first check in getModification has been unfolded to not call i if
  %there is no need. As it can be seen, the code is embarrasingly
  %repetitive and cumbersome. To unfold all the function is too tricky
  if ~isempty(ss.perturbations.pos.startTs); dpos = getModification(ss, t, 'perturbations', 1, dpos); end
  if ~isempty(ss.perturbations.vel.startTs); dvel = getModification(ss, t, 'perturbations', 2, dvel); end
  if ~isempty(ss.modify_p);                  [dpos, dvel] = ss.modify_p(ss, dpos, dvel, t);   end
  if ~isempty(ss.enforcements.pos.startTs);  dpos = getModification(ss, t, 'enforcements',  1, dpos); end
  if ~isempty(ss.enforcements.vel.startTs);  dvel = getModification(ss, t, 'enforcements',  2, dvel); end
  if ~isempty(dm);
    if ~isempty(ss.perturbations.m.startTs); dm = getModification(ss, t, 'perturbations', 3, dm); end
    if ~isempty(ss.modify_m);                dm = ss.modify_m(ss, dm, t, 'm'); end
    if ~isempty(ss.enforcements.m.startTs);  dm = getModification(ss, t, 'enforcements',  3, dm); end
  end
  if ~isempty(dk);
    if ~isempty(ss.perturbations.k.startTs); dk = getModification(ss, t, 'perturbations', 4, dk); end
    if ~isempty(ss.modify_k);                dk = ss.modify_k(ss, dk, t, 'k'); end
    if ~isempty(ss.enforcements.k.startTs);  dk = getModification(ss, t, 'enforcements',  4, dk); end
  end
  if ~isempty(dr);
    if ~isempty(ss.perturbations.r.startTs); dr = getModification(ss, t, 'perturbations', 5, dr); end
    if ~isempty(ss.modify_r);                dr = ss.modify_r(ss, dr, t, 'r'); end
    if ~isempty(ss.enforcements.r.startTs);  dr = getModification(ss, t, 'enforcements',  5, dr); end
  end
  if ~isempty(dc);
    if ~isempty(ss.perturbations.c.startTs); dc = getModification(ss, t, 'perturbations', 6, dc); end
    if ~isempty(ss.modify_c);                dc = ss.modify_c(ss, dc, t, 'c'); end
    if ~isempty(ss.enforcements.c.startTs);  dc = getModification(ss, t, 'enforcements',  6, dc); end
  end
  if ~isempty(du);
    if ~isempty(ss.perturbations.u.startTs); du = getModification(ss, t, 'perturbations', 7, du); end
    if ~isempty(ss.modify_u);                du = ss.modify_u(ss, du, t, 'u'); end
    if ~isempty(ss.enforcements.u.startTs);  du = getModification(ss, t, 'enforcements',  7, du); end
  end


%this function calculates the vars which are modified at the present T.
%This is very tricky, since modified vars are characterized by their
%indexes over all vars, while dynamical vars have their own indexes. To
%join both lists of indexes, we use the cumsum of the logical array
%which determines the indexes of dynamical vars.
%This function is also tricky because it is very, very generic: it
%handles both perturbations (which are added to the vars) and
%enforcements (which replace them), for every possible var (pos, vel,
%m, k, r, c).
%Arguments:
%  ss: BasicSystem object
%  t: current time
%  modif: string: either 'perturbations' or 'enforcements'
%  k_dynvar: index of the var to be modified
%  derivative: derivatives of the vars: they are the target to
%              modification
function derivative = getModification(ss, t, modif, k_dynvar, derivative)
  %string of the var: 'pos', 'vel', 'm', 'k', 'r', 'c'
  dynvar       = ss.dynVars{k_dynvar};
  %if there are not modifications, let's exit
  %THIS LINE IS UNFOLDED IN THE CALLER FUNCTION TO SAVE TIME; DYNAMICAL
  %FIELDS ARE RELATIVELY EXPENSIVE
%  if isempty(ss.(modif).(dynvar).startTs); return; end
  %determine the vars modifications within a valid time window
  withinTime = (ss.(modif).(dynvar).startTs<=t) & (ss.(modif).(dynvar).endTs>=t);
  %if there are not modifications at this T, let's exit
  if ~any(withinTime); return; end %if isempty(indexesVar); return; end
  %determine the indexes of vars within valid time window  
  indexesVar = ss.(modif).(dynvar).indexes(withinTime);
  %string identifying the var's logical array, which address dynamical
  %vars to manipulate them
  dynamical_x  = ss.doubledFlagVars{k_dynvar};
  %find modified vars which are not dynamical
  mod_butnot_dyn = ~ss.(dynamical_x)(indexesVar);
  %prune these vars from our selection
  indexesVar(mod_butnot_dyn) = [];
  %if there are not modifications to dynamical vars, let's exit
  if isempty(indexesVar); return; end
  %string identifying cumsum of vars' logical array. It is useful to
  %find indexes of dynamical vars
  cumsumvar = ss.doubledCumsumFlagVars{k_dynvar};
  %get indexes in dynamical array by re-indexing through the cumsum-ed version of dynamical_X 
  indexes_in_dyn_array = ss.(cumsumvar)(indexesVar);
  %if no var is both dynamical and modified, let's exit
  if isempty(indexes_in_dyn_array); return; end
  %find the modifications
  modifications = ss.(modif).(dynvar).values(withinTime,:);
  modifications(mod_butnot_dyn,:) = [];
  %apply the modifications
  switch modif(1)
    case 'p' %perturbations: 
      if ss.severalPerts.(dynvar) %add possibly more than one perturbation to each modified var 
        if size(derivative,2)==1 %if the var is scalar, it is easy
          derivative = derivative+sum(sparse(indexes_in_dyn_array, 1:numel(indexes_in_dyn_array), modifications, size(derivative,1), numel(indexes_in_dyn_array)), 2);
        else %if the var is vectorial, it is a bit more tricky and computationally expensive. It is just like adding spring forces in points
          derivative = derivative+(sparse(indexes_in_dyn_array, 1:numel(indexes_in_dyn_array), 1, size(derivative,1), numel(indexes_in_dyn_array))*modifications);
        end
      else %only one perturbation per modified var: far easier
        derivative(indexes_in_dyn_array,:) = derivative(indexes_in_dyn_array,:)+modifications;
      end
    case 'e' %enforcements: replace the var. If there are more than one enforcement, arbitrarily select one of them
      derivative(indexes_in_dyn_array,:) = modifications;
  end
