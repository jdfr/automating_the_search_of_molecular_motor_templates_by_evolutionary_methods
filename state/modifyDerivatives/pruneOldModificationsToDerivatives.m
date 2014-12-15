%delete modifications which are too old
%
% checkEmptyMods: flag to say whether to check if there are no
%                 modifications at all
% synchronizeDynFlags: flag to say whether the dynamical_X flags are to be
%                      updated accordingly to the removal of modifications
% setToSynchronize: indexes in ss.dynVars to selectively reset 
%                   only some flagvars
% emptyMods: flag saying whether there are not any modifications at all
% changedDynFlags: flag saying whether any dynamical_X flag has changed
%
function [ss emptyMods changedDynFlags] = pruneOldModificationsToDerivatives(ss, checkEmptyMods, synchronizeDynFlags, setToSynchronize)
  emptyMods       = checkEmptyMods;
  for k=1:numel(ss.modificationVars)
    mod = ss.modificationVars{k};
    for m=1:numel(ss.dynVars)
      dyn = ss.dynVars{m};
      tooOld = ss.(mod).(dyn).endTs<=ss.t;
      if any(tooOld)
        for n=1:numel(ss.modificationFields) %erase fields
          fld = ss.modificationFields{n};
          ss.(mod).(dyn).(fld)(tooOld,:) = [];
        end
      end
      emptyMods = emptyMods && isempty(ss.(mod).(dyn).endTs);
    end
  end
  
  changedDynFlags = false(numel(setToSynchronize), 1);
  if synchronizeDynFlags
    if ~exist('setToSynchronize', 'var')
      setToSynchronize = 1:numel(ss.dynVars);
    end
    for z=1:numel(setToSynchronize)
      m = setToSynchronize(z);
      dyn = ss.dynVars{m};
      flgvar = ss.doubledFlagVars{m};
      defaultvalue = ss.doubledFlagDefaultVals{m};
      newFlag = repmat(defaultvalue, size(ss.(flgvar)));
      for k=1:numel(ss.modificationVars)
        mod = ss.modificationVars{k};
        newFlag(ss.(mod).(dyn).indexes) = true;
      end
      changedDynFlags(m) = any(xor(newFlag, ss.(flgvar)));
      ss.(flgvar) = newFlag;
    end
  end

end
