%delete modifications (both 'perturbations' and 'enforcements' related to a
%set of derivatives for a set of dynamical variables (selectedIndexes).
%
% checkEmptyMods: flag to say whether to check if there are no
%                 modifications in the selected categories
% synchronizeDynFlags: flag to say whether the dynamical_X flags are to be
%                      updated accordingly to the removal of modifications
% emptyMods: flag saying whether there are not any modifications in the
%            selected categories
% changedDynFlags: flag saying whether any dynamical_X flag has changed
%
%BEWARE: the argument selectedDynVars should only contains dynvars for
%points or for springs, but not for both of them.
%   Examples:
%     {'pos', 'm'}: OK: only dynvars for points
%     {'k', 'r'}:   OK: only dynvars for springs
%     {'m', 'r'}:   NOT OK!!!!!
function [ss emptyMods changedDynFlags] = pruneModificationsToDerivatives(ss, checkEmptyMods, synchronizeDynFlags, selectedDynVars, selectedIndexes)
  emptyMods       = checkEmptyMods;
  changedDynFlags = false;
  for k=1:numel(ss.modificationVars)
    mod = ss.modificationVars{k};
    for m=1:numel(selectedDynVars)
      dyn = selectedDynVars{m};
      removeMods = ismember(ss.(mod).(dyn).indexes, selectedIndexes);
      if any(removeMods)
        for n=1:numel(ss.modificationFields)
          fld = ss.modificationFields{n};
          ss.(mod).(dyn).(fld)(removeMods,:) = [];
        end
        if synchronizeDynFlags
          flgvar = ss.doubledFlagVars{m};
          newFlag = false(size(ss.(flgvar)));
          newFlag(ss.(mod).(dyn).indexes) = true;
          changedDynFlags = any(xor(newFlag, ss.(flgvar)));
          ss.(flgvar) = newFlag;
        end
      end
      emptyMods = emptyMods && isempty(ss.(mod).(dyn).endTs);
    end
  end
end
