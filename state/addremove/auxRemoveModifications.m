%recalculate indexes after remove, and make up the dynamics'
%modifications
%THIS IS A PRIVATE HELPER FUNCTION
function ss = auxRemoveModifications(ss, indexesAfterRemove, indexDeleted, specificDynVars)
  %remove related dynamics' modifications, and re-index them
  for k=1:numel(ss.modificationVars)
    modif = ss.modificationVars{k}; %'pertubations', 'enforcements'
    for m=1:numel(specificDynVars)
      dynvar = specificDynVars{m}; %'k', 'r', 'c'
      %only if there are modifications
      if ~isempty(ss.(modif).(dynvar).indexes)
        %determine the modifications to erase
        toErase = ismember(ss.(modif).(dynvar).indexes, indexDeleted);
        %remove the earmarked modifications
        for n=1:numel(ss.modificationFields)
          fld = ss.modificationFields{n};
          ss.(modif).(dynvar).(fld)(toErase,:) = [];
        end
        %re-index the remaining indexes
        ss.(modif).(dynvar).indexes = indexesAfterRemove(ss.(modif).(dynvar).indexes);
      end
    end
  end
end
