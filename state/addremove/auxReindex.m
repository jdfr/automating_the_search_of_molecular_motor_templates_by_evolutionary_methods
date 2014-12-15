%recursive re-indexer, for nested sets of indexes
function indexes = auxReindex(newIndexes, indexes)
if iscell(indexes) % if the var is a cell array, fix each cell
  indexes = cellfun(@(x) auxReindex(newIndexes, x), indexes, 'UniformOutput', false);
else
  indexes = newIndexes(indexes);
end
