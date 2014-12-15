function indexesAfterRemove = findIndexesAfterRemove(numIndexesBeforeDeletion, indexDeleted)
  %recalculate indexes after remove
  remaining = true(numIndexesBeforeDeletion,1);
  remaining(indexDeleted) = false;
  indexesAfterRemove = cumsum(remaining);
  indexesAfterRemove(indexDeleted) = 0;
