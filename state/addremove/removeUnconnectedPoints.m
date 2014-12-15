function ss = removeUnconnectedPoints(ss)

toDelete = true(size(ss.pos, 1), 1);
toDelete(ss.springEnds) = false;
toDelete = find(toDelete);

if ~isempty(toDelete)
  ss = removePoints(ss, toDelete);
end