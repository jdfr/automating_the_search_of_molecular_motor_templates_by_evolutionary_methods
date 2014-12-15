%remove springs; also returns an array which contains, in the n-th position,
%the new index for the spring previously using index n. Of course, positions
%for deleted points should not be used
%
% optionally, it accepts a variable amount of spring index matrices, and
% re-index them
function [ss indexesAfterRemove varargout] = removeAtoms(ss, indexDeletedAtoms, deleteSprings, deleteUnconnectedPoints, varargin)
  if ~exist('deleteSprings',            'var'); deleteSprings            = true; end
  if ~exist('deleteUnconnectedPoints',  'var'); deleteUnconnectedPoints  = true; end
  numIndexesBeforeDeletion = size(ss.atomPoints,1);
  %remove atom mark from deleted atoms' points and springs
  for k=1:numel(indexDeletedAtoms)
    atompoints  = ss.atomPoints(indexDeletedAtoms(k), :);
    ss.pointAtoms(atompoints)   = cellfun(@(atomset) atomset(atomset~=indexDeletedAtoms(k)), ss.pointAtoms(atompoints),   'UniformOutput', false);
    atomsprings = ss.atomSprings(indexDeletedAtoms(k), :);
%     ss.springAtoms(atomsprings) = cellfun(@(atomset) atomset(atomset~=indexDeletedAtoms(k)), ss.springAtoms(atomsprings), 'UniformOutput', false);
  end
  if deleteSprings
    affectedSprings = unique(ss.atomSprings(indexDeletedAtoms, :));
  end
  %remove the atom vars associated to these points
  for k=1:numel(ss.atomVars)
    name = ss.atomVars{k};
    ss.(name)(indexDeletedAtoms,:) = [];
  end
  %find new values for surviving indexes
  indexesAfterRemove = findIndexesAfterRemove(numIndexesBeforeDeletion, indexDeletedAtoms);
  %re-index affected vars
  for k=1:numel(ss.atomIndexVars)
    name = ss.atomIndexVars{k};
    ss.(name) = auxReindex(indexesAfterRemove, ss.(name));
  end
  %re-index given vars
  for k=1:numel(varargin)
    varargout{k} = auxReindex(indexesAfterRemove, varargin{k});
  end
  %delete springs
  if deleteSprings
    ss = removeSprings(ss, affectedSprings);
  end
  %delete unconnected points (only if also deleting springs!!!)
  if deleteSprings && deleteUnconnectedPoints
%     affectedPoints  = unique(ss.atomPoints(indexDeletedAtoms, :));
%     affectedPoints(ismember(affectedPoints, ss.springEnds)) = [];
%     ss = removePoints(ss, affectedPoints);
    ss = removeUnconnectedPoints(ss);
  end
end
