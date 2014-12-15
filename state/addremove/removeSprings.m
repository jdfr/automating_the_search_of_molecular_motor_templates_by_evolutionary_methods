%remove springs; also returns an array which contains, in the n-th position,
%the new index for the spring previously using index n. Of course, positions
%for deleted points should not be used
%
% optionally, it accepts a variable amount of spring index matrices, and
% re-index them
function [ss indexesAfterRemove varargout] = removeSprings(ss, indexDeletedSprings, varargin)
  numIndexesBeforeDeletion = size(ss.springEnds,1);
  %remove the spring vars associated to these points
  for k=1:numel(ss.springVars)
    ss = reassignField(ss, ss.springVars{k}, 2, [], indexDeletedSprings);
  end
  %find new values for surviving indexes
  indexesAfterRemove = findIndexesAfterRemove(numIndexesBeforeDeletion, indexDeletedSprings);
%   %remove the programmed modifications, and re-index the remaining ones
%   ss = auxRemoveModifications(ss, indexesAfterRemove, indexDeletedSprings, ss.springDynVars);
  %re-index affected vars
  reindexer = myauxReindex(indexesAfterRemove);
  for k=1:numel(ss.springIndexVars)
    if ~isempty(ss.(ss.springIndexVars{k}))
      ss = reassignField(ss, ss.springIndexVars{k}, 0, reindexer);
      %ss = reassignField(ss, ss.springIndexVars{k}, 2, reindexer, indexesAfterRemove);
    end
  end
  %re-index given vars
  for k=1:numel(varargin)
    varargout{k} = auxReindex(indexesAfterRemove, varargin{k});
  end
end

function fun = myauxReindex(indexes)
fun =  @(x)auxReindex(indexes, x);
end