%remove points; also returns an array which contains, in the n-th position,
%the new index for the point previously using index n. Of course, positions
%for deleted points should not be used
%
% optionally, it accepts a variable amount of point index matrices, and
% re-index them
function [ss indexesAfterRemove varargout] = removePoints(ss, indexDeletedPoints, varargin)
  numIndexesBeforeDeletion = size(ss.pos,1);
  %remove the point vars associated to these points
  for k=1:numel(ss.pointVars)
    ss   = reassignField(ss, ss.pointVars{k}, 2, [], indexDeletedPoints);
  end
  %find new values for surviving indexes
  indexesAfterRemove = findIndexesAfterRemove(numIndexesBeforeDeletion, indexDeletedPoints);
%   %remove the programmed modifications, and re-index the remaining ones
%   ss = auxRemoveModifications(ss, indexesAfterRemove, indexDeletedPoints, ss.pointDynVars);
  %re-index affected vars
  reindexer = myauxReindex(indexesAfterRemove);
  for k=1:numel(ss.pointIndexVars)
    if ~isempty(ss.(ss.pointIndexVars{k}))
      ss   = reassignField(ss, ss.pointIndexVars{k}, 0, reindexer);
      %ss   = reassignField(ss, ss.pointIndexVars{k}, 2, reindexer, indexesAfterRemove);
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