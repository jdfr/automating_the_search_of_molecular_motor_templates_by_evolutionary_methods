function [ss varargout] = copySprings(ss, idxs)

numps = numel(idxs);
func  = copyFun(idxs);
%output indexes
if nargout>1; varargout{1} = size(ss.k,1)+(1:numps)'; end %indexes

names = ss.springVars;
for k=1:numel(names)
  ss = reassignField(ss, names{k}, 1, func);
end

function fun = copyFun(idx)
fun = @(values) values(idx,:);
