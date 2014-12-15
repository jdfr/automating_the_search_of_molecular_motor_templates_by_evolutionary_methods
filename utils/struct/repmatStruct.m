%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%makes repmat to all fields in a struct
function str = repmatStruct(str, siz)
fnames = fieldnames(str);
for k=1:numel(fnames)
  str.(fnames{k}) = repmat(str.(fnames{k}), siz);
end
