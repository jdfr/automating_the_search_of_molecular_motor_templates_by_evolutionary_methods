%get the BasicSystem's fields which are lists of other fields. If
%modeStruct is true, the result is delivered as a struct. If is false, the
%results are delivered as a cell array of the form
%    {'fieldSetName1', fieldSetName1, 'fieldSetName2', fieldSetName2, ...}
function out = getBasicSystemFieldSets(bs, modeStruct)
if ~exist('bs', 'var')
  n = zeros(0,1);
  bs = BasicSystem(1);
end
if ~exist('modeStruct', 'var')
  modeStruct = false;
end
nf = numel(bs.fieldSets);
contents = cell(1, nf);
for k=1:nf
  contents{k} = bs.(bs.fieldSets{k});
end
out = reshape([bs.fieldSets; cellfun(@(x){x}, contents, 'Uniformoutput', false)], [], 1);
if modeStruct
  out = struct(out{:});
end
