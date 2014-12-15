function STRUCT = makeEmptyStructure(strCell)

aux = [reshape(strCell, 1, []); repmat({[]}, 1, numel(strCell))];
STRUCT = struct(aux{:});
