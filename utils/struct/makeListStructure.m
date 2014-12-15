function STRUCT = makeListStructure(strCell)
STRUCT = reshape([strCell; arra2cellMX(1:numel(strCell))], 1, []);
STRUCT = struct(STRUCT{:});