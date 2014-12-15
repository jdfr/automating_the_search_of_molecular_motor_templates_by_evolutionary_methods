function result = array2cell(array)
result = cell(size(array));
for k=1:numel(array)
  result{k} = array(k);
end;
%result = arrayfun(@(x)x, array, 'uniformoutput', false);