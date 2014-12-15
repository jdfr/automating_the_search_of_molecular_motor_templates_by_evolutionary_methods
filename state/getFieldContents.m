function res = getFieldContents(ss, name)

if iscell(name)
  switch numel(name)
    case 3
      res = ss.(name{1}).(name{2}).(name{3});
    case 2
      res = ss.(name{1}).(name{2});
    otherwise
      error('This field is too much nested (%s)!!!!!', any2str(name));
  end
else
  res = ss.(name);
end
