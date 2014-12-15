function ss2 = copyFields(ss1, ss2, names)
for k=1:numel(names)
  name = names{k};
  if iscell(name)
    switch numel(name)
      case 2
        ss2.(name{1}).(name{2})           = ss1.(name{1}).(name{2});
      case 3
        ss2.(name{1}).(name{2}).(name{3}) = ss1.(name{1}).(name{2}).(name{3});
      otherwise
        error('This field is too much nested (%s)!!!!!', any2str(name));
    end
  else
    ss2.(name) = ss1.(name);
  end
end