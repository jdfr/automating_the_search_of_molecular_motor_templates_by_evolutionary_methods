function ss = reassignField(ss, name, command, newval, selection)
%commands:
% command==0 => simple assignment
% command==1 => add rows at the end
% command==2 => copy data to some selected rows
% command==3 => copy data from some rows to other selected rows
isfun      = isa(newval, 'function_handle');
if (nargin>4) && (command ~= 3)
  selection = selection(selection~=0);
end


if iscell(name)
  switch numel(name)
    case 2
      switch command
        case 0
          if isfun
            ss.(name{1}).(name{2}) = newval(ss.(name{1}).(name{2}));
          else
            ss.(name{1}).(name{2}) = newval;
          end
        case 1
          if isfun
            ss.(name{1}).(name{2}) = [ss.(name{1}).(name{2}); newval(ss.(name{1}).(name{2}))];
          else
            ss.(name{1}).(name{2}) = [ss.(name{1}).(name{2}); newval];
          end
        case 2
          if isfun
            ss.(name{1}).(name{2})(selection,:) = newval(ss.(name{1}).(name{2})(selection,:));
          else
            ss.(name{1}).(name{2})(selection,:) = newval;
          end
        case 3
            ss.(name{1}).(name{2})(selection,:) = ss.(name{1}).(name{2})(newval,:);
        otherwise
          unrecognized;
      end
    case 3
      switch command
        case 0
          if isfun
            ss.(name{1}).(name{2}).(name{3}) = newval(ss.(name{1}).(name{2}).(name{3}));
          else
            ss.(name{1}).(name{2}).(name{3}) = newval;
          end
        case 1
          if isfun
            ss.(name{1}).(name{2}).(name{3}) = [ss.(name{1}).(name{2}).(name{3}); newval(ss.(name{1}).(name{2}).(name{3}))];
          else
            ss.(name{1}).(name{2}).(name{3}) = [ss.(name{1}).(name{2}).(name{3}); newval];
          end
        case 2
          if isfun
            ss.(name{1}).(name{2}).(name{3})(selection,:) = newval(ss.(name{1}).(name{2}).(name{3})(selection,:));
          else
            ss.(name{1}).(name{2}).(name{3})(selection,:) = newval;
          end
        case 3
            ss.(name{1}).(name{2}).(name{3})(selection,:) = ss.(name{1}).(name{2}).(name{3})(newval,:);
        otherwise
          unrecognized;
      end
    otherwise
      error('This field is too much nested (%s)!!!!!', any2str(name));
  end
else
  switch command
    case 0
      if isfun
        ss.(name) = newval(ss.(name));
      else
        ss.(name) = newval;
      end
    case 1
      if isfun
        ss.(name) = [ss.(name); newval(ss.(name))];
      else
        ss.(name) = [ss.(name); newval];
      end
    case 2
      if isfun
        ss.(name)(selection,:) = newval(ss.(name)(selection,:));
      else
        ss.(name)(selection,:) = newval;
      end
    case 3
        ss.(name)(selection,:) = ss.(name)(newval,:);
    otherwise
      unrecognized;
  end
end


function unrecognized
error('Unrecognized command!!!');