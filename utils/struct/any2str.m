%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%converts almost everything to a string form
function str = any2str(c,forPrintf,nospace)
if nargin<3
  nospace = false;
  if nargin<2
    forPrintf = false;
  end
end
if ischar(c) 
  str = mat2str(c);
  if forPrintf
    str = strrep(str, '\', '\\');
  end
elseif isnumeric(c) || islogical(c)
    str=mat2str(c);
    if nospace
      str(str==' ') = ',';
    end
elseif isstruct(c)
  fns=fieldnames(c);
  if isempty(fns)
    str = 'struct';
  else
    supercell = struct2cell(c);
    strs=cell(size(fns));
    siz = size(c);
    if forPrintf
      apos = '''''';
    else
      apos = '''';
    end
    if nospace
      sep = ',';
    else
      sep = ', ';
    end
    for k=1:numel(fns)
      %this works as long as the struct has at most 2 dimensions
      strs{k} = [apos fns{k} apos sep any2str(reshape(supercell(k,:,:), siz), forPrintf, nospace) sep];
    end
    strs{end} = strs{end}(1:(end-numel(sep)));
    str = ['struct(' strs{:} ')'];
  end
%     if isempty(c)
%         str = 'STRUCT<>';
%     else
%         fns=fieldnames(c);
%         strs=cell(size(fns));
%         for k=1:numel(fns)
%             ministrs = cell(size(c));
%             for m=1:numel(c)
%                 ministrs{m} = [any2str(c(m).(fns{k}), forPrintf, nospace), ', '];
%             end
%             strs{k} = [fns{k} '=' ministrs{:}];
%         end
%         strs{end} = strs{end}(1:(end-2));
%         str = ['STRUCT<' strs{:} '>'];
%     end
elseif iscell(c)
    if isempty(c)
        str = '{}';
    else
        if nospace
          sep = ',';
        else
          sep = ', ';
        end
        cs = cellfun(@(x)[any2str(x, forPrintf, nospace) sep], c, 'uniformoutput', false);
        for k=1:size(cs, 1)
          cs{k,end}(end-numel(sep)+1) = ';';
        end
        cs{end} = cs{end}(1:(end-numel(sep)));
        cs = cs';
        str = ['{' cs{:} '}'];
    end
elseif ishandle(c)
    str = 'HANDLE';
elseif isobject(c)
    str='OBJECT';
elseif isa(c, 'function_handle')
    str=func2str(c);
    if str(1)~='@'
      str = ['@' str];
    end
else
    str='UNKNOWN';
end
