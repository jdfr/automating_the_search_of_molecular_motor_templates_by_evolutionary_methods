%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%assigns the fields of str2 to a certain range in the fields of str1
function str1 = assignSubStruct(str1, indexes1, mode1, str2, indexes2, mode2, fnames)

useIndexes1 = ~isempty(indexes1);
if ~useIndexes1
  mode1='n';
end

useIndexes2 = exist('indexes2', 'var') && ~isempty(indexes2);
if ~useIndexes2
  mode2='n';
end

mode = [mode1(1) mode2(1)];

if ~exist('fnames', 'var')
  fnames = fieldnames(str1);
end

switch mode
  case 'nn'
    for k=1:numel(fnames)
      str1.(fnames{k})                            = str2.(fnames{k});
    end
  case 'ne'
    for k=1:numel(fnames)
      str1.(fnames{k})                            = str2.(fnames{k})(indexes2);
    end
  case 'ni'
    for k=1:numel(fnames)
      str1.(fnames{k})                            = str2.(fnames{k})(indexes2(1):indexes2(end));
    end
  case 'en'
    for k=1:numel(fnames)
      str1.(fnames{k})(indexes1)                  = str2.(fnames{k});
    end
  case 'ee'
    for k=1:numel(fnames)
      str1.(fnames{k})(indexes1)                  = str2.(fnames{k})(indexes2);
    end
  case 'ei'
    for k=1:numel(fnames)
      str1.(fnames{k})(indexes1)                  = str2.(fnames{k})(indexes2(1):indexes2(end));
    end
  case 'in'
    for k=1:numel(fnames)
      str1.(fnames{k})(indexes1(1):indexes1(end)) = str2.(fnames{k});
    end
  case 'ie'
    for k=1:numel(fnames)
      str1.(fnames{k})(indexes1(1):indexes1(end)) = str2.(fnames{k})(indexes2);
    end
  case 'ii'
    for k=1:numel(fnames)
      str1.(fnames{k})(indexes1(1):indexes1(end)) = str2.(fnames{k})(indexes2(1):indexes2(end));
    end
  otherwise
    error('mode1 and mode2 must be either ''intensive'' or ''extensive''!!!');
end
