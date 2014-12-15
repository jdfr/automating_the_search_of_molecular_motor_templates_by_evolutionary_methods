%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = DAConvertGenome(in, geneTables, genSepChar)
% Build a genome string from a list of genes, or vice-versa. The second
% argument

if nargin<3
  genSepChar = '';
  %genSepChar = sprintf('\n');
else
  %genSepChar = sprintf(genSepChar);
end

if iscell(in)
  geneAbrevs = fieldnames(geneTables.abrev2num);
  %out = cellfun(@(x)gen2str(x, geneAbrevs, genSepChar), reshape(in, 1, []), 'UniformOutput', false);
  out = cell(1, numel(in));
  for k=1:numel(out)
    out{k} = gen2str(in{k}, geneAbrevs, genSepChar);
  end
  out = horzcat(out{:});
elseif ischar(in)
  genNameExpr = '[a-zA-Z][a-zA-Z0-9_]*';
  genExpr = ['{[\n\t ]*' genNameExpr '[^}]*}'];
  [s e] = regexp(in, genExpr);
  numgenes = size(s,2);
  out = repmat({struct('type', {[]}, 'args', {[]})}, size(s));
  for k=1:numgenes
    genstr = in(s(k):e(k));
    [sn en] = regexp(genstr, genNameExpr, 'once');
    name = genstr(sn:en);
    out{k}.type = geneTables.abrev2num.(name);
    args = genstr(en+1:end-1);
    %out{k}.args = cellfun(@processStr, regexp(args, '[^, ]*', 'match'), 'UniformOutput', false);
    args = regexp(args, '[^, ]*', 'match');
    for z=1:numel(args)
      args{z} = processStr(args{z});
    end
    out{k}.args = args;
  end
else
  error('No pillo lo que dices, pisha');
end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function str = gen2str(gen, geneAbrevs, genSepChar)
    name = geneAbrevs{gen.type};
    if isfield(gen, 'args')
      args = gen.args;
      for m=1:numel(args)
        args{m} = arg2str(args{m});
      end
%       args = cellfun(@arg2str, reshape(gen.args, 1, []), 'UniformOutput', false);
      args = horzcat(args{:});
    else
      args = '';
    end
    if ~isempty(args)
      str = ['{' name ',' args(1:end-1) '}' genSepChar];
    else
      str = ['{' name '}' genSepChar];
    end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = processStr(str)
[num ok] = str2num(str); %#ok<ST2NM>
if ok
  out = num;
else
  out = str;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = arg2str(arg)
if isempty(arg)
  str = '[],';
else
  str = [num2str(arg, '%20.20g') ',']; %also works fine for strings
end
end
