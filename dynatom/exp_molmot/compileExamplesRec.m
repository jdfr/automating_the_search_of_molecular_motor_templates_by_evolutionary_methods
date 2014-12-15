function [varargout] = compileExamplesRec(conf, fld, basedir, save, varargin)

if isempty(conf)
  conf = molMot3DConf(false);
end
if isempty(fld)
  fld = 'offsetAbs';
end
dosave = ~isempty(save) && ischar(save);

results = [];

if ischar(basedir)
  nc = fprintf('Processing folder %s...', basedir);

  ds = dir(basedir);

  [names dirs] = arrayfun(@(x)deal(x.name, x.isdir), ds, 'uniformoutput', false);

  dirs = find(cell2mat(dirs));

  fprintf(repmat('\b', 1, nc));

  if all(ismember({'poblacion.txt', 'arbol.txt', 'timestats.txt'}, names)) %&& (~any(strcmp(savename, names)))
    pth = basedir;
    seps = find(pth==filesep);
    pth = pth(seps(end-1)+1:seps(end)-1);
    if ~ ( dosave && exist([save pth '.mat'], 'file') )
      fprintf('Processing simulation %s...\n', basedir);
      fprintf('   Loading simulation...\n');
      poph                = conf.fun.loadSim(basedir);
      [idx idx]           = max(poph.(fld));
      fprintf('   Processing...\n');
      results             = struct('path', basedir, 'mode', poph.params, 'inds', idx, 'inSim', struct, 'structs', {{[]}});
      pophflds = fieldnames(poph);
      for k=1:poph.numFieldsPop
        results.inSim.(pophflds{k}) = poph.(pophflds{k})(idx);
      end
      results.structs{1}  = conf.fun.replayGenome(poph, 'poph', idx, varargin{:}); 
      if dosave
        fprintf('   Saving...\n');
        saveR(save, pth, results);
      end
    end
  else
    for k=1:numel(dirs)
      name = names{dirs(k)};
      if ~any(strcmp(name, {'.', '..'}))
        if nargout>0
          results = [results; compileExamplesRec(conf, fld, [basedir filesep name], save, varargin{:})]; %#ok<AGROW>
        else
          compileExamplesRec(conf, fld, [basedir filesep name], save, varargin{:})
        end
      end
    end
  end

elseif iscell(basedir)
  for k=1:numel(basedir)
    if nargout>0
      results = [results; compileExamplesRec(conf, fld, basedir{k}, save, varargin{:})]; %#ok<AGROW>
    else
      compileExamplesRec(conf, fld, basedir{k}, save, varargin{:})
    end
  end
else
  error('basedir not understood!!!!');
end

if nargout>0
  varargout{1} = results;
end

function saveR(saveStr, pth, results) %#ok<INUSD>
save([saveStr pth], 'results');
