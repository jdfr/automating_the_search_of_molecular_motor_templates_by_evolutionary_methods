function [comparison, compdata, toCompare] = checkSavedExamplesWithResultsRec(basedir, conf, path1, firstCall)

if isempty(path1)
  path1 = 'P_%s.mat';
end

if ~exist('firstCall', 'var')
  firstCall = true;
end

comparison = {};
compdata   = {};
k=1;

toCompare = {'fitness', 'offsetAbs', 'offsetRel', 'atpChg1', 'atpChg2', 'nsprings1', 'spectralGap', 'cpuTime'};



if ischar(basedir)
  nc = fprintf('Processing folder %s...', basedir);

  ds = dir(basedir);

  [names dirs] = arrayfun(@(x)deal(x.name, x.isdir), ds, 'uniformoutput', false);

  dirs = find(cell2mat(dirs));

  fprintf(repmat('\b', 1, nc));

  if all(ismember({'poblacion.txt', 'arbol.txt', 'timestats.txt'}, names)) %&& (~any(strcmp(savename, names)))
    path = basedir;
    seps = find(path==filesep);
    path = path(seps(end-1)+1:seps(end)-1);
    
    pth1  = sprintf(path1, path);
    
    if exist(pth1, 'file')
    
      name1 = who('-file', pth1);
      name1 = name1{1};
      load(pth1);
      res1  = eval(name1);
      clear(name1);
      
      ninds = numel(res1.structs);
      esfield2   = isfield(res1.inSim, toCompare);
      
      params     = res1.mode;
      if ischar(params)
        load([basedir filesep '..' filesep 'estado.mat']);
        %poph = conf.fun.loadSim(basedir);
        %params = poph.params;
      end
      useFO       = params.evparams.useFirstOrder;
      useHOH      = strcmpi(params.evparams.walker.d3.rotationMode, 'handoverhand');
      fitnessOK   = strcmpi(params.evparams.fitnessCalc, 'mode3');
      comparison{k} = [repmat({path useFO useHOH fitnessOK}, ninds, 1) array2cellMX(res1.inds) cell(ninds, 2*numel(toCompare))];
      for z=1:ninds
        zz = 6;
        for q=1:numel(toCompare)
          if isfield(res1.structs{z}.stats, toCompare)
            comparison{k}{z, zz} = res1.structs{z}.stats.(toCompare{q});
          else
            comparison{k}{z, zz} = eps(0);
          end
          zz = zz+1;
          if esfield2(q)
            comparison{k}{z, zz} = res1.inSim.(toCompare{q})(z);
          else
            comparison{k}{z, zz} = eps(0);
          end
          zz = zz+1;
        end
      end

      clear res1 res2;
      
    end
  else
    for k=1:numel(dirs)
      name = names{dirs(k)};
      if ~any(strcmp(name, {'.', '..'}))
        [comparison2] = checkSavedExamplesWithResultsRec([basedir filesep name], conf, path1, false);
        comparison = [comparison; comparison2];
      end
    end
  end

elseif iscell(basedir)
  for k=1:numel(basedir)
    [comparison2] = checkSavedExamplesWithResultsRec(basedir{k}, conf, path1, false);
    comparison = [comparison; comparison2];
  end
else
  error('basedir not understood!!!!');
end


if firstCall
  comparison = vertcat(comparison{:});

  compdata   = cell2arrayMX(comparison(:,6:end));

  compdata   = compdata(:,1:2:end)-compdata(:,2:2:end);
end


