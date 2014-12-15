function [comparison, compdata, toCompare, data] = checkSavedExamplesWithPoph(conf, nsims, path1)

data = compileExamples('data');

if isempty(nsims)
  nsims = (1:numel(data))';
end

if isempty(path1)
  path1 = 'Z_%s.mat';
end

comparison = cell(numel(nsims),1);

toCompare = {'fitness', 'offsetAbs', 'offsetRel', 'atpChg1', 'atpChg2', 'nsprings1', 'spectralGap', 'cpuTime'};

for k=1:numel(nsims)
  nsim  = nsims(k);
  path  = data{nsim}{3};
  poph  = conf.fun.loadSim(path);
  ninds = numel(data{nsim}{5});
  sps   = find(path==filesep); 
  path  = path(sps(end-1)+1:sps(end)-1);

  esfield2   = isfield(poph, toCompare);

  pth1  = sprintf(path1, path);
  fprintf('in #%03d of #%03d, loading %s...\n', k, numel(nsims), pth1);
  name1 = who('-file', pth1);
  name1 = name1{1};
  load(pth1);
  res1  = eval(name1);
  clear(name1);
  
  fprintf('in #%03d of #%03d, comparing %s with poph records...\n', k, numel(nsims), pth1);
  comparison{k} = [repmat({path}, ninds, 1) array2cellMX(res1.inds) cell(ninds, 2*numel(toCompare))];
  for z=1:ninds
    zz = 3;
    for q=1:numel(toCompare)
      if isfield(res1.structs{z}.stats, toCompare)
        comparison{k}{z, zz} = res1.structs{z}.stats.(toCompare{q});
      else
        comparison{k}{z, zz} = eps(0);
      end
      zz = zz+1;
      if esfield2(q)
        comparison{k}{z, zz} = poph.(toCompare{q})(res1.inds(z));
      else
        comparison{k}{z, zz} = eps(0);
      end
      zz = zz+1;
    end
  end
  
  clear res1 res2;
end

comparison = vertcat(comparison{:});

compdata = cell2arrayMX(comparison(:,3:end));

compdata = compdata(:,1:2:end)-compdata(:,2:2:end);

