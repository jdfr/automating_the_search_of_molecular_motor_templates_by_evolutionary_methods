function [comparison, compdata, toCompare, data] = checkSavedExamples(nsims, path1, path2)

data = compileExamples('data');

if isempty(nsims)
  nsims = (1:numel(data))';
end

if isempty(path1)
  path1 = '%s.mat';
end

if isempty(path2)
  path2 = 'Z_%s.mat';
end

comparison = cell(numel(nsims),1);

toCompare = {'fitness', 'offsetAbs', 'offsetRel', 'atpChg1', 'atpChg2', 'cpuTime'};

for k=1:numel(nsims)
  nsim  = nsims(k);
  path  = data{nsim}{3};
  ninds = numel(data{nsim}{5});
  sps   = find(path==filesep); 
  path  = path(sps(end-1)+1:sps(end)-1);
  
  pth1  = sprintf(path1, path);
  fprintf('in #%03d of #%03d, loading %s...\n', k, numel(nsims), pth1);
  name1 = who('-file', pth1);
  name1 = name1{1};
  load(pth1);
  res1  = eval(name1);
  clear(name1);
  
  pth2  = sprintf(path2, path);
  fprintf('in #%03d of #%03d, loading %s...\n', k, numel(nsims), pth2);
  name2 = who('-file', pth2);
  name2 = name2{1};
  load(pth2);
  res2  = eval(name2);
  clear(name2);
  
  fprintf('in #%03d of #%03d, comparing %s and %s...\n', k, numel(nsims), pth1, pth2);
  comparison{k} = [repmat({path}, ninds, 1) array2cellMX(res1.inds) array2cellMX(res2.inds) cell(ninds, 2*numel(toCompare))];
  for z=1:ninds
    zz = 4;
    for q=1:numel(toCompare)
      comparison{k}{z, zz} = res1.structs{z}.stats.(toCompare{q});
      zz = zz+1;
      comparison{k}{z, zz} = res2.structs{z}.stats.(toCompare{q});
      zz = zz+1;
    end
  end
  
  clear res1 res2;
end

comparison = vertcat(comparison{:});

compdata = cell2arrayMX(comparison(:,2:end));

compdata = compdata(:,1:2:end)-compdata(:,2:2:end);

