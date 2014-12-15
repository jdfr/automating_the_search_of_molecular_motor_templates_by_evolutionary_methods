function [allcpus simulationdirs infos numevs medevtime] = getAllCPUTime(basedir, simulationdirs)

if ~exist('basedir', 'var')
  basedir = ['resultados' filesep 'molmot' filesep];
end
if ~exist('simulationdirs', 'var')
  d = dir(basedir);
  toUse = false(size(d));
  for k=1:numel(d)
    toUse(k) = d(k).isdir && d(k).name(1)=='Y';
    if toUse(k)
      nm = [basedir d(k).name filesep];
      md = dir(nm);
      nm2 = [];
      for z=1:numel(md)
        if md(z).isdir && (md(z).name(1)~='.')
          nm2 = md(z).name;
          break;
        end
      end
      if isempty(nm2)
        toUse(k) = false;
      else
        d(k).name = [nm nm2];
      end
    end
  end
  d = d(toUse);
  simulationdirs = {d.name}';
end

allcpus = zeros(numel(simulationdirs), 2);
bc = 0 ;
tim = clock;
outi = nargout>2;
if outi
  infos = cell(size(simulationdirs));
  numevs = zeros(size(simulationdirs));
  medevtime = zeros(size(simulationdirs));
end
for k=1:numel(simulationdirs)
  if bc>0; fprintf(repmat('\b', 1, bc)); end; bc = fprintf('Now loading simulation %d of %d: %s, elapsed time: %s', k, numel(simulationdirs), simulationdirs{k}, mat2str(etime(clock, tim)));
  info = reapTimeStatsNew(simulationdirs{k}, false);
  if isempty(info(end).totalClusterReaped) || isempty(info(end).cputimeReaped)
    error('jarl!!');
  end
  allcpus(k, 1) = info(end).cputimeReaped;
  allcpus(k, 2) = info(end).totalClusterReaped;
  if outi
    infos{k} = info;
    load([simulationdirs{k} filesep '..' filesep 'estado.mat']);
    numevs(k) = params.tam_poblacion + sum([info.numStrings]);
    medevtime(k) = info(end).totalClusterReaped / numevs(k);
    clear params;
  end
end

fprintf('\n');