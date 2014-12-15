function [idxs rec nav structss] = showGenealogy(conf, poph, idx, onlyMuts)

d3   = size(poph.params.evparams.ss.pos,2)==3;

if ~d3
  error('This is only prepared for 3D!!!!');
end

rec  = makeNewRecorder(poph.params.evparams.ss, 'MemRecorder', 1, 1, 100, 100);

idxs = reshape(trackPopulation(poph, poph.individual(idx), poph.generation(idx), 0, 'all'), [], 1);

[ids ids] = sort(poph.generation(idxs)); %this is not really necessary
idxs      = idxs(ids);

if (nargout<2) && (~onlyMuts)
  return;
end

fname = @(gid)[poph.basedir filesep sprintf('lastFrames%04d_1.mat', gid)];

genomes = cell(size(idxs));

for k=1:numel(idxs)
  id = idxs(k);
  load(fname(poph.generation(id)));
  if numel(res.indexNewG)==numel(res.genomes)
    ind = res.indexNewG==poph.individual(id);
    if any(ind)
      genomes{k} = res.genomes{ind};
    else
      genomes{k} = genomes{k-1};
    end
  else
    ind = poph.individual(id);
    genomes{k} = res.genomes{ind};
  end
end

original    = [true; cellfun(@(a,b)~isequalwithequalnans(a,b), genomes(1:end-1), genomes(2:end))];
generations = poph.generation(idxs);

if onlyMuts
  generations = generations(original);
  genomes     =     genomes(original);
  idxs        =        idxs(original);
end

if (nargout<2)
  return;
end

% poss = zeros(poph.params.genome.fold.numPoints+2, 3, numel(idxs));

for k=1:numel(idxs)
  ss     = showGenome(genomes{k}, poph.params);
  ss.pos = bsxfun(@minus, ss.pos, sum(bsxfun(@times, ss.pos, ss.m/sum(ss.m))));
%   poss(:,:,k) = ss.pos;%reshape(ss.pos, size(ss.pos,1), 3, 1);
  ss.t   = generations(k);
  ss.allStateVars{end+1} = 'str';
  ss.str = [conf.fun.getInfo(poph, idxs(k)) sprintf('\nChange with respect to ancestor: %s', poph.tree.change{idxs(k)})];
  rec    = recordAllState(rec, ss);
end

if nargout<3
  return
end

% dfss = diff(poss,[],3);
% dfss=dfss;
nav = navigateSimulation(PophPlotter('axisWindow', getAxisWindow(rec.notdyn{1}.pos, 0.1), 'axisEqual', true, 'axisSquare', false, 'circleFaces', 10, 'showBalls', 'none', 'showRowContacts', false, 'showCM', false), rec);

if nargout<4
  return
end

structss = arrayfun(@(x)conf.fun.replayGenome(poph, 'poph', x, true, true, 100), idxs, 'uniformoutput', false);