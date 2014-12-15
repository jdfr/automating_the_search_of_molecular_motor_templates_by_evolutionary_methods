function out = getLineageOfChanges(poph, individual, mode, output)

indexesALL = trackPopulation(poph, poph.individual(individual), poph.generation(individual), min(poph.generation), 'all');
gs         = poph.generation(indexesALL);
[gs gs]    = sort(gs, 'descend');
indexesALL = indexesALL(gs);

if ~iscell(mode)
  mode = {mode};
end

keys = cell(size(mode));

for k=1:numel(keys)
  switch mode{k}
    case 'fitness'
      roundoff            = eps(1e6);
      keyIndividuals      = [indexesALL( abs(diff(poph.fitness(indexesALL))) > roundoff ); indexesALL(end)];
    case 'genome'
      [index index]       = ismember(poph.idx(indexesALL), poph.tree.idxD);
      keyIndividuals      = ~cellfun(@(x)strcmp(x,'='), poph.tree.change(index));
      keyIndividuals(end) = true;
      keyIndividuals      = indexesALL(keyIndividuals);
    otherwise
      error('mode==%s not understood!!!!', any2str(mode));
  end
  keys{k} = keyIndividuals;
end

keyIndividuals = keys{1};
for k=2:numel(keys)
  keyIndividuals = intersect(keyIndividuals, keys{k});
end

switch output
  case 'idxpoph'
    out = keyIndividuals;
  case 'evaluated'
    out = cell(size(keyIndividuals));
    for k=1:numel(keyIndividuals)
      fprintf('Let''s go with individual #%d\n', k);
      n = keyIndividuals(k);
      out{k} = evaluarMolMot(poph.params.genome2str(poph.genome{n}), poph.params.evparams, poph.rndSeedDev{n}, poph.rndSeedEv{n}, true, true);
    end
  otherwise
    error('output==%s not understood!!!!', any2str(output));
end