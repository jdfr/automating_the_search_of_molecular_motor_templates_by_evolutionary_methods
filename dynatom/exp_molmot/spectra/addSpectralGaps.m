function addSpectralGaps(conf, savename, basedir)

if ischar(basedir)
  nc = fprintf('Processing folder %s...', basedir);

  ds = dir(basedir);

  [names dirs] = arrayfun(@(x)deal(x.name, x.isdir), ds, 'uniformoutput', false);

  dirs = find(cell2mat(dirs));

  fprintf(repmat('\b', 1, nc));

  if all(ismember({'poblacion.txt', 'arbol.txt', 'timestats.txt'}, names)) && (~any(strcmp(savename, names)))
    fprintf('Processing simulation %s...\n', basedir);
    fprintf('   Loading...\n');
    poph = conf.fun.loadSim(basedir);
    fprintf('   Processing...\n');
    poph = addSpectralGapToPoph(poph);
    fprintf('   Saving...\n');
    save([basedir filesep savename], 'poph');
  end

  for k=1:numel(dirs)
    name = names{dirs(k)};
    if ~any(strcmp(name, {'.', '..'}))
      addSpectralGaps(conf, savename, [basedir filesep name]);
    end
  end

elseif iscell(basedir)
  for k=1:numel(basedir)
    addSpectralGaps(conf, savename, basedir{k})
  end
else
  error('basedir not understood!!!!');
end
