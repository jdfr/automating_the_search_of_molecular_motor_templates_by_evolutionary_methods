function datos = analyzeSimulations(conf, anadir, basedir, makePNGs)

if isempty(conf)
  conf = molMot3DConf(false);
end

if ~exist('makePNGs', 'var')
  makePNGs = true;
end

datos = [];

if ischar(basedir)
  nc = fprintf('Processing folder %s...', basedir);

  ds = dir(basedir);

  [names dirs] = arrayfun(@(x)deal(x.name, x.isdir), ds, 'uniformoutput', false);

  dirs = find(cell2mat(dirs));

  fprintf(repmat('\b', 1, nc));

  if all(ismember({'poblacion.txt', 'arbol.txt', 'timestats.txt'}, names)) %&& (~any(strcmp(savename, names)))
    fprintf('Processing simulation %s...\n', basedir);
    fprintf('   Loading params...\n');
    load([basedir filesep '..' filesep 'estado.mat']);
    datitos.basedir     = basedir;
    datitos.useFO       = params.evparams.useFirstOrder;
    datitos.useHOH      = strcmpi(params.evparams.walker.d3.rotationMode, 'handoverhand');
    datitos.fitnessOK   = strcmpi(params.evparams.fitnessCalc, 'mode3');
    datitos.selPACO     = strcmpi(params.selection.mode, 'PACO');
    datitos.fromTerclus = ~isempty(strfind(basedir, 'terclus'));
    fprintf('   Loading simulation...\n');
    poph                = conf.fun.loadSim(basedir);
    poph.basedir        = basedir;
    subana              = sprintf('%s%sHOH=%d_FO=%d_FIT=%d_SEL=%d', anadir, filesep, datitos.useHOH, datitos.useFO, datitos.fitnessOK, datitos.selPACO);
    if ~exist(subana, 'dir')
      mkdir(subana);
    end
    datitos.maxGen      = max(poph.generation);
    if ~isfield(poph, 'ahora')
      seps              = find(basedir==filesep);
      poph.ahora        = basedir((seps(end-1)+1):(seps(end)-1));
    end
    datitos.ahora       = poph.ahora;
    fprintf('   Processing...\n');
    datitos = showHistogram(datitos, subana, basedir, poph, makePNGs);
    datitos = showHistory(  datitos, subana, basedir, poph, makePNGs);
    datos   = [datos; datitos];
    %save([basedir filesep savename], 'poph');
  else
    for k=1:numel(dirs)
      name = names{dirs(k)};
      if ~any(strcmp(name, {'.', '..'}))
        datos = [datos; analyzeSimulations(conf, anadir, [basedir filesep name], makePNGs)]; %#ok<AGROW>
      end
    end
  end

elseif iscell(basedir)
  for k=1:numel(basedir)
    datos = [datos; analyzeSimulations(conf, anadir, basedir{k}, makePNGs)]; %#ok<AGROW>
  end
else
  error('basedir not understood!!!!');
end


function datitos = showHistogram(datitos, subana, basedir, poph, makePNGs)

if ~makePNGs
  return
end

doVisible = 'off';
doSave    = true;
doClose   = true;

nlabs = 10;

flds = {'fitness', 'offsetAbs'};

for k=1:numel(flds)
  fld = flds{k};

  h = figure('Visible', doVisible);
  fld2 = 'generation';
  [n c] = hist3([poph.(fld), poph.(fld2)], [100 100]);
  nlog = log10(n+1);
  ticks = linspace(0, max(nlog(:)), nlabs);
  ticknums = (10.^ticks)-1;

  imagesc(c{2}, c{1}, nlog);
  axis xy;
  xlabel(fld2);
  ylabel(fld);
  colorbar('YTick', ticks, 'YTickLabel', ticknums);

  if doSave;
    nam = ['HISTOGRAM_' fld 'VS' fld2 '.png'];
    nam1 = [basedir filesep nam];
    nam2 = [subana  filesep poph.ahora '_' nam];
    saveas(h, nam1, 'png');
    copyfile(nam1, nam2);
  end;
  if doClose; close(h); end;
end

poph.atpChgMin = min(poph.atpChg1, poph.atpChg2);
poph.atpChgMax = max(poph.atpChg1, poph.atpChg2);
flds = {{'offsetAbs', 'atpChgMin'}, {'offsetAbs', 'atpChgMax'}, {'offsetAbs', 'spectralGap'}};

for k=1:numel(flds)
  fld  = flds{k}{1};
  fld2 = flds{k}{2};

  h = figure('Visible', doVisible);
  [n c] = hist3([poph.(fld), poph.(fld2)], [100 100]);
  nlog = log10(n+1);
  ticks = linspace(0, max(nlog(:)), nlabs);
  ticknums = (10.^ticks)-1;

  imagesc(c{2}, c{1}, nlog);
  axis xy;
  xlabel(fld2);
  ylabel(fld);
  colorbar('YTick', ticks, 'YTickLabel', ticknums);

  if doSave;
    nam = ['HISTOGRAM_' fld 'VS' fld2 '.png'];
    nam1 = [basedir filesep nam];
    nam2 = [subana  filesep poph.ahora '_' nam];
    saveas(h, nam1, 'png');
    copyfile(nam1, nam2);
  end;
  if doClose; close(h); end;
end


function datitos = showHistory(datitos, subana, basedir, poph, makePNGs)

doVisible = 'off';
doSave    = true;
doClose   = true;

flds = {'fitness', 'offsetAbs'};
otherflds = {'atpChg1', 'atpChg2', 'offsetRel', 'spectralGap'};

for k=1:numel(flds)
  fld = flds{k};
  
  mxfld  = ['mx_'  fld];
  mxfldi = ['mxi_' fld];
  [datitos.(mxfld) idx] = max(poph.(fld));
  
  for z=1:numel(otherflds)
    ofld = otherflds{z};
    datitos.([mxfld '_' ofld]) = poph.(ofld)(idx);
  end
  
  datitos.(mxfldi) = idx;
  idxs = reshape(trackPopulation(poph, poph.individual(idx), poph.generation(idx), 0, 'all'), [], 1);

  [ids ids] = sort(poph.generation(idxs)); %this is not really necessary
  idxs      = idxs(ids);
  
  mxidxsfld = ['mxhist_' fld];
  
  datitos.(mxidxsfld) = idxs;
  
  if makePNGs
    for z=1:numel(flds)
      fld2use = flds{z};

      h = figure('Visible', doVisible);
      fld2 = 'generation';
      hh = stairs(poph.(fld2)(idxs), poph.(fld2use)(idxs));
      set(hh, 'Color', 'b', 'LineWidth', 2);
      grid on;

      xlabel(fld2);
      ylabel(fld2use);

      if doSave;
        nam = ['HISTORY_' fld2use '_FOR_mx_' fld '.png'];
        nam1 = [basedir filesep nam];
        nam2 = [subana  filesep poph.ahora '_' nam];
        saveas(h, nam1, 'png');
        copyfile(nam1, nam2);
      end;
      if doClose; close(h); end;
    end
  end

end
