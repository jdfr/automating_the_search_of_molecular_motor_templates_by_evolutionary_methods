function ret = analyzePoph(poph, fld, doSave, doVisible, doClose, doPrint)

if ~exist('fld', 'var')
  fld = 'offsetAbs';
end
if ~exist('doSave', 'var')
  doSave = false;
end
if ~exist('doVisible', 'var')
  doVisible = true;
end
if ~exist('doClose', 'var')
  doClose = false;
end
if ~exist('doPrint', 'var')
  doPrint = false;
end

if iscell(fld)
  ret = cellfun(@(fld)analyzePoph(poph, fld, doSave, doVisible, doClose, doPrint), fld, 'UniformOutput', false);
  ret = [ret{:}];
  return
end

if doPrint
  fprintf('Processing %s...\n', poph.basedir);
end

if doVisible
  doVisible = 'on';
else
  doVisible = 'off';
end

[mxi mx] = max(poph.(fld));

f = fopen([poph.basedir filesep fld '_best.txt'], 'w');
s=sprintf('Best %s: %s -> %s\n', fld, mat2str(mx), mat2str(mxi));
fprintf(f, s);
fclose(f);


if strcmp(fld, 'offsetAbs')
  if isfield(poph.params.evparams.walker, 'd3') && isfield(poph.params.evparams, 'useFirstOrder')
    dat = sprintf('rotation=%s, useFirstOrder=%s\n', poph.params.evparams.walker.d3.rotationMode, mat2str(poph.params.evparams.useFirstOrder));
  else
    dat = '';
  end
  ret = [poph.basedir sprintf('\n') ...
         dat ...
         s];
else
  ret = '';
end

nlabs = 10;

h = figure('Visible', doVisible);
fld2 = 'generation';
[n c] = hist3([poph.(fld), poph.(fld2)], [100 100]);
nlog = log10(n+1);
ticks = linspace(0, max(nlog(:)), nlabs);
ticknums = (10.^ticks)-1;

imagesc(c{1}, c{2}, nlog');
xlabel(fld);
ylabel(fld2);
colorbar('YTick', ticks, 'YTickLabel', ticknums);

if doSave; saveas(h, [poph.basedir filesep fld 'VS' fld2 '.png'], 'png'); end;
if doClose; close(h); end;

fld2 = 'sum_adv';
if isfield(poph, fld2) && ~all(isnan(poph.(fld2)))
  h = figure('Visible', doVisible);
  [n c] = hist3([poph.(fld2), poph.(fld)], [100 100]);
  nlog = log10(n+1);
  ticks = linspace(0, max(nlog(:)), nlabs);
  ticknums = (10.^ticks)-1;
  imagesc(c{1}, c{2}, nlog');
  fld3 = fld2; fld3(fld3=='_') = ' ';
  xlabel(fld3);
  ylabel(fld);
  colorbar('YTick', ticks, 'YTickLabel', ticknums);

  if doSave; saveas(h, [poph.basedir filesep fld 'VS' fld2 '.png'], 'png'); end;
  if doClose; close(h); end;
end

fld2 = 'num_adv';
if isfield(poph, fld2) && ~all(isnan(poph.(fld2)))
  h = figure('Visible', doVisible);
  [n c] = hist3([poph.(fld2), poph.(fld)], [100 100]);
  nlog = log10(n+1);
  ticks = linspace(0, max(nlog(:)), nlabs);
  ticknums = (10.^ticks)-1;
  imagesc(c{1}, c{2}, nlog');
  fld3 = fld2; fld3(fld3=='_') = ' ';
  xlabel(fld3);
  ylabel(fld);
  colorbar('YTick', ticks, 'YTickLabel', ticknums);

  if doSave; saveas(h, [poph.basedir filesep fld 'VS' fld2 '.png'], 'png'); end;
  if doClose; close(h); end;
end
