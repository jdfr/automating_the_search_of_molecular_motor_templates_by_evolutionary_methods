function [srt autovalues spectra] = compareSpectralGapsInResults(results, mode, attr, showText)

if ~exist('attr', 'var')
  attr = 'offsetRel';
end
if ~exist('showText', 'var')
  showText = false;
end

spectra  = cell(size(results.structs));
autovalues = cell(size(results.structs));

for k=1:numel(results.structs);
  switch mode
    case 'dev'
      ss = results.structs{k}.ssDev;
    case 'legs'
      ss = results.structs{k}.ssLegs;
    otherwise
      ss = [];
  end
  [autovals autovals autovals] = ANMdecomposition(ss);
  
  autovals = diag(autovals);
  
  bigEnough      = autovals>1e-12;
  bigEnough(1:3) = false;

  autovals = autovals(bigEnough);

  autovalues{k} = autovals;
  spectra{k} = log10(autovals/autovals(1));
end

nums = cellfun('prodofsize', spectra);

individuals = results.inds;

if ~isempty(attr)
  [attrs srt] = sort(cellfun(@(x)x.stats.(attr), results.structs));
  individuals = individuals(srt);
  autovalues  = autovalues(srt);
  spectra     = spectra(srt);
  nums        = nums(srt);
else
  srt         = 1:numel(individuals);
end

ys = cell(size(spectra));
for k=1:numel(ys)
  ys{k} = zeros(nums(k),1)+k;
end
ys = vertcat(ys{:});

spectra = vertcat(spectra{:});

figure;
if isempty(attr)
  plot(spectra,ys,'+');
else
  plot(spectra,ys,'+',   -attrs, (1:numel(individuals))', '*');
  %plot(spectra,ys,'+');
  if showText
    hold on;
    for k=1:numel(individuals)
      text(0, k, [num2str(attrs(k), '%f') '\rightarrow'], 'HorizontalAlignment', 'right');
    end
  end
end
grid on;
