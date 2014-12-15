function individuals = compareSpectralGaps(poph, individuals, removeIdenticals, maxeps, attr, showText)

if ~exist('removeIdenticals', 'var')
  removeIdenticals = false;
end
if ~exist('maxeps', 'var')
  maxeps = 1e-3;%eps(1e3);
end
if ~exist('attr', 'var')
  attr = 'offsetRel';
end
if ~exist('showText', 'var')
  showText = false;
end

spectra = poph.spectra(individuals);

nums = cellfun('prodofsize', spectra);

if removeIdenticals
  toRemove = false(size(nums));
  
  for k=1:numel(spectra)
    if toRemove(k)
      continue;
    end
    %fprintf('%d...\n', k);
    nk = nums(k);
    toInspect = (nk==nums) & (~toRemove);
    toInspect(1:k) = false;
    toInspect = find(toInspect);
    sk = spectra{k};
    for kk=1:numel(toInspect)
      tk   = toInspect(kk);
      stk  = spectra{tk};
      rels = (abs(sk-stk)./min(sk, stk));
      rels(isnan(rels)) = 0;
      toRemove(tk) = all(rels<maxeps);
    end
  end
  
  spectra(toRemove)     = [];
  nums(toRemove)        = [];
  if islogical(individuals)
    individuals         = find(individuals);
  end
  individuals(toRemove) = [];
  
  fprintf('Drawing %d of %d spectra...\n', numel(toRemove)-sum(toRemove), numel(toRemove));
end

if ~isempty(attr)
  [attrs srt] = sort(poph.(attr)(individuals));
  individuals = individuals(srt);
  spectra     = spectra(srt);
  nums        = nums(srt);
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
