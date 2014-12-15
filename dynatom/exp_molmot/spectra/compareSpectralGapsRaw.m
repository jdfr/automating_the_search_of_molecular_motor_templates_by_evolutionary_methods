function [autovalues spectra] = compareSpectralGapsRaw(sss, positions, labels)

spectra  = cell(size(sss));
autovalues = cell(size(sss));

for k=1:numel(sss);
  ss = sss{k};
  
  [autovals autovals autovals] = ANMdecomposition(ss);
  
  autovals = diag(autovals);
  
  bigEnough      = autovals>1e-12;
  bigEnough(1:3) = false;

  autovals = autovals(bigEnough);

  autovalues{k} = autovals;
  spectra{k} = log10(autovals/autovals(1));
end

nums = cellfun('prodofsize', spectra);

ys = cell(size(spectra));
for k=1:numel(ys)
  ys{k} = zeros(nums(k),1)+k;
end
ys = vertcat(ys{:});

spectra = vertcat(spectra{:});

figure;
zys=ys;
if not(exist('labels', 'var'))
  labels = {'B', 'A'};
end
if not(exist('positions', 'var'))
  positions = [2 0; 1 1;];
end
for k=1:size(positions, 1)
  zys(zys==positions(k,1)) = positions(k,2); 
end
d=0.1; 
lns = nan(numel(spectra)*6, 2); 
lns(1:6:end,:) = [spectra-d, zys]; 
lns(2:6:end,:) = [spectra+d, zys]; 
lns(4:6:end,:) = [spectra, zys-d]; 
lns(5:6:end,:) = [spectra, zys+d]; 
figure; 
h=line(lns(:,1),lns(:,2),'Color', 0*ones(1,3), 'LineWidth', 1.5); 
fnt = 'Times New Roman'; 
set(gca, 'GridLineStyle', 'none', 'LineWidth', 1, 'FontName', fnt, ...
    'FontSize', 28, 'box', 'off', 'Layer', 'top', ...
    'XLim', [-0.5 6.5], 'YLim', [-0.5 6.5], 'XTick', [0:6], 'YTick', [0 1], ...
    'YTickLabel', labels);
set(gcf, 'Color', 'w');
xlabel('$\log_{10}(\lambda_{i}/\lambda_{1})$', 'Interpreter','latex', 'FontName', fnt, 'FontSize', 42); 
%figure; h=line(spectra,ys,'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 20, 'MarkerEdgeColor', 0*ones(1,3)); fnt = 'Times New Roman'; set(gca, 'GridLineStyle', 'none', 'LineWidth', 1, 'FontName', fnt, 'FontSize', 28, 'box', 'off', 'Layer', 'top', 'XLim', [-0.5 6.5], 'YLim', [0.5 2.5], 'XTick', [0:6], 'YTick', []); set(gcf, 'Color', 'w');
