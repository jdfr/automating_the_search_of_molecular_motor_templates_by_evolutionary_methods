function makeGraphEvoRun(pnames)
conf = molMot3DConf(false);

%zn=tic; cput = cputime; datos = analyzeSimulations([], 'resultados\_parapaper\ana', 'resultados\_parapaper'); toc(zn), zcput=cputime-cput
% tit = 'Monomer configuration'; zsel =  ~[datos.useHOH] & ~[datos.useFO] & [datos.fitnessOK]; zmxoff = [datos.mx_offsetAbs]; zmxoff = zmxoff(zsel); fnt = 'Times New Roman'; figure; zp=zmxoff; hist(zp, 20:10:120); h = findobj(gca,'Type','patch'); set(h,'FaceColor',0.4*ones(3,1),'EdgeColor','w', 'LineWidth', 1); xlabel('best distance covered', 'FontName', fnt, 'FontSize', 42); set(gcf, 'Color', 'w'); set(gca, 'box', 'off', 'FontName', fnt, 'fontSize', 28, 'layer', 'top'); uicontrol('Style', 'text', 'fontName', fnt, 'fontSize', 42, 'BackgroundColor', 'w', 'Units', 'normalized', 'Position', [0.6 0.7 0.3 0.25], 'String', tit);
% tit = 'Dimer configuration';   zsel =   [datos.useHOH] & ~[datos.useFO] & [datos.fitnessOK]; zmxoff = [datos.mx_offsetAbs]; zmxoff = zmxoff(zsel); fnt = 'Times New Roman'; figure; zp=zmxoff; hist(zp, 20:10:120); h = findobj(gca,'Type','patch'); set(h,'FaceColor',0.4*ones(3,1),'EdgeColor','w', 'LineWidth', 1); xlabel('best distance covered', 'FontName', fnt, 'FontSize', 42); set(gcf, 'Color', 'w'); set(gca, 'box', 'off', 'FontName', fnt, 'fontSize', 28, 'layer', 'top'); uicontrol('Style', 'text', 'fontName', fnt, 'fontSize', 42, 'BackgroundColor', 'w', 'Units', 'normalized', 'Position', [0.6 0.7 0.3 0.25], 'String', tit);

%zsel = ~[datos.useFO] & [datos.fitnessOK]; basedirs = {datos.basedir}; basedirs = basedirs(zsel); tic; makeGraphEvoRun(basedirs); toc
showLumpedSpectra(pnames, conf);
return




%namedir = 'Y2010M06D30h11m34s16';
namedir = 'Y2010M06D01h16m07s06';
%namedir = 'Y2010M07D12h11m42s54';

poph = conf.fun.loadSim(['resultados\molmot\' namedir '\1_0.01_1']);

ming = min(poph.generation);
maxg = max(poph.generation);

maxgs = zeros(maxg-ming+1,1);
mings = maxgs;
meangs = maxgs;
stdgs = maxgs;

g = ming;
fld = 'offsetAbs';
for k=1:numel(maxgs);
  sel = poph.(fld)(poph.generation==g);
  sel = sel(~isnan(sel));
  maxgs(k) = max(sel);
  mings(k) = min(sel);
  meangs(k) = mean(sel);
  stdgs(k) = std(sel, 1);
  g=g+1;
end;

figure;
hold on;
[X1 Y1] = stairs(ming:maxg, meangs+stdgs);
[X2 Y2] = stairs(ming:maxg, meangs-stdgs);
[F V] = poly2fv([X1(:); X2(end:-1:1)], [Y1(:); Y2(end:-1:1)]);
h = patch(struct('faces', F, 'vertices', V), 'EdgeColor', 'none', 'FaceColor', 0.9*ones(1,3));
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

% h = line(X1,Y1,               'Color', 0.5*ones(1,3), 'LineWidth', 2);
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
% h = line(X2,Y2,               'Color', 0.5*ones(1,3), 'LineWidth', 2);
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

stairs(ming:maxg, maxgs,  'Color', 0.0*ones(1,3),  'LineWidth', 2);
stairs(ming:maxg, meangs, 'Color', 0.6*ones(1,3), 'LineWidth', 2);
stairs(ming:maxg, mings,  'Color', 0.4*ones(1,3), 'LineWidth', 2);
grid on;
set(gcf, 'Color', ones(1,3));
fnt = 'Times New Roman';
set(gca, 'GridLineStyle', 'none', 'LineWidth', 3, 'FontName', fnt, 'FontSize', 28, 'Layer', 'top');
xlabel('generations', 'FontName', 'Times New Roman', 'FontSize', 42);
ylabel('absolute distance covered', 'FontName', fnt, 'FontSize', 42);

h = legend({'max', 'mean', 'min'}, 'Location', 'NorthWest');

pp = get(h, 'Position');

set(h, 'FontName', fnt, 'FontSize', 28, 'layer', 'top', 'Color', 'none', 'box', 'off');

%set(h, 'Position', [pp(1:2), pp(3:4)*1.5]);

%legend('boxoff');


function showLumpedSpectra(pnames, conf)

%pnames = {'Y2010M05D25h16m26s04', 'Y2010M05D26h13m26s34', 'Y2010M05D26h16m37s57', 'Y2010M05D26h16m48s12', 'Y2010M05D27h13m01s42', 'Y2010M05D27h13m01s51', 'Y2010M05D28h10m16s23', 'Y2010M05D28h10m16s38', 'Y2010M05D29h10m32s12', 'Y2010M05D29h10m35s53', 'Y2010M05D30h10m03s07', 'Y2010M05D30h10m05s41'}';

pophs = cell(size(pnames));
for k=1:numel(pophs)
  %load(['resultados\molmot\', pnames{k}, '\1_0.01_1\pophSG.mat']);
  poph = conf.fun.loadSim(pnames{k});
  pophs{k} = poph;
end

% namegap = 'trdSpecGap'; for k=1:numel(pophs); poph = pophs{k}; figure; ftn = poph.offsetAbs; ftn(isnan(ftn)) = -10; hist3([ftn, poph.(namegap)], [100 100]); grid on; xlabel('offsetAbs'); ylabel(namegap); end;
% 
% namegap = 'trdSpecGap';  [lumpeda lumpedg] = cellfun(@(x)deal(x.offsetAbs, x.(namegap)), pophs, 'uniformoutput', false); lumpeda = vertcat(lumpeda{:}); lumpedg = vertcat(lumpedg{:}); lumpeda(isnan(lumpeda)) = -10; figure; hist3([lumpeda(lumpeda>80), lumpedg(lumpeda>80)], [100 100]); grid on; xlabel('offsetAbs'); ylabel(namegap);
% 
% namegap = 'trdSpecGap';  [lumpeda lumpedg] = cellfun(@(x)deal(x.offsetAbs, x.(namegap)), pophs, 'uniformoutput', false); lumpeda = vertcat(lumpeda{:}); lumpedg = vertcat(lumpedg{:}); lumpeda(isnan(lumpeda)) = -10; figure; hist3([lumpeda, lumpedg], [100 100]); grid on; xlabel('offsetAbs'); ylabel(namegap);
% 
% namegap = 'trdSpecGap';  [lumpedao lumpedgo] = cellfun(@(x)deal(x.offsetAbs(x.tree.idxA(1:numel(x.idx))==0), x.(namegap)((x.tree.idxA(1:numel(x.idx))==0))), pophs, 'uniformoutput', false); lumpedao = vertcat(lumpedao{:}); lumpedgo = vertcat(lumpedgo{:}); lumpedao(isnan(lumpedao)) = -10; figure; hist3([lumpedao, lumpedgo], [100 100]); grid on; xlabel('offsetAbs'); ylabel(namegap);
% 
namegap = 'spectralGap'; 'trdSpecGap';  [lumpedao lumpedgo] = cellfun(@(x)deal(x.offsetAbs(x.tree.idxA(1:numel(x.idx))==0), x.(namegap)((x.tree.idxA(1:numel(x.idx))==0))), pophs, 'uniformoutput', false); lumpedao = vertcat(lumpedao{:}); lumpedgo = vertcat(lumpedgo{:}); lumpedao(isnan(lumpedao)) = -10; figure;                            hist(lumpedgo, linspace(0, 3.5, 100)); fnt = 'Times New Roman'; set(gca, 'GridLineStyle', 'none', 'LineWidth', 1, 'FontName', fnt, 'FontSize', 28, 'box', 'off', 'Layer', 'top', 'XLim', [0 3.5]); ylabel('number of structures', 'FontName', fnt, 'FontSize', 42); xlabel('third spectral gap: $\log_{10}(\lambda_{4}/\lambda_{1})$', 'Interpreter','latex', 'FontName', fnt, 'FontSize', 42); h = findobj(gca,'Type','patch'); set(h,'FaceColor',0.5*ones(1,3),'EdgeColor','w', 'LineWidth', 1); set(gcf, 'Color', 'w');
numel(lumpedgo)
prctile(lumpedgo,[1 99])
uicontrol('Style', 'text', 'fontName', fnt, 'fontSize', 42, 'BackgroundColor', 'w', 'Units', 'normalized', 'Position', [0.5 0.5 0.3 0.2], 'String', 'random structures');
set(gca, 'Position', [0.1300    0.1500    0.7750    0.8150])
%namegap = 'trdSpecGap';  [lumpedao lumpedgo] = cellfun(@(x)deal(x.offsetAbs(x.tree.idxA(1:numel(x.idx))~=0), x.(namegap)((x.tree.idxA(1:numel(x.idx))~=0))), pophs, 'uniformoutput', false); lumpedao = vertcat(lumpedao{:}); lumpedgo = vertcat(lumpedgo{:}); lumpedao(isnan(lumpedao)) = -10; figure; lumpedgo(lumpedgo>3.5)=[]; hist(lumpedgo, linspace(0, 3.5, 100)); fnt = 'Times New Roman'; set(gca, 'GridLineStyle', 'none', 'LineWidth', 1, 'FontName', fnt, 'FontSize', 28, 'box', 'off', 'Layer', 'top', 'XLim', [0 3.5]); ylabel('number of structures', 'FontName', fnt, 'FontSize', 42); xlabel('third spectral gap: $\log_{10}(\lambda_{4}/\lambda_{1})$', 'Interpreter','latex', 'FontName', fnt, 'FontSize', 42); h = findobj(gca,'Type','patch'); set(h,'FaceColor',0.5*ones(1,3),'EdgeColor','w', 'LineWidth', 1); set(gcf, 'Color', 'w');
return
namegap = 'spectralGap'; 'trdSpecGap';  [lumpedao lumpedgo] = cellfun(@(x)deal(x.offsetAbs(x.tree.idxA(1:numel(x.idx))~=0), x.(namegap)((x.tree.idxA(1:numel(x.idx))~=0))), pophs, 'uniformoutput', false); lumpedao = vertcat(lumpedao{:}); lumpedgo = vertcat(lumpedgo{:}); lumpedao(isnan(lumpedao)) = -10; figure; lumpedgo(lumpedao<50)=[];  hist(lumpedgo, linspace(0, 3.5, 100)); fnt = 'Times New Roman'; set(gca, 'GridLineStyle', 'none', 'LineWidth', 1, 'FontName', fnt, 'FontSize', 28, 'box', 'off', 'Layer', 'top', 'XLim', [0 3.5]); ylabel('number of structures', 'FontName', fnt, 'FontSize', 42); xlabel('third spectral gap: $\log_{10}(\lambda_{4}/\lambda_{1})$', 'Interpreter','latex', 'FontName', fnt, 'FontSize', 42); h = findobj(gca,'Type','patch'); set(h,'FaceColor',0.5*ones(1,3),'EdgeColor','w', 'LineWidth', 1); set(gcf, 'Color', 'w');
numel(lumpedgo)
uicontrol('Style', 'text', 'fontName', fnt, 'fontSize', 42, 'BackgroundColor', 'w', 'Units', 'normalized', 'Position', [0.5 0.5 0.3 0.2], 'String', 'evolved structures'); %with repetition
set(gca, 'Position', [0.1300    0.1500    0.7750    0.8150])

namegap = 'spectralGap'; 'trdSpecGap';  [lumpedao lumpedgo lumpednotchanged] = cellfun(@(x)deal(x.offsetAbs(x.tree.idxA(1:numel(x.idx))~=0), x.(namegap)((x.tree.idxA(1:numel(x.idx))~=0)), ~strcmp('=', x.tree.change(x.tree.idxA(1:numel(x.idx))~=0))), pophs, 'uniformoutput', false); lumpedao = vertcat(lumpedao{:}); lumpedgo = vertcat(lumpedgo{:}); lumpednotchanged = vertcat(lumpednotchanged{:}); lumpedao(isnan(lumpedao)) = -10; figure; lumpedgo((lumpedao<50) | lumpednotchanged)=[];  hist(lumpedgo, linspace(0, 3.5, 100)); fnt = 'Times New Roman'; set(gca, 'GridLineStyle', 'none', 'LineWidth', 1, 'FontName', fnt, 'FontSize', 28, 'box', 'off', 'Layer', 'top', 'XLim', [0 3.5]); ylabel('number of structures', 'FontName', fnt, 'FontSize', 42); xlabel('third spectral gap: $\log_{10}(\lambda_{4}/\lambda_{1})$', 'Interpreter','latex', 'FontName', fnt, 'FontSize', 42); h = findobj(gca,'Type','patch'); set(h,'FaceColor',0.5*ones(1,3),'EdgeColor','w', 'LineWidth', 1); set(gcf, 'Color', 'w');
numel(lumpedgo)
uicontrol('Style', 'text', 'fontName', fnt, 'fontSize', 42, 'BackgroundColor', 'w', 'Units', 'normalized', 'Position', [0.5 0.5 0.3 0.2], 'String', 'evolved structures'); %without repetition
set(gca, 'Position', [0.1300    0.1500    0.7750    0.8150])