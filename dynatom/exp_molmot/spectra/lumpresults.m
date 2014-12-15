function res = lumpresults


conf = molMot3DConf(false);

bc=0;


% basedir = 'resultados\molmot\';
% baseline = 'Y2010M05D25h16m26s04';
% d = dir(basedir);
% goon = false;
% names   = {};
% g = 1;
% for k=1:numel(d)
%   goon = goon || (d(k).isdir && strcmp(d(k).name, baseline) );
%   if goon && d(k).isdir && (d(k).name(1)=='Y')
%     dd = dir([basedir d(k).name]);
%     for z=1:numel(dd)
%       if dd(z).isdir && (dd(z).name(1)~='.')
%         bd = [basedir d(k).name filesep dd(z).name];
%         if bc>0; fprintf(repmat('\b', 1, bc)); end
%         bc = fprintf('Going for <<%s>>...', bd);
%         poph = conf.fun.loadSim(bd);
%         names{g,1} = bd;
%         maxgens(g,1) = max(poph.generation);
%         bestAbs(g,1) = max(poph.offsetAbs);
%         g=g+1;
%       end
%     end
%   end
% end

names = { ...
'resultados\molmot\Y2010M05D25h16m26s04\1_0.01_1', ...% 3D random defective
'resultados\molmot\Y2010M05D25h16m30s29\1_0.01_1', ...% 3D random defective
'resultados\molmot\Y2010M05D26h13m26s34\1_0.01_1', ...% 3D random defective
'resultados\molmot\Y2010M05D26h16m37s57\1_0.01_1', ...% 3D random defective
'resultados\molmot\Y2010M05D26h16m48s12\1_0.01_1', ...% 3D random defective
'resultados\molmot\Y2010M05D27h13m01s42\1_0.01_1', ...% 3D random defective 
'resultados\molmot\Y2010M05D27h13m01s51\1_0.01_1', ...% 3D random defective
'resultados\molmot\Y2010M05D31h19m32s43\1_0.01_1', ...% 3D random simplest
'resultados\molmot\Y2010M05D31h19m56s22\1_0.01_1', ...% 3D random neighsPassive
'resultados\molmot\Y2010M06D01h16m07s06\1_0.01_1', ...% 3D random neighsPassive
'resultados\molmot\Y2010M06D01h16m07s13\1_0.01_1', ...% 3D random simplest
'resultados\molmot\Y2010M06D02h19m00s49\1_0.01_1', ...% 3D random neighsPassive
'resultados\molmot\Y2010M06D02h19m03s42\1_0.01_1', ...% 3D random simplest
'resultados\molmot\Y2010M06D04h10m05s44\1_0.01_1', ...% 3D random simplest 3spg:0-0.25
'resultados\molmot\Y2010M06D04h10m06s02\1_0.01_1', ...% 3D random neighsPassive 3spg:0-0.25
'resultados\molmot\Y2010M06D05h09m06s05\1_0.01_1', ...% 3D random neighsPassive 3spg:0-0.25
'resultados\molmot\Y2010M06D05h09m11s10\1_0.01_1', ...% 3D random simplest 3spg:0-0.25
'resultados\molmot\Y2010M06D06h08m20s14\1_0.01_1', ...% 3D random simplest 3spg:0-0.25
'resultados\molmot\Y2010M06D06h08m38s56\1_0.01_1', ...% 3D random neighsPassive 3spg:0-0.25
'resultados\molmot\Y2010M06D07h07m39s36\1_0.01_1', ...% 3D random simplest 3spg:2.5-inf
'resultados\molmot\Y2010M06D07h12m19s21\1_0.01_1', ...% 3D random neighsPassive 3spg:2.5-inf
'resultados\molmot\Y2010M06D08h08m40s23\1_0.01_1', ...% 3D random simplest 1spg:2.5-inf
'resultados\molmot\Y2010M06D08h15m22s39\1_0.01_1', ...% 3D random neighsPassive 1spg:2.5-inf
'resultados\molmot\Y2010M06D24h16m19s35\1_0.01_1', ...% 3D random hoh simplest attachmode=1
'resultados\molmot\Y2010M06D24h18m29s44\1_0.01_1', ...% 3D random hoh simplest attachmode=1
'resultados\molmot\Y2010M06D25h08m14s05\1_0.01_1', ...% 3D random hoh simplest attachmode=0
'resultados\molmot\Y2010M06D28h20m39s05\1_0.01_1', ...% 3D random hoh simplest attachmode=0 spchg(softChange)
'resultados\molmot\Y2010M06D29h01m19s44\1_0.01_1', ...% FO 3D random hoh simplest attachmode=0 spchg
'resultados\molmot\Y2010M06D29h15m56s08\1_0.01_1', ...% FO 3D random hoh simplest attachmode=0 spchg
'resultados\molmot\Y2010M06D30h11m34s16\1_0.01_1', ...%    3D random hoh dist2Near simplest attachmode=0 spchg
'resultados\molmot\Y2010M06D30h11m34s25\1_0.01_1', ...% FO 3D random hoh dist2Near simplest attachmode=0 spchg
'resultados\molmot\Y2010M07D02h13m49s44\1_0.01_1', ...%    3D random hoh dist2Near simplest attachmode=0 spchg
'resultados\molmot\Y2010M07D02h13m56s22\1_0.01_1', ...% FO 3D random hoh dist2Near simplest attachmode=0 spchg
...
'resultados\molmot\Y2010M06D09h09m17s13\1_0.01_1', ...% 3D random neighsPassive 8spg:0-0.5 %max 10 space units
'resultados\molmot\Y2010M06D09h12m02s11\1_0.01_1', ...% 3D random neighsPassive 8spg:0-0.5 %max 6 space units
'resultados\molmot\Y2010M06D15h11m37s06\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max 8 space units
'resultados\molmot\Y2010M06D15h12m00s36\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max [42] space units
'resultados\molmot\Y2010M06D15h13m09s36\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max 11 space units
'resultados\molmot\Y2010M06D16h19m12s33\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 16 space units
'resultados\molmot\Y2010M06D16h19m49s24\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 15 space units
'resultados\molmot\Y2010M06D16h21m27s20\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 11 space units
'resultados\molmot\Y2010M06D18h01m55s53\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 13 space units
'resultados\molmot\Y2010M06D18h02m11s13\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 11 space units
'resultados\molmot\Y2010M06D18h03m37s19\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 16 space units
'resultados\molmot\Y2010M06D19h09m32s02\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 10 space units
'resultados\molmot\Y2010M06D19h09m46s08\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max [24] space units
'resultados\molmot\Y2010M06D19h10m48s02\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 fitness plain %max 15 space units
'resultados\molmot\Y2010M06D20h17m11s56\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max 7 space units
'resultados\molmot\Y2010M06D20h17m15s30\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max 12 space units
'resultados\molmot\Y2010M06D20h18m05s50\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max [31] space units
'resultados\molmot\Y2010M06D22h03m05s56\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max 15 space units
'resultados\molmot\Y2010M06D22h03m20s46\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max [53] space units
'resultados\molmot\Y2010M06D22h04m20s28\1_0.01_1', ...% 3D random simplest 8spg:0-0.5 %max 8 space units
}';

maxgens = zeros(size(names));
bestAbs = zeros(size(names));

for k=1:numel(names)
  bd = names{k};
  if bc>0; fprintf(repmat('\b', 1, bc)); end
  bc = fprintf('Processing <<%s>>...', bd);
  poph = conf.fun.loadSim(bd);
  maxgens(k) = max(poph.generation);
  bestAbs(k) = max(poph.offsetAbs);
end

res = struct('name', {names}, 'maxgen', maxgens, 'maxdisp', bestAbs);

fprintf('\n');


fnt = 'Times New Roman'; figure; zp=res.maxdisp; zp(30:49) = []; hist(zp, 10); h = findobj(gca,'Type','patch'); set(h,'FaceColor',0.4*ones(3,1),'EdgeColor','w', 'LineWidth', 1); xlabel('best distance covered', 'FontName', fnt, 'FontSize', 42); set(gcf, 'Color', 'w'); set(gca, 'box', 'off', 'FontName', fnt, 'fontSize', 28, 'layer', 'top');


