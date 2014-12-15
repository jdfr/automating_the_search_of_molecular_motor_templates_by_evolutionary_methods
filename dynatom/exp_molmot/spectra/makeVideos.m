function makeVideos(toMovie, nfile)

maxiwindow = [-2.431845952794960e+02,3.237891523556778e+02,-2.587753300365620e+02,3.081984175986120e+02,-2.789713804342808e+02,2.880023672008930e+02];
%maxiwindow = [-2.5e+02,3.3e+02,-2.6e+02,3.1e+02,-2.8e+02,2.9e+02];

resA = [1200 1200];
resB = [600  600];

p = ElegantMolPlotter('axisWindow', maxiwindow, 'axisEqual', true, 'axisSquare', false, 'rowSegment', -200:200, 'moveCM', true, 'AVIWidthHeight', resA, 'AVImethod', 'IMAGES');
%p = ElegantMolPlotter('axisWindow', [-2.5e+02,3.3e+02,-2.6e+02,3.1e+02,-2.8e+02,2.9e+02], 'axisEqual', true, 'axisSquare', false, 'rowSegment', -200:200, 'moveCM', true, 'AVIWidthHeight', [1200 1200], 'AVImethod', 'IMAGES', 'AVIfilename', '_render\\aaa%04d.png', 'handleImage', @(filename, img, numFrame) imwrite(imresize(uint8(sum(double(img), 3)/3), 'OutputSize', [600 600], 'Method', 'lanczos2', 'Antialiasing', true), sprintf(filename, numFrame), 'png'));

cmd1 = 'poph.params.evparams.walker.row.spchg = struct(''doit'', true, ''constT'', 2);';
cmd2 = 'poph.params.evparams.walker.row.spchg = struct(''doit'', true, ''constT'', 1);';

individuals = {
  001,  'Y2010M07D02h13m56s22',                    18913, false(1,3), [0 0 0],      {'manual', [-1609.98744355839 -568.416948616827 -44.0252559519688], 'manual', [2.93220857994939 -7.94145730972313 -0.193084525156027], 'manual', [0.328009282908032 -0.944632443447937 0.0089138718713917], 'manual', 1.47918473917376, 'orthographic'}; ...
  002,  'Y2010M07D02h13m49s44',                    18966, false(1,3), [0 0 0],      {'manual', [524.332196235205 -930.721021860282 2399.86809105261], 'manual', [2.27882896297268 -22.3440627174956 -0.0343458249214397], 'manual', [-0.0737350073748244 -0.93790573648451 -0.338963092617593], 'manual', 1.27480534639727, 'orthographic'}; ...
  003,  'Y2010M05D27h13m01s51',                    27969, true(1,3),  [0 0 0],      {'manual', [-289.473139938819 293.637443039248 -2416.85267668576], 'manual', [19.7217389193902 13.150729200089 1.92769036373692], 'manual', [0.0144896584971116 0.993449435982313 0.113350200454444], 'manual', 1.4, 'orthographic'}; ...
  004,  'Y2010M05D31h19m32s43',                    29861, false(1,3), [-2.5 0 0],   {'manual', [-82.7415765620083 441.620651086652 2090.80066250783], 'manual', [18.6599697740972 12.4853912093597 -2.02206292718975], 'manual', [0.00971029977049371 0.979663889260798 -0.200410514087395], 'manual', 1.5, 'orthographic'}; ...
  005,  'Y2010M06D24h18m29s44',                    39726, false(1,3), [0 0 0],      {'manual', [-1858.93546064298 -756.816567054286 -10.4877753571449], 'manual', [9.6375122263948 -15.0642812276464 1.03445135674494], 'manual', [0.368942059493325 -0.929449611887302 0.00227501634428739], 'manual', 1.4584470217686, 'orthographic'}; ...
  006,  'Y2010M06D22h03m20s46',                    39909, true(1,3),  [0 0 0],      {'manual', [-343.446679525868 260.913394484723 -1344.18198750208], 'manual', [9.80940723650247 7.88742910280804 2.64702571525161], 'manual', [0.045360794991399 0.983886706429992 0.172943190638116], 'manual', 1.3563557302448, 'orthographic'}; ...
  007,  'Y2010M07D16h16m24s02',                    27780, true(1,3),  [0 0 0],      {'manual', [-642.964396430556 -924.361845874865 -1939.91617692022], 'manual', [4.72783363637917 -11.4430497383306 -0.0679797514712766], 'manual', [0.129093977388964 -0.913151386522491 0.386638449063139], 'manual', 1.45, 'orthographic'}; ...
  008, {'Y2010M07D02h13m56s22', '_constT2', cmd1}, 18913, false(1,3), [0 0 0],      {'manual', [-1609.98744355839 -568.416948616827 -44.0252559519688], 'manual', [2.93220857994939 -7.94145730972313 -0.193084525156027], 'manual', [0.328009282908032 -0.944632443447937 0.0089138718713917], 'manual', 1.47918473917376, 'orthographic'}; ...
  009, {'Y2010M07D02h13m49s44', '_constT2', cmd1}, 18966, false(1,3), [0 0 0],      {'manual', [524.332196235205 -930.721021860282 2399.86809105261], 'manual', [2.27882896297268 -22.3440627174956 -0.0343458249214397], 'manual', [-0.0737350073748244 -0.93790573648451 -0.338963092617593], 'manual', 1.27480534639727, 'orthographic'}; ...
  010, {'Y2010M05D27h13m01s51', '_constT2', cmd1}, 27969, true(1,3),  [0 0 0],      {'manual', [-289.473139938819 293.637443039248 -2416.85267668576], 'manual', [19.7217389193902 13.150729200089 1.92769036373692], 'manual', [0.0144896584971116 0.993449435982313 0.113350200454444], 'manual', 1.4, 'orthographic'}; ...
  011, {'Y2010M05D31h19m32s43', '_constT2', cmd1}, 29861, false(1,3), [-2.5 0 0],   {'manual', [-82.7415765620083 441.620651086652 2090.80066250783], 'manual', [18.6599697740972 12.4853912093597 -2.02206292718975], 'manual', [0.00971029977049371 0.979663889260798 -0.200410514087395], 'manual', 1.5, 'orthographic'}; ...
  012, {'Y2010M06D24h18m29s44', '_constT2', cmd1}, 39726, false(1,3), [0 0 0],      {'manual', [-1858.93546064298 -756.816567054286 -10.4877753571449], 'manual', [9.6375122263948 -15.0642812276464 1.03445135674494], 'manual', [0.368942059493325 -0.929449611887302 0.00227501634428739], 'manual', 1.4584470217686, 'orthographic'}; ...
  013, {'Y2010M06D22h03m20s46', '_constT2', cmd1}, 39909, true(1,3),  [0 0 0],      {'manual', [-343.446679525868 260.913394484723 -1344.18198750208], 'manual', [9.80940723650247 7.88742910280804 2.64702571525161], 'manual', [0.045360794991399 0.983886706429992 0.172943190638116], 'manual', 1.3563557302448, 'orthographic'}; ...
  014, {'Y2010M07D16h16m24s02', '_constT2', cmd1}, 27780, true(1,3),  [0 0 0],      {'manual', [-642.964396430556 -924.361845874865 -1939.91617692022], 'manual', [4.72783363637917 -11.4430497383306 -0.0679797514712766], 'manual', [0.129093977388964 -0.913151386522491 0.386638449063139], 'manual', 1.45, 'orthographic'}; ...
  015, {'Y2010M07D02h13m56s22', '_constT1', cmd2}, 18913, false(1,3), [0 0 0],      {'manual', [-1609.98744355839 -568.416948616827 -44.0252559519688], 'manual', [2.93220857994939 -7.94145730972313 -0.193084525156027], 'manual', [0.328009282908032 -0.944632443447937 0.0089138718713917], 'manual', 1.47918473917376, 'orthographic'}; ...
  016, {'Y2010M07D02h13m49s44', '_constT1', cmd2}, 18966, false(1,3), [0 0 0],      {'manual', [524.332196235205 -930.721021860282 2399.86809105261], 'manual', [2.27882896297268 -22.3440627174956 -0.0343458249214397], 'manual', [-0.0737350073748244 -0.93790573648451 -0.338963092617593], 'manual', 1.27480534639727, 'orthographic'}; ...
  017, {'Y2010M05D27h13m01s51', '_constT1', cmd2}, 27969, true(1,3),  [0 0 0],      {'manual', [-289.473139938819 293.637443039248 -2416.85267668576], 'manual', [19.7217389193902 13.150729200089 1.92769036373692], 'manual', [0.0144896584971116 0.993449435982313 0.113350200454444], 'manual', 1.4, 'orthographic'}; ...
  018, {'Y2010M05D31h19m32s43', '_constT1', cmd2}, 29861, false(1,3), [-2.5 0 0],   {'manual', [-82.7415765620083 441.620651086652 2090.80066250783], 'manual', [18.6599697740972 12.4853912093597 -2.02206292718975], 'manual', [0.00971029977049371 0.979663889260798 -0.200410514087395], 'manual', 1.5, 'orthographic'}; ...
  019, {'Y2010M06D24h18m29s44', '_constT1', cmd2}, 39726, false(1,3), [0 0 0],      {'manual', [-1858.93546064298 -756.816567054286 -10.4877753571449], 'manual', [9.6375122263948 -15.0642812276464 1.03445135674494], 'manual', [0.368942059493325 -0.929449611887302 0.00227501634428739], 'manual', 1.4584470217686, 'orthographic'}; ...
  020, {'Y2010M06D22h03m20s46', '_constT1', cmd2}, 39909, true(1,3),  [0 0 0],      {'manual', [-343.446679525868 260.913394484723 -1344.18198750208], 'manual', [9.80940723650247 7.88742910280804 2.64702571525161], 'manual', [0.045360794991399 0.983886706429992 0.172943190638116], 'manual', 1.3563557302448, 'orthographic'}; ...
  021, {'Y2010M07D16h16m24s02', '_constT1', cmd2}, 27780, true(1,3),  [0 0 0],      {'manual', [-642.964396430556 -924.361845874865 -1939.91617692022], 'manual', [4.72783363637917 -11.4430497383306 -0.0679797514712766], 'manual', [0.129093977388964 -0.913151386522491 0.386638449063139], 'manual', 1.45, 'orthographic'}; ...
  };

sims = individuals(:,2);
inds = individuals(:,3);
vars = individuals(:,4);
sfht = individuals(:,5);
cams = individuals(:,6);

p.cylinderDetail = 40;

conf = molMot3DConf(false);

if ~exist('nfile', 'var')
  nfile  = 'gait';
end
nomdir = ['_render\\%s.' nfile];
p.handleImage = handleImageFUN(resB);

for k=1:numel(toMovie)
  tic;
  z               = toMovie(k);
  p.toVaryPoint   = vars{z};
  p.shiftCamera   = sfht{z};
  ind             = inds{z};
  indstr          = mat2str(ind);
  setCameraProps(p, 'before', cams{z});
  ndir            = sprintf(nomdir, indstr);
  if exist(ndir, 'dir')
    rmdir(ndir, 's');
  end
  mkdir(ndir);
  ndir = strrep(ndir, '\', '\\');
  p.AVIfilename   = [ndir filesep filesep indstr '.' nfile '.%04d.png'];
  if ischar(sims{z})
    load(['Z_' sims{z}]);
    q = find(results.inds==ind);
    structs = results.structs{q};
    clear results;
  else
    path = ['resultados\molmot\' sims{z}{1} '\1_0.01_1'];
    poph = conf.fun.loadSim(path);
    eval(sims{z}{3});
    structs = conf.fun.replayGenome(poph, 'poph', ind, [false true], true, 100);
    clear poph;
    results = struct('path', path, 'mode', ['modified: <<' sims{z}{3} '>>'], 'inds', ind, 'structs', {{structs}});
    save(['Z_' sims{z}{1:2} ], 'results');
    clear results;
  end
  TS = allTimeSteps(structs.ssAfter.rec);
  TS = TS(abs(TS-round(TS))<eps(1e3));
  drawSimulation(p, structs.ssAfter.rec, TS);
  close all hidden;
  pth = [ndir filesep indstr '.' nfile];
  makeVideo([pth '.%04d.png'], [pth '.mp4']);
  toc
end


function fun = handleImageFUN(newsize)
  fun = @(varargin)handleImage(newsize, varargin{:});

function handleImage(newsize, filename, img, numFrame)
  img = uint8(sum(double(img), 3)/3);
  img = imresize(img, 'OutputSize', newsize, 'Method', 'lanczos2', 'Antialiasing', true);
  imwrite(img, sprintf(filename, numFrame), 'png');
%@(filename, img, numFrame) imwrite(imresize(uint8(sum(double(img), 3)/3), 'OutputSize', [600 600], 'Method', 'lanczos2', 'Antialiasing', true), sprintf(filename, numFrame), 'png');