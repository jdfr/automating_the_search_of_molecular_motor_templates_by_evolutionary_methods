function [H, eigVecsH, eigValsH, K, eigVecsK, eigValsK] = testANM(ss, useUniformK)
%function [H, eigVecsH, eigValsH, K, eigVecsK, eigValsK, TM, eigVecsTM, eigValsTM] = testANM(ss, useUniformK)

if ~exist('useUniformK', 'var')
  useUniformK = true;
end

if exist('ss', 'var')
  ss = fuseRepeatedSprings(ss, false);
  df = ss.pos(ss.springEnds(:,1),:)-ss.pos(ss.springEnds(:,2),:);
  ss.r = realsqrt(sum(df.*df, 2));
else
  ss = makecosa;
end

if useUniformK
  ss.k(1:end) = median(ss.k);
end

[H,  eigVecsH,  eigValsH]  = ANMdecomposition(ss);
[K,  eigVecsK,  eigValsK]  = GNMdecomposition(ss);
%[TM, eigVecsTM, eigValsTM] =  TMdecomposition(ss);

pos = ss.pos;

np = size(pos,1);

facH = 5;10;
facK = 10;

facesK = 15;

axiswindow = getAxisWindow(pos, 0.1);

p_def = {'axisWindow', axiswindow, 'axisEqual', true, 'axisSquare', false};%, 'circleFaces', 50}; 
%p_def = {'axisWindow', [-10 60 -30 40], 'axisEqual', false, 'axisSquare', true};%, 'circleFaces', 50}; 
p_fun = @EigenPlotter;

bothsigns = true;

p = p_fun();
set(p, p_def{:});
rec = makeNewRecorder(ss, 'MemRecorder', 10000, 10000, 100, 100);
ss.t = 0;
for V=1:np;
  ss.stick.eigTxt = 'base structure';
  ss.stick.eigVec1 = [];
  ss.stick.eigVectorss = [];
%   ss.stick.eigVector = [];
  ss.pos = pos;
  rec = recordAllState(rec, ss);
  ss.t = ss.t+1;
  ss.stick.eigVec1 = facH*reshape(eigVecsH(:,V), size(ss.pos, 2), [])'; 
%   ss.stick.eigVector = realsqrt(sum(realpow(ss.stick.eigVec1/facH, 2), 2));
  ss.pos = ss.pos + ss.stick.eigVec1;
%   ss.stick.eigVectorss = realsqrt(sum(realpow(ss.stick.eigVec1(ss.springEnds(:,1),:)-ss.stick.eigVec1(ss.springEnds(:,2),:), 2), 2));
  ss.stick.eigVectorss = ss.stick.eigVec1/facH;
  ss.stick.eigTxt = sprintf('warped by eigenvalue #%d: %s, flavor a', V, mat2str(eigValsH(V,V)));
  rec = recordAllState(rec, ss);
  ss.t = ss.t+1;
  
  if ~bothsigns
    continue;
  end
  
  ss.stick.eigTxt = 'base structure';
  ss.stick.eigVec1 = [];
  ss.stick.eigVectorss = [];
%   ss.stick.eigVector = [];
  ss.pos = pos;
  rec = recordAllState(rec, ss);
  ss.t = ss.t+1;
  ss.stick.eigVec1 = -facH*reshape(eigVecsH(:,V), size(ss.pos, 2), [])'; 
%   ss.stick.eigVector = realsqrt(sum(realpow(ss.stick.eigVec1/facH, 2), 2));
  ss.pos = ss.pos + ss.stick.eigVec1;
%   ss.stick.eigVectorss = realsqrt(sum(realpow(ss.stick.eigVec1(ss.springEnds(:,1),:)-ss.stick.eigVec1(ss.springEnds(:,2),:), 2), 2));
  ss.stick.eigVectorss = ss.stick.eigVec1/facH;
  ss.stick.eigTxt = sprintf('warped by eigenvalue #%d: %s, flavor b', V, mat2str(eigValsH(V,V)));
  rec = recordAllState(rec, ss);
  ss.t = ss.t+1;
end
nav = navigateSimulation(p, rec);
set(gcf, 'name', 'anisotropic vibrational modes (ANM)');

ss.stick = rmfield(ss.stick, 'eigVectorss');

p = p_fun();
set(p, p_def{:});
p.change = true;
rec = makeNewRecorder(ss, 'MemRecorder', 10000, 10000, 100, 100);
ss.stick.eigTxt = 'base structure';
ss.stick.eigVec1 = [];
ss.stick.eigVec2 = [];
% ss.stick.eigVector = [];
ss.pos = pos;
ss.t = 0;
rec = recordAllState(rec, ss);
ss.t = ss.t+1;
if size(ss.pos,2)<3
  angles   = linspace(0, 2*pi, facesK+1)';
  senos    = sin(angles)*facK;
  cosenos  = cos(angles)*facK;
  supermat = repmat({[[cosenos, senos]; [nan nan]]}, size(ss.pos, 1), 1);
end
for V=1:np;
  eigVecsHMags = reshape(eigVecsH(:,V), size(ss.pos, 2), [])';
  eigVecsHMags = realsqrt(sum(realpow(eigVecsHMags, 2), 2));
  if size(ss.pos,2)<3
    mimats = supermat;
    for z=1:size(ss.pos)
      mimats{z}      = mimats{z}*eigVecsHMags(z);
      mimats{z}      = bsxfun(@plus, mimats{z}, ss.pos(z,:));%[mimats{z}(:,1)+ss.pos(z,1), mimats{z}(:,2)+ss.pos(z,2)];
    end
    ss.stick.eigVec1 = vertcat(mimats{:});
  else
    ss.stick.eigVec1 = eigVecsHMags;
  end
  ss.stick.eigVector = eigVecsHMags;
  ss.stick.eigTxt = sprintf('warped by eigenvalue #%d: %s', V, mat2str(eigValsH(V,V)));
  rec = recordAllState(rec, ss);
  ss.t = ss.t+1;
end
nav = navigateSimulation(p, rec);
set(gcf, 'name', 'anisotropic vibrational magnitudes (ANM)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = p_fun();
set(p, p_def{:});
p.change = true;
rec = makeNewRecorder(ss, 'MemRecorder', 10000, 10000, 100, 100);
ss.stick.eigTxt = 'base structure';
ss.stick.eigVec1 = [];
ss.stick.eigVec2 = [];
ss.stick.eigVector = [];
ss.pos = pos;
ss.t = 0;
rec = recordAllState(rec, ss);
ss.t = ss.t+1;
if size(ss.pos,2)<3
  angles   = linspace(0, 2*pi, facesK+1)';
  senos    = sin(angles)*facK;
  cosenos  = cos(angles)*facK;
  supermat = repmat({[[cosenos, senos]; [nan nan]]}, size(ss.pos, 1), 1);
end
for V=1:np;
  signo1 = eigVecsK(:,V)>0;
  if size(ss.pos,2)<3
    mimats = {supermat(signo1); supermat(~signo1)};
    ks = [1 1];
    for z=1:size(ss.pos)
      kz = (~signo1(z))+1;
      ksz = ks(kz);
      mimats{kz}{ksz}      = mimats{kz}{ksz}*eigVecsK(z,V);
      mimats{kz}{ksz}      = bsxfun(@plus, mimats{z}, ss.pos(z,:));%[mimats{kz}{ksz}(:,1)+ss.pos(z,1), mimats{kz}{ksz}(:,2)+ss.pos(z,2)];
      ks(kz) = ks(kz)+1;
    end
    ss.stick.eigVec1 = vertcat(mimats{1}{:});
    ss.stick.eigVec2 = vertcat(mimats{2}{:});
  else
    ss.stick.eigVec1 = eigVecsK(signo1,V);
    ss.stick.eigVec2 = eigVecsK(~signo1,V);
  end
  ss.stick.eigVector = eigVecsK(:,V);
  ss.stick.eigTxt = sprintf('warped by eigenvalue #%d: %s', V, mat2str(eigValsK(V,V)));
  rec = recordAllState(rec, ss);
  ss.t = ss.t+1;
end
nav = navigateSimulation(p, rec);
set(gcf, 'name', 'isotropic vibrational magnitudes (GNM)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p = p_fun();
% set(p, p_def{:});
% rec = makeNewRecorder(ss, 'MemRecorder', 10000, 10000, 100, 100);
% ss.t = 0;
% for V=1:np;
%   ss.stick.eigTxt = 'base structure';
%   ss.stick.eigVec1 = [];
%   ss.pos = pos;
%   rec = recordAllState(rec, ss);
%   ss.t = ss.t+1;
%   ss.stick.eigVec1 = facH*reshape(eigVecsTM(:,V), 2, [])'; 
%   ss.pos = ss.pos + ss.stick.eigVec1;
%   ss.stick.eigTxt = sprintf('warped by eigenvalue #%d: %s, flavor a', V, mat2str(eigValsTM(V,V)));
%   rec = recordAllState(rec, ss);
%   ss.t = ss.t+1;
%   
%   if ~bothsigns
%     continue;
%   end
%   
%   ss.stick.eigTxt = 'base structure';
%   ss.stick.eigVec1 = [];
% %   ss.stick.eigVector = [];
%   ss.pos = pos;
%   rec = recordAllState(rec, ss);
%   ss.t = ss.t+1;
%   ss.stick.eigVec1 = -facH*reshape(eigVecsTM(:,V), 2, [])'; 
%   ss.pos = ss.pos + ss.stick.eigVec1;
%   ss.stick.eigTxt = sprintf('warped by eigenvalue #%d: %s, flavor b', V, mat2str(eigValsTM(V,V)));
%   rec = recordAllState(rec, ss);
%   ss.t = ss.t+1;
% end
% nav = navigateSimulation(p, rec);
% set(gcf, 'name', 'anisotropic vibrational modes (TM)');
% 
% p = p_fun();
% set(p, p_def{:});
% p.change = true;
% rec = makeNewRecorder(ss, 'MemRecorder', 10000, 10000, 100, 100);
% ss.stick.eigTxt = 'base structure';
% ss.stick.eigVec1 = [];
% ss.stick.eigVec2 = [];
% % ss.stick.eigVector = [];
% ss.pos = pos;
% ss.t = 0;
% rec = recordAllState(rec, ss);
% ss.t = ss.t+1;
% angles   = linspace(0, 2*pi, facesK+1)';
% senos    = sin(angles)*facK;
% cosenos  = cos(angles)*facK;
% supermat = repmat({[[cosenos, senos]; [nan nan]]}, size(ss.pos, 1), 1);
% for V=1:np;
%   eigVecsTMMags = reshape(eigVecsTM(:,V), 2, [])';
%   eigVecsTMMags = realsqrt(sum(realpow(eigVecsTMMags, 2), 2));
%   mimats = supermat;
%   for z=1:size(ss.pos)
%     mimats{z}      = mimats{z}*eigVecsTMMags(z);
%     mimats{z}      = [mimats{z}(:,1)+ss.pos(z,1), mimats{z}(:,2)+ss.pos(z,2)];
%   end
%   ss.stick.eigVec1 = vertcat(mimats{:});
% %   ss.stick.eigVector = eigVecsHMags;
%   ss.stick.eigTxt = sprintf('warped by eigenvalue #%d: %s', V, mat2str(eigValsTM(V,V)));
%   rec = recordAllState(rec, ss);
%   ss.t = ss.t+1;
% end
% nav = navigateSimulation(p, rec);
% set(gcf, 'name', 'anisotropic vibrational magnitudes (TM)');

end

function ss = makecosa

conf.organism.geneTypes = [];
conf.organism.useShift = false;
conf.organism.useAdaptativeSprings = false;
conf.organism.atomsShareSprings = true;
conf.organism.atomConf = [];
conf.organism.tickWidth = 0.01;
conf.organism.u = 5;
conf.organism.rows = 2;
conf.organism.cols = 10;
conf.organism.xdesp = 5;
conf.organism.ydesp = 5;
conf.organism.springFixation.springToFix = 1;
conf.organism.springFixation.fixedForceDensity = 0;
conf.organism.springFixation.springSpec = springSpec('F', 1, 1);

conf = atomGridSquare(conf);

ss = conf.ssBase;

end

function [kirchoff, eiVecs, eiVals] = GNMdecomposition(ss)

ks = ss.k;

pos = ss.pos;

springEnds = ss.springEnds;

N = size(pos,1);

conns = sparse(springEnds(:), reshape(springEnds(:,[2,1]), [], 1), 1, N, N);

sumcons = sum(conns);

kirchoff = full(sparse(1:N, 1:N, sumcons, N,N) - conns);

[eiVecs, eiVals] = eig(kirchoff);%, 'nobalance');

end

function [lambda, eiVecs, eiVals] = TMdecomposition(ss)

ks = ss.k;

pos = ss.pos;

springEnds = ss.springEnds;

N = size(pos,1);

diffs = pos(springEnds(:,1),:)-pos(springEnds(:,2),:);

dists2 = sum(realpow(diffs, 2), 2);

lambda = zeros(2*N);

for k=1:size(springEnds, 1)
  rix = springEnds(k,1)*2-1;
  riy = springEnds(k,1)*2;
  rjx = springEnds(k,2)*2-1;
  rjy = springEnds(k,2)*2;
  diffsi     =  diffs(k,:);
  diffsj     = -diffsi;
  factorsij  = [rix         riy        rjx         rjy];
  sequenceij = [    diffsi                  diffsj    ]/dists2(k);
  factorsji  = [rjx         rjy        rix         riy];
  sequenceji = -sequenceij;
  lambda(rix, factorsij) =  sequenceij*diffsi(1);
  lambda(riy, factorsij) =  sequenceij*diffsi(2);
  lambda(rjx, factorsji) =  sequenceji*diffsj(1);
  lambda(rjy, factorsji) =  sequenceji*diffsj(2);
end

[eiVecs, eiVals] = eig(lambda);%, 'nobalance');

end


