function poph = addSpectralGapToPoph(poph)

ss = poph.params.evparams.ss;
if isfield(ss, 'atomPoints')
  ss = removeAtoms(ss, (1:size(ss.atomPoints, 1))', true, true);
end

if (size(poph.params.evparams.ss.pos, 2)==3)
  resType = struct('unit',       {{[]}}, ...
                   'getRasters', {@(res)res.genomes});
else
  resType = 2;
end
poph = postProcessPoph(poph, true, @(varargin)addSpectralGap(ss, varargin{:}), resType);

function poph = addSpectralGap(ss, poph, generation, thisGen, indexNewG, rasters) %#ok<INUSL>

if isempty(generation)
  poph.fstSpecGap     = zeros(size(poph.fitness));
  poph.scdSpecGap     = zeros(size(poph.fitness));
  poph.trdSpecGap     = zeros(size(poph.fitness));
  poph.maxSpecGap     = zeros(size(poph.fitness));
  poph.idxMaxSpecGap  = zeros(size(poph.fitness));
  poph.spectra        =  cell(size(poph.fitness));
  poph.minAutoval     = zeros(size(poph.fitness));
%   poph.collectivityA  =  cell(size(poph.fitness));
%   poph.collectivityB  =  cell(size(poph.fitness));
  return
end

for k=1:numel(thisGen)
  g = thisGen(k);
  
  if (size(poph.params.evparams.ss.pos, 2)==3) && iscell(rasters) && isa(rasters{k}, 'uint16')
    pos = convertCompressedPos(rasters{k}, poph.params.evparams, true);
    if size(pos,1)~=poph.params.genome.fold.numPoints
      error('This should never happen!!!');
      %pos = pos(poph.params.genome.fold.numPoints,:);
    end
    dst = poph.params.evparams.genome.L0;
    [sp1 sp2] = find(triu(distanceMatrixSquared(pos)<=(dst*dst), 1));
    r         = realsqrt(sum(realpow(pos(sp1,:)-pos(sp2,:), 2), 2));
    ss        = addSprings(addPoints(poph.params.evparams.ss, pos), [sp1 sp2], poph.params.evparams.walker.springK, r);
    [autovecs autovecs autovals] = ANMdecomposition(ss);
  else
    [autovecs autovecs autovals] = ANMdecomposition(combineBasicSystems(ss, rasters(k).dev, [ss.pointVars, ss.springVars]));
  end
  
  autovals = diag(autovals);
  
%         eiv     = st.stick.eigVectorss;
%         vec2    = realpow(eiv, 2);
%         vec2s   = vec2/sum(vec2);
%         ka      = 1/size(eiv, 1)*exp(-sum(vec2s.*log(vec2s)));
%         eivsp   = eiv(st.springEnds(:,1),:)-eiv(st.springEnds(:,2),:);
%         vecsp2  = realpow(eivsp, 2);
%         vecsp2s = vecsp2/sum(vecsp2);
%         kb      = 1/size(eivsp, 1)*exp(-sum(vecsp2s.*log(vecsp2s)));
%   collectivity = 
  
  bigEnough      = (autovals>1e-12) ;%| ((autovals>1e-4) & (collectivityB>0.1));
  bigEnough(1:3) = false;

  autovals       = autovals(bigEnough);
  minautoval     = autovals(1);

  autovals       = log10(autovals/autovals(1));
  
  gaps           = diff(autovals);
  
  poph.spectra{g}                            = autovals;
  poph.minAutoval(g)                         = minautoval;
  poph.fstSpecGap(g)                         = autovals(2);
  poph.scdSpecGap(g)                         = autovals(3);
  poph.trdSpecGap(g)                         = autovals(4);
  [poph.maxSpecGap(g) poph.idxMaxSpecGap(g)] = max(gaps);
end
