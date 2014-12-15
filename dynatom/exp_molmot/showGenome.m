function ss = showGenome(genome, params)

rndSeedDev  = [];
useKeyboard = true;
recDummy    = [];


if isstruct(genome)
  ss               = genome;
  legerror         = [];
else
  evparams         = transformCorrectParams(params.evparams);
  ss               = manipulateGenome3D(genome, rndSeedDev, evparams, useKeyboard, recDummy);
  [ss legerror ss] = prepareLegsSimulation(ss, evparams);
end

if isempty(legerror)
  for k=1:numel(ss.legs)
    leg = ss.legs(k);
        bnds                   = leg.ABCbinding;
        ABCw                   = leg.ABCweight;
        ABCb                   = ss.pos(bnds,:);
        ABCv                   = ss.vel(bnds,:);
        atpp                   = sum(bsxfun(@times, ABCb, ABCw))/sum(ABCw);
        atpv                   = sum(bsxfun(@times, ABCv, ABCw))/sum(ABCw);
        [ss ss.legs(k).atp]    = addPoints(ss, atpp, atpv);
        if leg.dynLengths
          dists                = realsqrt(sum(realpow(bsxfun(@minus, ABCb, atpp), 2), 2));
          rs                   = dists .* leg.dynLengthFactor;
        else
          rs                   = leg.atpLengths;
        end
        [ss ss.legs(k).atpb]   = addSprings(ss, [repmat(size(ss.pos,1), size(rs)), bnds], leg.atpK, rs);
  end
  if isfield(ss, 'updateSAP')
    ss.sap = ss.updateSAP(ss);
  end
end