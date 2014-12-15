function [ok spectrum]         = validChain(pos, conns, sp1, sp2, r, evparams)

chainr                         = realsqrt(sum(realpow(pos(1:end-1,:)-pos(2:end,:), 2), 2));
ok                             = all((sum(conns)+sum(conns,2)')>=evparams.genome.minConns) && ...
                                 all(chainr<=evparams.genome.maxChainDist) && ...
                                 all(r>=evparams.genome.minChainDist);

if ok && isfield(evparams.genome, 'gap') && evparams.genome.gap.useIt
  rotCutoff                    = evparams.spectralGap.rotCutoff;
  gapToRecord                  = evparams.genome.gap.toRecord;
  gapInterval                  = evparams.genome.gap.interval;
  [spectrum spectrum spectrum] = ANMdecomposition({pos, conns, sp1, sp2}, false);
  spectrum                     = diag(spectrum);
  bigEnough                    = spectrum>rotCutoff;
  bigEnough(1:3)               = false; %these ALWAYS ARE DISPLACEMENTS/ROTATIONS
  spectrumS                    = spectrum(bigEnough);
  spectrumS                    = log10(spectrumS/spectrumS(1));
  gap                          = spectrumS(gapToRecord+1);
  ok                           = (gap>=gapInterval(1)) && (gap<=gapInterval(2));
else
  spectrum                     = [];
end