function [genome cad maskNew] = mutateFold(idx, G, params, P, idxInP, maskNew)  %#ok<INUSD>

genome                     = G.genome{idx};

%{0.5, 0:10, poisspdf(0:10, 0.5)}
numMuts                    = pickOption(params.genome.poisson{[2 3]}, 1);

if numMuts==0
  maskNew                  = false;
  cad                      = '';
  return;
end

pos                        = convertCompressedPos(genome, params.evparams);
n                          = size(pos,1);

dst                        = params.evparams.genome.L0;

dstsq                      = distanceMatrixSquared(pos);
%dstsq((n+1):(n+1):(2*n))   = 0; %equivalent to dstsq(sub2ind(size(dstsq), 1:n-1, 2:n)) = 0;
%[sp1 sp2]                  = find(triu(dstsq<=(dst*dst), 1));
%triu(m,2): we get rid of the lower triangle, as well as the 
%first diagonal after the main one, since it contains the
%backbone links of the form (n, n+1), which we DO NOT want to change if we
%want to model molecules
[sp1 sp2]                  = find(triu(dstsq<=(dst*dst), 2)); 
r                          = realsqrt(dstsq(sp1+(sp2-1)*n));

springs                    = pickOption(1:numel(sp1), r, numMuts);
factors                    = affinTransform(rand(numMuts,1), params.genome.mutationFactorRange, [0 1]);

mutation                   = [sp1(springs), sp2(springs), factors];
genome                     = {genome, mutation};
maskNew                    = true;
cad                        = any2str(mutation);
cad(cad==' ')              = ',';
