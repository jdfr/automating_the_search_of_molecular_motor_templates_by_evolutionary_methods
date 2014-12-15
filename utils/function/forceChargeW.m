function out = forceChargeW(dims, npoints, pos, e, chargeParams)

vectors = bsxfun(@minus, reshape(cjto,[npoints 1 dims]), reshape(cjto,[1 npoints dims]));
dists = realsqrt(sum(diffs.^2,3));
vectors = bsxfun(@rdivide, vectors, dists);
intensity = nWeibull(chargeParams.wl, chargeParams.wk, dists);

force

forces = 

x1=;
x2=; xx= xx2 = xx.^2; xx3 = sum(xx2,3); xx4=realsqrt(xx3)


forces = zeros(size(pos));

%find distance matrix
dists = distanceMatrixSquared(pos);
%find intensity of the force for each pair of particles
intensity = chargeParams.K*nWeibull(chargeParams.wl, chargeParams.wk, dists);
%multiply each intensity by the particle's charge to compute the field
%intensity
field = bsxfun(@times, intensity, e);
field = bsxfun(@times, field, e');
%add field intensities for each particle and multiply by its charge.
out = sum(field)'.*e;

%unfolded code
out = sum(bsxfun(@times, chargeParams.K*nWeibull(chargeParams.wl, chargeParams.wk, distanceMatrix(pos)), e))'.*e;