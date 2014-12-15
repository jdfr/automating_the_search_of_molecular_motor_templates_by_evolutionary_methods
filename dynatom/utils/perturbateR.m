%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function rdyn = perturbateR(rdyn, r, t, indexes, newr)
rdyn.radapt(indexes) = true;
rdyn.rchg = true;
rdyn.r0(indexes) = newr;
rdyn.rtb(indexes) = t;
rdyn.rspan(indexes) = r(indexes)-newr;
end
