function [rdyn r] = clearOldPerturbations(ss, rdyn, r, t)
ra = rdyn.radapt;
convergence = (ss.k(ra).*abs(rdyn.r0(ra)-r(ra))./rdyn.r0(ra)) < 0.01;
if any(convergence)
  ra = find(ra);
  ra = ra(convergence);
  rdyn.radapt(ra) = false;
  r(ra) = rdyn.r0(ra);
  rdyn.rchg = any(rdyn.radapt);
end
end
