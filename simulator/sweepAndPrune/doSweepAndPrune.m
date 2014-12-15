function [candidatos sap] = doSweepAndPrune(pos, rad, sap)

sets   = sap.sets;
mode   = sap.mode;
filter = sap.filter;
idxsx  = sap.idxsx;
idxsy  = sap.idxsy;

[np nd] = size(pos);

if nd<3
  
  %construct limits
  limits             = [pos; pos];
  limits(1:np,1)     = limits(1:np,    1) + rad;
  limits(np+1:end,1) = limits(np+1:end,1) - rad;
  limits(1:np,2)     = limits(1:np,    2) + rad;
  limits(np+1:end,2) = limits(np+1:end,2) - rad;
  %limits             = sapLimits(pos, rad);

  %sweep in x
  [tmp tmp]          = sort(limits(idxsx,1));
  idxsx              = idxsx(tmp);
  %sweep in y
  [tmp tmp]          = sort(limits(idxsy,2));
  idxsy              = idxsy(tmp);

  %prune
  candidatos         = sweepAndPruneWithSetsCAll(sets, mode, filter, idxsx, idxsy);

  sap.idxsx        = idxsx;
  sap.idxsy        = idxsy;
  
else
  idxsz = sap.idxsz;

  %construct limits
  limits             = [pos; pos];
  limits(1:np,1)     = limits(1:np,    1) + rad;
  limits(np+1:end,1) = limits(np+1:end,1) - rad;
  limits(1:np,2)     = limits(1:np,    2) + rad;
  limits(np+1:end,2) = limits(np+1:end,2) - rad;
  limits(1:np,3)     = limits(1:np,    3) + rad;
  limits(np+1:end,3) = limits(np+1:end,3) - rad;
  %limits             = sapLimits(pos, rad);

  %sweep in x
  [tmp tmp]          = sort(limits(idxsx,1));
  idxsx              = idxsx(tmp);
  %sweep in y
  [tmp tmp]          = sort(limits(idxsy,2));
  idxsy              = idxsy(tmp);
  %sweep in y
  [tmp tmp]          = sort(limits(idxsz,2));
  idxsz              = idxsz(tmp);

  %prune
  candidatos         = sweepAndPruneWithSetsCAll(sets, mode, filter, idxsx, idxsy, idxsz);

  sap.idxsx        = idxsx;
  sap.idxsy        = idxsy;
  sap.idxsz        = idxsz;

end

