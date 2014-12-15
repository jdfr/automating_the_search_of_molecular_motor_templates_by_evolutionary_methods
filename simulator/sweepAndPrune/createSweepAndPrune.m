function sap = createSweepAndPrune(pos, rad, sets, mode, filter)

np = size(pos,1);

if nargin<3
  sets = [];%zeros(np, 1, 'uint8');
else
  sets = uint8(sets);
end

if nargin<4
  mode = true;
else
  mode = logical(mode);
end

if nargin<5
  filter = [];
else
  filter = int8(filter);
end

%construct limits
limits               = [pos; pos];
limits(1:np,1)       = limits(1:np,    1) + rad;
limits(np+1:end,1)   = limits(np+1:end,1) - rad;
limits(1:np,2)       = limits(1:np,    2) + rad;
limits(np+1:end,2)   = limits(np+1:end,2) - rad;

%sweep in x
[tmpx tmpx]          = sort(limits(:,1));
%sweep in y
[tmpy tmpy]          = sort(limits(:,2));

if size(pos,2)<3
  sap                = struct('sets', sets, 'mode', mode, 'filter', filter, 'idxsx', int32(tmpx), 'idxsy', int32(tmpy));
else
  limits(1:np,3)     = limits(1:np,    3) + rad;
  limits(np+1:end,3) = limits(np+1:end,3) - rad;
  [tmpz tmpz]        = sort(limits(:,3));
  sap                = struct('sets', sets, 'mode', mode, 'filter', filter, 'idxsx', int32(tmpx), 'idxsy', int32(tmpy), 'idxsz', int32(tmpz));
end


