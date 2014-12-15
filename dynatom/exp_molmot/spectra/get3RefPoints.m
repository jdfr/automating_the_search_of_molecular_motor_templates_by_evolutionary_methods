function points = get3RefPoints(ss)

threshold           = 1e-12;

[H, eiVecs, eiVals] = ANMdecomposition(ss, false);

eiVals              = diag(eiVals);
modes               = eiVals>threshold;
eiVecs              = eiVecs(:,modes);
eiVals              = eiVals(modes);
if any(diff(eiVals)<(-threshold))
  error('The normal modes are not ordered!!!!!');
end

mode1               = reshape(eiVecs(:,1), 3, [])';
mode2               = reshape(eiVecs(:,2), 3, [])';

d0                  = distanceMatrix(ss.pos);
dd01                = abs(distanceMatrix(ss.pos+mode1)-d0);

[idx idx]           = max(dd01(:));
[p1 p2]             = ind2sub(size(dd01), idx);
[p3 p3]             = max(abs(distanceMatrix(ss.pos(p1,:)+mode2(p1,:), ss.pos+mode2)-d0(p1,:)));
if p2==p3
  [p3 p3]           = max(abs(distanceMatrix(ss.pos(p2,:)+mode2(p2,:), ss.pos+mode2)-d0(p2,:)));
end
points              = [p1; p2; p3];

