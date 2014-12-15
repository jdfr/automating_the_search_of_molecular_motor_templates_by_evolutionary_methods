function [H, eiVecs, eiVals] = ANMdecomposition(ss, doFusion)

if nargin<2
  doFusion = true;
end

if isstruct(ss)
  if doFusion
    ss = fuseRepeatedSprings(ss, false);
  end

  pos = ss.pos;

  springEnds = ss.springEnds;

  N = size(pos,1);

  conns = sparse(springEnds(:), reshape(springEnds(:,[2,1]), [], 1), true, N, N);
elseif iscell(ss)
  pos = ss{1};
  N = size(pos,1);
  springEnds = [ss{3} ss{4}];
  %conns = ss{2};
  conns = sparse(springEnds(:), reshape(springEnds(:,[2,1]), [], 1), true, N, N);
else
  error('First argument not understood!!!!');
end

sumcons = sum(conns);

diffs = pos(springEnds(:, 2),:)-pos(springEnds(:,1),:);

dists = realsqrt(sum(realpow(diffs, 2), 2));

nd = size(pos,2);
d3 = nd>2;

H = repmat({zeros(nd, nd)}, N, N);
dffs = zeros(nd,nd,size(springEnds, 1));

dffXX = diffs(:,1).*diffs(:,1);
dffXY = diffs(:,1).*diffs(:,2);
dffYY = diffs(:,2).*diffs(:,2);
dffs(1,1,:) = dffXX(:);
dffs(1,2,:) = dffXY(:);
dffs(2,1,:) = dffXY(:);
dffs(2,2,:) = dffYY(:);

if d3
  dffZZ = diffs(:,3).*diffs(:,3);
  dffXZ = diffs(:,1).*diffs(:,3);
  dffYZ = diffs(:,2).*diffs(:,3);
  dffs(1,3,:) = dffXZ(:);
  dffs(2,3,:) = dffYZ(:);
  dffs(3,1,:) = dffXZ(:);
  dffs(3,2,:) = dffYZ(:);
  dffs(3,3,:) = dffZZ(:);
end

dffs = -bsxfun(@rdivide, dffs, reshape(dists, 1, 1, []));

for k=1:size(springEnds, 1)
  a = springEnds(k, 1);
  b = springEnds(k, 2);
%   if dists(k)>0
    H{a,b} = dffs(:,:,k);
    H{b,a} = H{a,b};
%   else
%     warning('ANMDecomposition:zeroLength', 'There is a connection between points %d and %d, but their locations are identical. The spring will be omitted from the analysis, but this is unsafe', a, b);
%   end
end

for a=1:N
  H{a,a} = -sum(reshape(horzcat(H{a, conns(a,:)}), nd, nd, sumcons(a)), 3);
end
  
H = cell2mat(H);

if all(dists>0)%(~(any(isnan(H(:))))) && (~(any(isinf(H(:)))))

  [eiVecs, eiVals] = eig(H);%, 'nobalance');
  
else
  
  eiVecs = nan;
  eiVals = nan;
  
end

end
