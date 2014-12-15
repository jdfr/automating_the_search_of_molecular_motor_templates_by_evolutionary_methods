%to be used before simulating if the the number of points or spring
%changes, or the connectivity of any spring changes
function ss = makeSpringEndsMatrix(ss)
  springEnds = ss.springEnds;
  %remakes the matrix of spring ends. It is a sparse matrix
  %numspringsXnumpoints; its (i,j)-th position is a placeholder for
  %quantities related to the i-th spring and the j-th point,
  %provided that the j-th point is one of the ends of the i-th spring
  ons = ones(ss.nsprings,1);
  ss.springsMatrix = sparse(springEnds(:), repmat(1:ss.nsprings,1,2), [ons -ons], ss.npoints, ss.nsprings);
  
%   ss.pSpring(springEnds(:,1)) = 1:size(springEnds,1);
%   ss.pSpring(springEnds(:,2)) = 1:size(springEnds,1);

%   spneighs = cell(size(ss.r));
%   for k=1:size(springEnds,1)
%     spm         = (springEnds(k,1)==springEnds) | (springEnds(k,2)==springEnds);
%     spm         = spm(:,1) | spm(:,2);
%     spm(k)      = false;
%     spm         = find(spm);
%     spneighs{k} = [repmat(k, size(spm)), spm];
%   end
%   spneighs = vertcat(spneighs{:});
%   ss.springsNeighs = sparse(spneighs(:,1), spneighs(:,2), true, numel(ss.r), numel(ss.r));
  
%   springsNeighs = abs(ss.springsMatrix);
%   springsNeighs = logical(springsNeighs'*springsNeighs);
%   springsNeighs(1:(size(springsNeighs,1)+1):numel(springsNeighs)) = false;
%   
%   ss.springsNeighs = springsNeighs;
  
%   if any(any( springsNeighs ~= ss.springsNeighs ))
%     error('This should never happen!!!');
%   end
end
