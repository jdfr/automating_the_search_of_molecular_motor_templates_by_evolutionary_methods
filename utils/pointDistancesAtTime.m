function [I1 ] = pointDistancesAtTime(rec, ts)

%distsMatrix = [];
I1 = zeros(numel(ts),1);
index = 1;
playSimulation(rec, ts, @distanceRecorder);


  function distanceRecorder(t, ss)
 %   distsMatrix = distanceMatrix(single(ss.pos));
%     quad = ss.pos.^2;
%     mu00 = size(ss.pos,1);
%     mu02 = sum(quad(:,1).*ss.m);
%     mu20 = sum(quad(:,2).*ss.m);
%     nu02 = mu02/(mu00.^(1+(1+1)/2));
%     nu20 = mu20/(mu00.^(1+(1+1)/2));
%     I1   = nu02+nu20;
    I1(index)   = sum(sum(bsxfun(@times, ss.pos.^2, ss.m)))/(size(ss.pos,1).^2);
    index=index+1;
    if (mod(index,1000)==0)
      fprintf('tardaaaaa %d\n', index);
    end
  end
end