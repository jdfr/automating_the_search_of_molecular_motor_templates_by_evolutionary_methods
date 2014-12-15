function [dpos, dvel] = magicallyMovedPoints(ss, dpos, dvel, t)
%Given several triplets of points (M,1,2), control the points M's kinematics:  
%   -before the interval, their movements are stuck to points 1's
%   -after  the interval, their movements are stuck to points 2's
%   -during the interval, they smoothly move along the lines from points 1
%    to points 2 
%
%dpos==vel, dvel==acel

if ~isempty(ss.pShifting_M)
  %get data
  pos = ss.pos;
  pM = ss.pShifting_M;
  if isfield(ss, 'shiftedPoints')
    shifted = ss.shiftedPoints;
    shiftingMatrix = ss.shiftingMatrix;
  else
    [shifted nevermind onlyShiftedIndexes] = unique(ss.pShifting_M);
    %this sparse matrix has 1 in the (i,j)-th position, meaning that the
    %i-th shifted point is in the j-th shifting
    % ss.shiftingMatrix = sparse(ss.pShifting_M, 1:numel(ss.pShifting_M), 1, size(ss.pos,1), numel(ss.pShifting_M));
    shiftingMatrix = sparse(onlyShiftedIndexes, 1:numel(ss.pShifting_M), 1, numel(shifted), numel(ss.pShifting_M));
  end
  p1 = ss.pShifting_1;
  p2 = ss.pShifting_2;
  T0 = ss.pShifting_T0;
  Tf = ss.pShifting_Tf;
  TD = ss.pShifting_TD;
  %calculate linear factors: 1 at t==T0, 0 at t==Tf
  factors  = (t-Tf)./TD;
% PREVIOUSLY, THE CODE WAS PREPARED TO HANDLE TRIPLETS OUTSIDE ITS LIFE
% CYCLE. HOWEVER, THE CODE WAS UNABLE TO HANDLE POINTS ATTACHED TO SEVERAL
% TRIPLETS. SINCE IT IS QUITE IMPORTANT TO ACHIEVE THIS, WE DROP THE
% ABILITY TO HANDLE TRIPLETS OUTSIDE ITS LIFE CYCLE (TOO EXPENSIVE TO DO
% IT PROPERLY). SO, WE ***RELY*** ON THE SYSTEM TO HALT THE SIMULATION
% WHENEVER A TRIPLET HAS FINISHED ITS LIFE CYCLE, TO WIPE IT OUT. THIS ALSO
% IS BETTER, SINCE THE MAKING OF THE INCIDENCE (SPARSE) MATRIX CAN BE
% OUTSOURCED
%   %decide the status for each point
%   moving   = (factors>0) & (factors<1);
%   pre      = factors>1;
%   post     = factors<0;
%   moving   = (~pre) & (~post);
%   %pM points in 'pre'  state: stuck to p1
%     selM = pM(pre);
%   if ~isempty(selM)
%     sel1 = p1(pre);
%     dpos(selM,:) = dpos(sel1,:);
%     dvel(selM,:) = dvel(sel1,:);
%   end
%   %pM points in 'post' state: stuck to p2
%     selM = pM(post);
%   if ~isempty(selM)
%     sel2 = p2(post);
%     dpos(selM,:) = dpos(sel2,:);
%     dvel(selM,:) = dvel(sel2,:);
%   end
%   %pM points in 'moving' state:
%     selM = pM(moving);
%   if ~isempty(selM)
%     sel2 = p2(moving);
%     sel1 = p1(moving);
% %     %tune the kinematics of colliding points (the influence is weighted by
% %     %mass)
% %     V   = ss.pShifting_V(moving,:);
% %     mM  = ss.m(selM);
% %     m2  = ss.m(sel2);
% %     mM2 = mM+m2;
% %     dpos_avg = bsxfun(@rdivide, ...
% %                         bsxfun(@times, dpos(selM,:), mM) ... 
% %                       + bsxfun(@times, dpos(sel2,:), m2), ... 
% %                       mM2);
% %     dvel_avg = bsxfun(@rdivide, ...
% %                         bsxfun(@times, dvel(selM,:), mM) ... 
% %                       + bsxfun(@times, dvel(sel2,:), m2), ... 
% %                       mM2);
% %     dpos(selM,:) = dpos_avg + V;
% %     dpos(sel2,:) = dpos_avg - V;
% %     dvel(selM,:) = dvel_avg;
% %     dvel(sel2,:) = dvel_avg;
%     %TODO: this seems to be unrealistic. Clearly, the behaviour must not be
%     %so asymmetric
%     %calculate velocities to move along the line between p1 and p2, during
%     %the interval [T0, Tf]. This is simply the derivative of
%     %       posM(t) = pos1(t)*factor(t) + pos2(t)*(1-factor(t))
%     %and the derivative is
%     %       dposM(t) = dpos1(t)*factor(t) + pos1(t)*dfactor(t) + dpos2(t)*(1-factor(t)) + pos2(t)*(1-dfactor(t))  
%     dpos(selM,:) = ...
%       (   bsxfun(@times,   dpos(sel1,:) - dpos(sel2,:), factors(moving))  ...
%         + bsxfun(@rdivide,  pos(sel1,:)  - pos(sel2,:), TD(moving)) ...
%         + dpos(sel2,:) ...
%       );% / 2;
%     %cancel accelerations
%     dvel(selM,:) = 0;
%   end


    triplet_dpos = (dpos(p2,:)+dpos(pM,:))/2 + bsxfun(@rdivide,  pos(p2,:)-pos(pM,:), (Tf-t+1e-3)./(-TD));
%     triplet_dpos = ...
%       (   bsxfun(@times,   dpos(p1,:) - dpos(p2,:), factors)  ...
%         + bsxfun(@rdivide,  pos(p1,:)  - pos(p2,:), TD) ...
%         + dpos(p2,:) ...
%       ) ;%/ 2;
    dpos(shifted,:) = shiftingMatrix*triplet_dpos;
    %cancel accelerations
    %dvel(shifted,:) = 0;
    dvel(shifted,:) = shiftingMatrix*((dvel(p2,:)+dvel(pM,:))/2);
end
