function [vels, acels, varargout] = collisionsWithRow(ss, vels, acels, r, t, allPoints) %#ok<INUSL>

%executes collisions between an infinite row of inamovible balls and
%mobile points
%
%these inamobible points are in the line y=0, with a ball at x=0, and balls
%uniformly spaced
%
% EXTREMELY IMPORTANT: THIS CODE WILL NOT WORK PROPERLY IF THIS CONDITION
% HOLDS AT ANY MOMENT, SINCE IT IS ASSUMED THAT A POINT ONLY CAN INTERACT
% WITH TWO BALLS (THOSE INSIDE ITS SEGMENT) AT ANY GIVEN MOMENT
%                    any(ballSep<=(ballAllr+ss.stick.allr))
%
% [NOT LONGER TRUE]
% EXTREMELY IMPORTANT: THIS CODE WILL NOT WORK PROPERLY IF THIS CONDITION
% HOLDS AT ANY MOMENT, SINCE IT IS ASSUMED THAT A POINT ONLY CAN INTERACT
% WITH A BALL AT ANY GIVEN MOMENT
%                    any(ballSep<=(2*ballAllr+ss.stick.allr))

  row            = ss.row;
  ballAllr       = row.ballAllr;
  pos            = ss.pos;
  allr           = ss.stick.allr;
  n2             = nargout>2;
  d3             = size(vels,2)>2;
  if nargin<6
    allPoints    = row.allPoints;
  end
  if allPoints
    %for each point, the distance at which the effect becomes active
    candAllr     = allr + ballAllr;
    %all points are affected, although only some of them will be able to
    %get caught by the sticky regions
    proposed     =             pos(:,2)<candAllr;
    if d3
    proposed     = proposed & (pos(:,3)<candAllr);
    end
    candAllr     = candAllr(proposed);
    proposed     = find(proposed);
  else
    %the affected points are those which can be caught by the sticky
    %regions
    activePoints = row.activePoints;
    %for each point, the distance at which the effect becomes active
    candAllr     = allr(activePoints) + ballAllr;
    proposed     =             pos(activePoints,2)<candAllr;
    if d3
    proposed     = proposed & (pos(activePoints,3)<candAllr);
    end
    candAllr     = candAllr(proposed);
    proposed     = activePoints(proposed);
  end
  attractive     = row.attractive(proposed);
  ballPenk       = row.ballPenk;
  rad            = ss.row.rad; %ss.stick.rad;
  br             = row.ballRad;
  nant           = isnan(t);
  %separation between balls
  ballSep        = row.ballSep;
  lockFreeze     = row.lockFreeze;
  if lockFreeze
    freezed      = row.locked;
  end
  if isempty(proposed)
    candidates   = zeros(0,1);
  else
    %fold X axis for candidates into segments of length "ballSep"
    foldedX      = mod(pos(proposed,1), ballSep);
    doblex       = [foldedX; ballSep - foldedX];
    %calculate distances between each point and the balls at the two
    %endpoints of its segment (in the segment, left balls are located at 
    %[0, 0] and right balls at [ballSep, 0])
    if d3
    dobley       = pos(proposed,[2 3]);
    dobley       = dobley.*dobley;
    dobley       = dobley(:,1)+dobley(:,2);
    else
    dobley       = pos(proposed,2);
    dobley       = dobley.*dobley;
    end
    dists        = realsqrt([dobley; dobley] + (doblex.*doblex));
    %filter definitive candidates
    cands        = reshape(dists<[candAllr; candAllr], [], 2);
    c1           = proposed(cands(:,1),1);
    c2           = proposed(cands(:,2),1);
    candidates   = [c1; c2];
  end
  if isempty(candidates)
    if nant
      varargout(1:3) = cell(1,3);
      return
    else
      nballs     = zeros(0,1);
      if n2
        ballRad  = zeros(0,1);
        radCands = zeros(0,1);
      end
    end
  else
    
    %get offset for points contacting balls at the right side of the segment
    nc1          = 1+sum(cands(:,1));%sum(c1);
    distsc       = dists(cands);
    %typeSpec == [A B], where A is a period and B is a sub-period length,
    %specifying the types of balls. For example, if A=5 and B=2, the pattern
    %will be:
    %     ... * * # # # * * # # # * * # # # * * # # # ...
    % balls of type (*) will have ballRad as specified by the user, while
    % balls of type (#) will have ballRad:=ballAllr, making them unable to be
    % sticky
    %get the index for each ball (right balls' index must be shifted by 1,
    %since they are *the same* as the left ball of the next segment)
    typeSpec             = row.typeSpec;
    switch row.attachMode
      %case 0: do nothing
      case 1 %attach only on the Y upside
        attractive(attractive) = (pos(proposed(attractive),2)>0);
      %other cases most complicated: attach site depends on the leg
    end
    nballs               = floor(pos(candidates,1)/ballSep);
    nballs(nc1:end)      = nballs(nc1:end)+1;
    actives              = [attractive(cands(:,1)); attractive(cands(:,2))] & (mod(nballs, typeSpec(1))<typeSpec(2));
    %THIS IS INEFFICIENT: WE ONLY NEED TO ENTER HERE IF:
    %       (~nant) && ( n2 || (~lockFreeze) || (~all(freezed(candidates(:)))) )    
    if (~nant) && ( n2 || (~lockFreeze) || (~all(freezed(candidates(:)))) )  
      ballRad            = zeros(size(nballs))+ballAllr;
      ballRad(actives)   = br;
      radCands           = zeros(size(candidates));
      radCands(actives)  = rad(candidates(actives));
      radCands(~actives) = allr(candidates(~actives));
      %In fact, the force should be
      %          -ballPenk.*(distsc-radCands-ballRad)./distsc
      %but we divide it one more time by distsc to include the normalization
      %factor for the diff vector
      forces    = -ballPenk.*(distsc-radCands-ballRad)./(distsc.*distsc);
      %get relative unit vectors in the direction from the ball to the
      %candidate point. 
      %For left balls, the relative vectors are the positions of the points in
      %the segment. For right balls, the X component must be shifted by
      %"ballSep" to the left before the position can be used 
      if d3
      diffs      = [[foldedX(cands(:,1)); foldedX(cands(:,2))], pos(candidates,[2 3])];
      else
      diffs      = [[foldedX(cands(:,1)); foldedX(cands(:,2))], pos(candidates,2)];
      end
      diffs(nc1:end, 1) = diffs(nc1:end, 1)-ballSep;
      %make the vectors unitary, and multiply by the force in order to obtain
      %the force vector
      diffs(:,1) = diffs(:,1).*(forces);
      diffs(:,2) = diffs(:,2).*(forces);
      if d3
      diffs(:,3) = diffs(:,3).*(forces);
      end
%       %add acels
%       acels(candidates,:) = acels(candidates,:) + diffs;
%       %now, there is a special case: points contacting both left and irght
%       %balls. As per Matlab indexing, only the last index is used in
%       %assignment, so the first one must also be added. First, identify those
%       %points:
%       bothballsI = proposed(bothballs(:,1));
%       acels(bothballsI,:) = acels(bothballsI,:) + diffs(indexboth(:,1),:);
      acels(c1,:) = acels(c1,:) + diffs(1:nc1-1,:);
      acels(c2,:) = acels(c2,:) + diffs(nc1:end,:);
    end
  end  
  
  if nant
  %         HAY QUE TENER EN CUENTA LOS SIGUIENTES CASOS:
  %             -AMBOS REPULSIVOS: NO HAY QUE BLOQUEAR NINGUNO
  %             -AMBOS ATRACTIVOS: HAY QUE BLOQUEAR EL MAS CERCANO
  %             -UNO DE CADA: HAY QUE BLOQUEAR EL ATRACTIVO.
    notlocked   = ~row.locked(candidates);
    %provisional list of newly locked points
    toBeLocked  = (notlocked) & actives;
    bothballs  = cands(:,1) & cands(:,2);
    indexboth  = reshape(cumsum(cands(:)), [], 2);
    indexboth  = indexboth(bothballs, :);
    %select points that are binded to two balls, both of them active
    sels = notlocked(indexboth(:,1)) & actives(indexboth(:,1)) & actives(indexboth(:,2));
    if any(sels)
      indexboth  = indexboth(sels,:);
      firstIsCloser = distsc(indexboth(:,1))<=distsc(indexboth(:,2));
      toBeLocked(indexboth(firstIsCloser,2))  = false;
      toBeLocked(indexboth(~firstIsCloser,1)) = false;
    end

    varargout{1} = candidates;
    varargout{2} = candidates(toBeLocked);
    varargout{3} = nballs(toBeLocked);
  else
    if lockFreeze
      %freezed            = row.locked;
      %fk                 = row.freezeK;
      acels(freezed,:)   = 0;%dirtyAcels(freezed,:)-vels(freezed,:)*ss.u;
      vels(freezed,:)    = 0;%vels(freezed,:)*fk;
    else
      uselocked = row.locked;
      lockednball = row.lockednball;
      uselocked(candidates(lockednball(candidates)==nballs)) = false;
      nballs2 = lockednball(uselocked);
      if ~isempty(nballs2)
        diffs = pos(uselocked,:);
        diffs(:,1) = diffs(:,1)-nballs2*ballSep;
        ddists2 = sum(realpow(diffs, 2), 2);
        %In fact, the force should be
        %  -ballPenk.*(realsqrt(ddists2)-rad(uselocked)-br)./realsqrt(ddists2)
        %but we divide it one more time by realsqrt(ddists2) (thereby dividing
        %it by ddists2) to include the normalization factor for the diff vector
        forces    = -ballPenk.*(realsqrt(ddists2)-rad(uselocked)-br)./ddists2;
        diffs(:,1) = diffs(:,1).*forces;
        diffs(:,2) = diffs(:,2).*forces;
        if d3
        diffs(:,3) = diffs(:,3).*forces;
        end
        acels(uselocked,:) = acels(uselocked,:) + diffs;
        %dirtyAcels(uselocked,:) = dirtyAcels(uselocked,:) + diffs;
      end
    end
    
    if n2
      if lockFreeze
        uselocked = row.locked;
        nballs2 = row.lockednball(uselocked);
      end
      if isempty(nballs2)
        varargout(1:4) = {candidates, ...
                          nballs, ...
                          ballRad, ...
                          radCands};
      else
        varargout(1:4) = {[candidates; find(uselocked)], ...
                          [nballs; nballs2], ...
                          [ballRad; repmat(br, size(nballs2))], ...
                          [radCands; rad(uselocked)]};
      end
    end
  end



