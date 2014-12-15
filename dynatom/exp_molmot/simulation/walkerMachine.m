function [ss r rdyn veryChanged] = walkerMachine(ss, r, rdyn, t)

veryChanged = false;

% colPending = true;

nevermind = ss.vel;
[nevermind, nevermind, contactps, newlocked, nballlocked] = collisionsWithRow(ss, nevermind, nevermind, nevermind, nan);%, true);
contacts = false(size(ss.m));
contacts(contactps) = true;
news     = false(size(ss.m));
news(newlocked) = true;
if ~isempty(newlocked)
  ss.row.locked(newlocked) = true;
  ss.row.lockednball(newlocked) = nballlocked;
  veryChanged = true;
end

%d3 = size(ss.pos,2)>2;

for k=1:numel(ss.legs)
  leg = ss.legs(k);
  switch leg.state(1)
    case 's' %sticky
%       if colPending
%         %detect collisions of sticky toes and the row
%         [colPending contacts] = getContacts(colPending, ss);
%       end
      %when attachMode~=0, not all contacts being toes in state sticky
      %are also newly locked, since they might be touching the filament
      %in the wrong place. Thus, we must use just the newly locked points
      %to be sure
      if ss.row.correct.attach
        cond = any(news(leg.toes));
      else
        cond = any(contacts(leg.toes));
      end
      if cond
        %if a collision is detected, change state and trigger
        %conformational change
        bnds                   = leg.ABCbinding;
        ABCw                   = leg.ABCweight;
        ABCb                   = ss.pos(bnds,:);
        ABCv                   = ss.vel(bnds,:);
        atpp                   = sum(bsxfun(@times, ABCb, ABCw), 1)/sum(ABCw);
        atpv                   = sum(bsxfun(@times, ABCv, ABCw), 1)/sum(ABCw);
        [ss ss.legs(k).atp]    = addPoints(ss, atpp, atpv);
        dists                  = realsqrt(sum(realpow(bsxfun(@minus, ABCb, atpp), 2), 2));
        softChg = isfield(ss.row, 'spchg') && ss.row.spchg.doit;
        if softChg
          ss.r                 = r;
          ss.rdyn              = rdyn;
          rs                   = dists;
          if leg.dynLengths
            newr               = dists .* leg.dynLengthFactor;
          else
            newr               = leg.atpLengths;
          end
        else
          if leg.dynLengths
            rs                 = dists .* leg.dynLengthFactor;
          else
            rs                 = leg.atpLengths;
          end
        end
        [ss newsprs]           = addSprings(ss, [repmat(size(ss.pos,1), size(rs)), bnds], leg.atpK, rs);
        if softChg
          ss.rdyn              = perturbateR(ss.rdyn, ss.r, t, newsprs, newr);
          ss.rdyn.rt(newsprs)  = ss.row.spchg.constT;
        end
        r                      = ss.r;
        rdyn                   = ss.rdyn;
        ss.legs(k).atpb        = newsprs;
        ss                     = prepareSimulationStructure(ss);
        ss.legs(k).atpt        = t;
        ss.legs(k).state       = 'atpChange';
        ss.legs(k).atpCount    = leg.atpCount+1;
        veryChanged            = true;
        if ss.record.toes1
          %this is a lightweight recording
          idx                  = ss.rc.it;
          ps1                  = ss.pos(ss.legs(1).toes,:);
          ps2                  = ss.pos(ss.legs(2).toes,:);
          ss.rc.Pt(idx,:)      = [sum(ps1,1)/size(ps1,1), sum(ps2,1)/size(ps2,1)];
          ss.rc.Ft(idx)        = k;
          ss.rc.Tt(idx)        = t;
          ss.rc.it             = ss.rc.it+1;
        end
      end
    case 'a' %ATP conformational change
      if (t-leg.atpt)>leg.timeAtp
        ss.legs(k).state                          = 'unsticky';
        ss.row.attractive(leg.toes)               = false;
        %remove persistent contacts
        ss.row.locked(leg.toes)                   = false;
        veryChanged            = true; %only true if we wanna view it in BallPlotter
      end
    case 'u' %unsticky
%       if colPending
%         [colPending contacts] = getContacts(colPending, ss);
%       end
      if ~any(contacts(leg.toes))
        %no contact with the row, so we can change state and undo the
        %conformational change
        ss.legs(k).state                          = 'relaxation';
        ss.legs(k).atpt                           = t;
        ss.legs(k).relaxCount                     = leg.relaxCount+1;
        k1                                        = mod(k,numel(ss.legs))+1;
        if ~isempty(ss.legs(k).atpb)
          [ss indexesAfterRemove]                 = removeSprings(ss, ss.legs(k).atpb);
          if ~isempty(ss.legs(k1).atpb)
            ss.legs(k1).atpb                      = indexesAfterRemove(ss.legs(k1).atpb);
          end
        end
        r                                         = ss.r;
        rdyn                                      = ss.rdyn;
        if ~isempty(ss.legs(k).atp)
          ss.row.activePoints                     = ss.row.activePoints(ss.row.activePoints~=ss.legs(k).atp); %just in case
          [ss indexesAfterRemove]                 = removePoints(ss, ss.legs(k).atp);
          if ~isempty(ss.legs(k1).atp)
            ss.legs(k1).atp                       = indexesAfterRemove(ss.legs(k1).atp);
          end
        end
        ss.legs(k).atpb                           = [];
        ss.legs(k).atp                            = [];
        ss                                        = prepareSimulationStructure(ss);
        veryChanged                               = true;
      end
    case 'r' %relaxation conformational change
      if (t-leg.atpt)>leg.timeRelax
        ss.legs(k).state                          = 'sticky';
        ss.row.attractive(leg.toes)               = true;
        veryChanged                               = true;
      end
    otherwise
      error('Unknown leg state: %s!!!!', any2str(leg.state));
  end
end

if ss.record.cm
  idx                         = ss.rc.i;
  if idx>size(ss.rc.P,1)
    ss.rc.P(ceil(end*1.1),1) = 0;
    ss.rc.F(ceil(end*1.1),1) = false;
  end
%   ss.rc.P(idx,:)             = sum(ss.pos/size(ss.pos,1));
%   w                          = ss.m/sum(ss.m);
  %ps                         = ss.pos;
  ps                         = ss.pos(:,1);
  ss.rc.P(idx)               = sum(ps)/numel(ps);
  %ss.rc.P(idx,1)             = sum(ps(:,1))/size(ps,1);%, ps(:,2).*w]); %[sum(ss.pos(:,1).*w) sum(ss.pos(:,2).*w)];
  contacts(ss.row.locked)    = true;
%   ts                         = [ss.legs(1).toes(:); ss.legs(2).toes(:)];
%   ss.rc.F(idx)               = any(contacts(ts));
  ss.rc.F(idx)               = any(contacts(ss.legs(1).toes)) || any(contacts(ss.legs(2).toes));
  ss.rc.i                    = ss.rc.i+1;
end

