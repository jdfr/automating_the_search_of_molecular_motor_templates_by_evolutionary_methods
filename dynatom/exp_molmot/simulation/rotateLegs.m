%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ss toes contactPoints] = rotateLegs(ss, evparams, toes)

if size(ss.pos,2)<3

  %center of mass of toes of leg 1
  t1 = toes{1};
  w1 = ss.m(t1)/sum(ss.m(t1));
  cm1 = sum(bsxfun(@times, ss.pos(t1,:), w1));

  %center of mass of toes of leg 2
  t2 = toes{2};
  w2 = ss.m(t2)/sum(ss.m(t2));
  cm2 = sum(bsxfun(@times, ss.pos(t2,:), w2));

  ss.pos = bsxfun(@minus, ss.pos, cm1);   %displace cm1 to the origin
  t12 = cm2-cm1;                          %get the vector from cm1 to cm2
  t12d = realsqrt(sum(realpow(t12, 2)));  %get distance
  t12u = t12/t12d;                        %get unit vector
  
  %rotate the structure so the line from cm1 to cm2 is horizontal
  %(t12u=[cos(alpha), sin(alpha)], and we want to rotate the structure by
  %the amount -alpha)
  ss.pos = ss.pos*[t12u(1), -t12u(2); t12u(2), t12u(1)];

  %if most of the structure is below, flip it
  if sum(ss.pos(:,2))<0
    ss.pos(:,2) = -ss.pos(:,2);
  end

  %get the lowest toe in each leg
  [nevermind mt1] = min(ss.pos(t1,2));
  mt1 = t1(mt1);
  [nevermind mt2] = min(ss.pos(t2,2));
  mt2 = t2(mt2);

  if isfield(evparams.walker.correct, 'rotation') && evparams.walker.correct.rotation
    ss.pos = bsxfun(@minus, ss.pos, ss.pos(mt1,:));   %displace mt1 to the origin
  else
    ss.pos = bsxfun(@minus, ss.pos, mt1);   %displace mt1 to the origin (ERRONEOUS)
  end
  t12 = ss.pos(mt2,:)-ss.pos(mt1,:);      %get the vector from mt1 to mt2
  t12d = realsqrt(sum(realpow(t12, 2)));  %get distance
  t12u = t12/t12d;                        %get unit vector

  %rotate the structure so the line from mt1 to mt2 is horizontal
  %(t12u=[cos(alpha), sin(alpha)], and we want to rotate the structure by
  %the amount -alpha)
  ss.pos = ss.pos*[t12u(1), -t12u(2); t12u(2), t12u(1)];

  %now, offset the structure so the lowest toes are just in contact with the
  %row
  ss.pos(:,2) = ss.pos(:,2)-ss.pos(mt1,2)+evparams.walker.row.ballAllr+evparams.walker.pointAllr-evparams.walker.initialDepth;

  contactPoints = [mt1; mt2];
else

  switch evparams.walker.d3.rotationMode
    case 'inchworm'
      %this should leave the walking toes in an inchworm position. It is
      %essentially the 3D version of the 2D transformations
      
      %center of mass of toes of leg 1
      t1 = toes{1};
      w1 = ss.m(t1)/sum(ss.m(t1));
      cm1 = sum(bsxfun(@times, ss.pos(t1,:), w1));

      %center of mass of toes of leg 2
      t2 = toes{2};
      w2 = ss.m(t2)/sum(ss.m(t2));
      cm2 = sum(bsxfun(@times, ss.pos(t2,:), w2));

      ss.pos = bsxfun(@minus, ss.pos, cm1);   %displace cm1 to the origin
      t12 = cm2-cm1;                          %get the vector from cm1 to cm2
      t12d = realsqrt(sum(realpow(t12, 2)));  %get distance
      t12u = t12/t12d;                        %get unit vector
      
      %rotate the structure so the line from cm1 to cm2 coincides with the
      %X axis
      rot = vrrotvec2mat(vrrotvec(t12u', [1; 0; 0]))';
      ss.pos = ss.pos*rot;
      
      %rotate the structure over the X axis so the center of mass is on the
      %YX plane, in a situation analogue to being "up" in 2D. To get this
      %done, we must get the vector of the center of mass of the structure,
      %and then project it on the YZ plane. The resulting vector must be
      %rotated to be parallel to the Y vector
      w    = ss.m/sum(ss.m);
      cm   = sum(bsxfun(@times, ss.pos, w));
      cmYZ = [0; cm(2); cm(3)];
      rot  = vrrotvec2mat(vrrotvec(cmYZ, [0; 1; 0]))';
      ss.pos = ss.pos*rot;
      
      %now, the structure is oriented: -Y means "down", +Y means "up".
      %Let's get the lowest toe in each leg
      [nevermind mt1] = min(ss.pos(t1,2));
      mt1 = t1(mt1);
      [nevermind mt2] = min(ss.pos(t2,2));
      mt2 = t2(mt2);
      
      ss.pos = bsxfun(@minus, ss.pos, ss.pos(mt1,:));   %displace mt1 to the origin

      t12 = ss.pos(mt2,:)-ss.pos(mt1,:);      %get the vector from mt1 to mt2
      t12d = realsqrt(sum(realpow(t12, 2)));  %get distance
      t12u = t12/t12d;                        %get unit vector

      %rotate the structure so the line from mt1 to mt2 coincides with the
      %X axis
      rot = vrrotvec2mat(vrrotvec(t12u', [1; 0; 0]))';
      ss.pos = ss.pos*rot;
      
      %now, offset the structure so the lowest toes are just in contact 
      %with the X axis
      ss.pos(:,2) = ss.pos(:,2)-ss.pos(mt1,2)+evparams.walker.row.ballAllr+evparams.walker.pointAllr-evparams.walker.initialDepth;
      
      contactPoints = [mt1; mt2];
      
    case 'handoverhand'
      %this should leave the walking toes hand over hand
      
      %This scheme is far less canalizing than the inchworm. We must do
      %some choices for the orientation of the structure, and, since we
      %cannot devise a good heuristic, we do them mostly on an arbitrary
      %way
      
      if isa(evparams.walker.d3.hoh.pivot,       'function_handle'); pivot  = evparams.walker.d3.hoh.pivot(ss);       else pivot  = evparams.walker.d3.hoh.pivot;       end %select pivoting    point (default should be the last   one)
      if isa(evparams.walker.d3.hoh.pinpoint,    'function_handle'); pinpt  = evparams.walker.d3.hoh.pinpoint(ss);    else pinpt  = evparams.walker.d3.hoh.pinpoint;    end %select pinpoint    point (default should be the first  one)
      if isa(evparams.walker.d3.hoh.orientPoint, 'function_handle'); orient = evparams.walker.d3.hoh.orientPoint(ss); else orient = evparams.walker.d3.hoh.orientPoint; end %select orientating point (default should be the second one)
      
      doMirroring = (~isfield(evparams.walker.d3.hoh, 'duplicateBy')) || strcmp(evparams.walker.d3.hoh.duplicateBy, 'mirroring');
      
      if all(toes{1}~=pinpt) && all(toes{2}~=pinpt)
        contactPoints = 'the shape is not suited for this disposition';
        return
      end
      
      ss.pos = bsxfun(@minus, ss.pos, ss.pos(pivot,:));   %displace pinpoint to the origin
      t12 = ss.pos(pinpt,:)-ss.pos(pivot,:);              %get the vector from pinpoint to pivot
      t12d = realsqrt(sum(realpow(t12, 2)));              %get distance
      t12u = t12/t12d;                                     %get unit vector
      
      %rotate the structure so the line from pivot to pinpt coincides 
      %with the Y axis
      rot = vrrotvec2mat(vrrotvec(t12u', [0; 1; 0]))';
      ss.pos = ss.pos*rot;
      
      %now, in the inchworm approach, we used the center of mass for
      %orientating the structure. That approach was sound because good
      %inchwormers were expected to have a center of mass far from the line
      %between the centers of the legs. Now, this is not the case, so we
      %will use an arbitrary point

      %rotate the structure over the Y axis so the "orient point" is 
      %placed in a precise way. To get this done, we must consider the 
      %arbitrarily chosen "orient point", and then project it on 
      %the XZ plane. The resulting vector must be rotated to be parallel to
      %the X vector.
      opos   = ss.pos(orient,:);
      oXZ    = [opos(1); 0; opos(3)];
      rot    = vrrotvec2mat(vrrotvec(oXZ, [1; 0; 0]))';
      ss.pos = ss.pos*rot;
      
      %record the orientation
      ss.hoh.orientedPos = ss.pos;
      
      pos1 = ss.pos;
      pos2 = ss.pos;
      
      if ~doMirroring
        %the second structure is obtained by rotating 180º around 
        %the Y axis rather than mirroring, as mirroring can easily 
        %yield false dimers
        rot = [-1 0 0; 0 1 0; 0 0 -1];
        pos2 = pos2*rot;
      end
      
      %now that the structure is unequivocally oriented, we will swing it
      %in the XY plane (over the Z axis) by a prespecified amount
      
      swing = evparams.walker.d3.hoh.swingAngle;
      rot1  = [cos(swing), -sin(swing), 0;  sin(swing),  cos(swing), 0; 0 0 1];
      rot2  = [cos(swing),  sin(swing), 0; -sin(swing),  cos(swing), 0; 0 0 1];
      
      pos1 = pos1*rot1';
      pos2 = pos2*rot2';
      
      %now, we will gently swing it in the YZ plane (over the X axis), so
      %that it will become within contact range of the microfilament
      
      dist = evparams.walker.row.ballAllr+evparams.walker.pointAllr-evparams.walker.initialDepth;
      distpivot = t12d * cos(swing); %sin(swing);
      
      %cos_swing_half = realsqrt(1-dist*dist/4/distpivot/distpivot);%dist/2/distpivot;
      %cos_swing = realsqrt((1+cos_swing_half)/2);%2*cos_swing_half*cos_swing_half-1;
      %sin_swing = realsqrt(1-cos_swing.*cos_swing);
      sin_swing = dist/distpivot;
      cos_swing = realsqrt(1-sin_swing*sin_swing);
      
      if doMirroring
        rot  = [1, 0, 0; 0, cos_swing, -sin_swing; 0,  sin_swing, cos_swing];
        pos1 = pos1*rot;
        pos2 = pos2*rot;

        pos2(:,3) = -pos2(:,3);
      else
        rot1 = [1, 0, 0; 0, cos_swing, -sin_swing; 0,  sin_swing, cos_swing];
        rot2 = [1, 0, 0; 0, cos_swing,  sin_swing; 0, -sin_swing, cos_swing];
        pos1 = pos1*rot1;
        pos2 = pos2*rot2;
      end
      
      %displace the structure to its intended final position
      
      %distdesp = t12d * cos(swing); % sin(swing);
      
      pos1(:,2) = pos1(:,2) - pos1(pinpt,2);
      pos2(:,2) = pos2(:,2) - pos2(pinpt,2);
      
      %unite the two legs in one superstructure
      
      switch evparams.walker.d3.hoh.fusionMode
        case 'justOnePoint'
          np                         = size(ss.pos,1);
          spEnds                     = ss.springEnds;
          newSps                     = spEnds+np;
          newSps(newSps==(pivot+np)) = pivot;
          [ss is]                    = copyPoints(ss, (1:np)');
          ss.pos(1:np,:)             = pos1;
          ss.pos(is,:)               = pos2;
          distp                      = realsqrt(sum(realpow(ss.pos(pivot,:)-ss.pos(pivot+np,:), 2)));
          if distp>1e-10
            error('This should never happen!!!!');
          end
          [ss is]                    = copySprings(ss, (1:size(spEnds,1))');
          ss.springEnds(is,:)        = newSps;
          for k=1:numel(toes)
            if (numel(toes{k})==1) && (toes{k}==pivot)
              contactPoints          = sprintf('There is a set of toes composed of just the point %s, but it is also the pivot (%s), so the leg structure makes no sense!!!', any2str(toes{k}), any2str(pivot));
              return
            else
              toes{k}                = toes{k}(toes{k}~=pivot);
            end
          end
          [ss ss.hoh.newIndexes]     = removePoints(ss, pivot+np, (np+1):(2*np));
          %[ss toes{:}]               = removePoints(ss, pivot+np, toes{:});
          
          contactPoints              = [pinpt ss.hoh.newIndexes(pinpt)];
          
        otherwise
          error('option not supported: evparams.walker.d3.hoh.fusionMode=%s', evparams.walker.d3.hoh.fusionMode);
      end
      
      %to visualize the result:
      %zss = ss; drawSS(BallPlotter('axisWindow', getAxisWindow(zss.pos, 0.1), 'axisEqual', true, 'axisSquare', false, 'circleFaces', 10), zss); for k='xyz'; eval(sprintf('%slabel(''%s'');', k, k)); end; line([0 0 nan 0 0 nan 0 100], [0 0 nan 0 100 nan 0 0], [0 100 nan 0 0 nan 0 0]);
      
      ss.hoh.pivot                   = pivot;
      ss.hoh.idxoffset               = np;
      
      %error('handoverhand mode still not implemented!!!!');
      
    otherwise
      
      error('evparams.walker.d3.rotationMode == %s not understood!!!', any2str(evparams.walker.d3.rotationMode));
      
  end

end