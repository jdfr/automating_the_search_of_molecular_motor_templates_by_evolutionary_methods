%the main code: calculate spring forces, and deal them to the points, plus
%drag forces. Then, divide by the point mass
function acels = calculateNaturalAcels(ss, dynamicVels, pos)
  %first, deal some local variables to cut access time to properties
  springEnds    = ss.springEnds;
  %pos           = ss.pos;
  vel           = dynamicVels;%ss.vel;
  m             = ss.m;
  k             = ss.k;
  r             = ss.r;
  c             = ss.c;
  u             = ss.u;
  springsMatrix = ss.springsMatrix;
  nd            = ss.ndims;
  %spring vectors
  springVectors = pos(springEnds(:,2),:)-pos(springEnds(:,1),:);
  %lengths of spring vectors
  springLengths = realsqrt(sum(springVectors.*springVectors, 2));
  %springs' displacements from rest length, signed
     %springDisplacements = r-springLengths;
  %hookean forces exerted by springs. They are for springs' start
  %points. Springs' end points are the same with switched sign
  hookeanForces = -k.*(r-springLengths);
  %spring unit vectors
  %  this may be needed in some cases
     springUnitVectors = zeros(size(springVectors));
     notzero = springLengths>1e-1;
     switch nd
       case 2
         springUnitVectors(notzero,1) = springVectors(notzero,1)./springLengths(notzero);
         springUnitVectors(notzero,2) = springVectors(notzero,2)./springLengths(notzero);
       otherwise
         springUnitVectors(notzero,:) = bsxfun(@rdivide, springVectors(notzero,:), springLengths(notzero));
     end
%   %  unsafe but faster alternative
%      %springUnitVectors = bsxfun(@rdivide, springVectors, springLengths);
%   %relative speeds of springs' ends
%   springSpeeds = vel(springEnds(:,2),:)-vel(springEnds(:,1),:);
%   %relative speeds of springs' ends along springVector direction:
%   %projection of springSpeeds into springVectors: dot product with unit
%   %vector
%   springDirectSpeeds = sum(springSpeeds.*springUnitVectors,2);
%   %damping forces for each spring. Note that they aren't negated
%   %as it may seem sound, because of the chosen reference frame
%   dampForces = c .* springDirectSpeeds;
  %vector force for each spring: hookean + damp forces
  switch nd
    case 2
      allF = hookeanForces;%+dampForces;
      springVectorForces = springUnitVectors;
      springVectorForces(:,1) = springVectorForces(:,1).*allF;
      springVectorForces(:,2) = springVectorForces(:,2).*allF;
    otherwise
      springVectorForces = bsxfun(@times, springUnitVectors, hookeanForces);%+dampForces);
  end
  %compute hookean forces for each dynamic point:
  %spring forces for first ends minus spring forces for second ends
  switch nd
    case 2
      acels = springsMatrix*springVectorForces - dynamicVels*u;
      md = m(ss.dynamical_p);
      acels(:,1) = acels(:,1)./md;
      acels(:,2) = acels(:,2)./md;
    otherwise
      acels = bsxfun(@rdivide, springsMatrix*springVectorForces - dynamicVels*u, m(ss.dynamical_p));
  end  
%   sps    = [2 5 2]';
% %   ps1    = springEnds(sps,1);
% %   ps2    = springEnds(sps,2);
% %   spvcs  = springVectors(sps,:);
% %   spls   = springLengths(sps,:);
%   perps1 = [-springUnitVectors(sps,2), springUnitVectors(sps,1)];%bsxfun(@rdivide, [-spvcs(:,2), spvcs(:,1)], spls);
%   torque = 10*cos(1*ss.t)*[1 0 0]';
%   fors1  = bsxfun(@times, perps1, torque);
%   %acels(1:end) = 0;
%   acels  = acels + bsxfun(@rdivide, springsMatrix(:,sps)*fors1, m(ss.dynamical_p));
  
%%%%ANISOTROPIC DAMP%%%%%
  %springPosCentroids = (pos(springEnds(:,2),:)+pos(springEnds(:,1),:))/2;
  springVelCentroids = (vel(springEnds(:,2),:)+vel(springEnds(:,1),:))/2;
  springUnitNormals  = [-springUnitVectors(:,2), springUnitVectors(:,1)];
  anisodamp          = sum(springUnitNormals.*springVelCentroids,2).*springLengths;
%   fprintf('Min: %s, Max: %s\n', mat2str(min(anisodamp)), mat2str(max(anisodamp)));
  dampSprings        = -springUnitNormals;
  dampSprings(:,1)   = dampSprings(:,1).*anisodamp;
  dampSprings(:,2)   = dampSprings(:,2).*anisodamp;
  dampacels          = abs(springsMatrix)*dampSprings;
  dampacels(:,1)     = dampacels(:,1)./md;
  dampacels(:,2)     = dampacels(:,2)./md;
  acels              = acels + dampacels;
  
%   acels = acels+(rand(size(acels))-0.5)*0.1;
  
%   rad = ss.rad; rdist = ss.stick.dist;
% toTest      = ss.delaunay(:,[1 2])<=size(ss.pos, 1);
% toTest      = toTest(:,1) & toTest(:,2);
% toTest(1)   = false;
% segs        = ss.delaunay(toTest, [1 2]);
% % segs        = segs(randperm(size(segs,1)),:);
% relative    = ss.pos(segs(:,2),:)-ss.pos(segs(:,1),:);
% sqrdD       = realsqrt(sum(relative.*relative, 2));
% 
% minD        = rad(segs(:,1)) + rad(segs(:,2));
% stickD      = minD + rdist(segs(:,1)) + rdist(segs(:,2));
% 
% toConsider    = (sqrdD<=stickD);
%   units      = relative(toConsider,:);
%   sqrdD      = sqrdD(toConsider);
%   segs       = segs(toConsider,:);
%   zz=15;
%   units(:,1) = zz*units(:,1)./sqrdD;
%   units(:,2) = zz*units(:,2)./sqrdD;
%   for k=size(units, 1):-1:1%randperm(size(units, 1))%1:size(units, 1)
%     acels(segs(k,1),:) = acels(segs(k,1),:)+units(k,:);
%     acels(segs(k,2),:) = acels(segs(k,2),:)-units(k,:);
% %     ss.vel(segs(k,1),:) = ss.vel(segs(k,1),:)+units(k,:);
% %     ss.vel(segs(k,2),:) = ss.vel(segs(k,2),:)-units(k,:);
%   end


  rad = ss.rad; rdist = ss.stick.dist;
toTest      = ss.delaunay(:,[1 2])<=size(ss.pos, 1);
toTest      = toTest(:,1) & toTest(:,2);
toTest(1)   = false;
segs        = ss.delaunay(toTest, [1 2]);
% segs        = segs(randperm(size(segs,1)),:);
relative    = ss.pos(segs(:,2),:)-ss.pos(segs(:,1),:);
sqrdD       = realsqrt(sum(relative.*relative, 2));

minD        = rad(segs(:,1)) + rad(segs(:,2));
stickD      = minD + rdist(segs(:,1)) + rdist(segs(:,2));

toConsider  = (sqrdD<=stickD);
units       = relative(toConsider,:);
sqrdD       = sqrdD(toConsider);
segs        = segs(toConsider,:);
dists       = sqrdD-minD(toConsider);
forces      = ss.stick.k*dists./sqrdD;
units(:,1) = units(:,1).*forces;
units(:,2) = units(:,2).*forces;
  for k=size(units, 1):-1:1%randperm(size(units, 1))%1:size(units, 1)
    acels(segs(k,1),:) = acels(segs(k,1),:)+units(k,:);%./ss.m(segs(k,1));
    acels(segs(k,2),:) = acels(segs(k,2),:)-units(k,:);%./ss.m(segs(k,2));
%     ss.vel(segs(k,1),:) = ss.vel(segs(k,1),:)+units(k,:);
%     ss.vel(segs(k,2),:) = ss.vel(segs(k,2),:)-units(k,:);
  end

  
%   %springPosCentroids = (pos(springEnds(:,2),:)+pos(springEnds(:,1),:))/2;
%   springVelCentroids = (vel(springEnds(:,2),:)+vel(springEnds(:,1),:))/2;
%   springUnitNormals  = [-springUnitVectors(:,2), springUnitVectors(:,1)];
%   dampSprings        = - bsxfun(@times, springUnitNormals, sum(springUnitNormals.*springVelCentroids,2).*springLengths);
%   acels  = acels + bsxfun(@rdivide, abs(springsMatrix)*dampSprings, m(ss.dynamical_p));
  
%   n = [11; 22];
%   f = 1;
%   f1 = 30;
%   f2 = f1;
%   acels(n,:) = bsxfun(@plus, acels(n,:), [f1*cos(f*ss.t) 0]);%[f1*sin(f*ss.t) f2*cos(f*ss.t)];
  
  %acels(ps2,:) = acels(ps2,:) + bsxfun(@rdivide, fors2, m(ps2));
%   %profilaxis, it pays to locate errors inside this function, before
%   %they spread over the code
%   if any(isnan(forces(:)))
%       error('Some point''s force is NaN');
%   end;
%   if ~isreal(forces)
%       error('Some point''s force is complex');
%   end
end
