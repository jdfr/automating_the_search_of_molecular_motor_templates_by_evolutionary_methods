%the main code: calculate spring forces, and deal them to the points, plus
%drag forces. Then, divide by the point mass
function [acels ss] = calculateAcels(ss, pos, vel, r, d3) %#ok<INUSL>
%springEnds, springsMatrix, pos, vel, m, k, r, u, au, rad, allr, penk, sap

m   = ss.m;
k   = ss.k;
u   = ss.u;
springEnds = ss.springEnds;
springsMatrix = ss.springsMatrix;

  springVectors = pos(springEnds(:,2),:)-pos(springEnds(:,1),:);
  %lengths of spring vectors
  springLengths = springVectors.*springVectors;
  if d3
    springLengths = realsqrt(springLengths(:,1)+springLengths(:,2)+springLengths(:,3));
  else
    springLengths = realsqrt(springLengths(:,1)+springLengths(:,2));
  end
  %springs' displacements from rest length, signed
     %springDisplacements = r-springLengths;
  %hookean forces exerted by springs. They are for springs' start
  %points. Springs' end points are the same with switched sign
  notzero = springLengths>1e-9;
  springUnitVectors = zeros(size(springVectors));
  spV = springVectors(notzero,:);
  spL = springLengths(notzero);
  springUnitVectors(notzero,1) = spV(:,1)./spL;
  springUnitVectors(notzero,2) = spV(:,2)./spL;
  if d3
  springUnitVectors(notzero,3) = spV(:,3)./spL;
  end

  %vector force for each spring: hookean + damp forces
  hookeanForces = -k.*(r-springLengths);
  springVectorForces = springUnitVectors;
  springVectorForces(:,1) = springVectorForces(:,1).*hookeanForces;
  springVectorForces(:,2) = springVectorForces(:,2).*hookeanForces;
  if d3
  springVectorForces(:,3) = springVectorForces(:,3).*hookeanForces;
  end
  %compute hookean forces for each dynamic point:
  %spring forces for first ends minus spring forces for second ends
  %IMPORTANT: THIS WORKS JUST AS LONG AS EACH SPRING HAS HIS OWN TWO
  %BALLS!!!!!!! OTHERWISE, USE springsMatrix
  acels = springsMatrix*springVectorForces - vel*u; %[vel(:,1).*u, vel(:,2).*u];
  
%%%%ANISOTROPIC DAMP%%%%%
%ATTENTION: THIS ISN'T GONNA FUNCTION PROPERLY IN 3D!!!!!!!!!!!!
if ss.useAnisodamp
  if d3
    error('anisodamping is not prepared for 3D simulations!!!');
  end
  %springPosCentroids = (pos(springEnds(:,2),:)+pos(springEnds(:,1),:))/2;
  springVelCentroids = (vel(springEnds(:,2),:)+vel(springEnds(:,1),:))/2;
  springUnitNormals  = [-springUnitVectors(:,2), springUnitVectors(:,1)];
  anisodamp          = sum(springUnitNormals.*springVelCentroids,2).*springLengths;
%   fprintf('Min: %s, Max: %s\n', mat2str(min(anisodamp)), mat2str(max(anisodamp)));
  dampSprings        = -springUnitNormals;
  dampSprings(:,1)   = dampSprings(:,1).*anisodamp;
  dampSprings(:,2)   = dampSprings(:,2).*anisodamp;
  %IMPORTANT: THIS WORKS JUST AS LONG AS EACH SPRING HAS HIS OWN TWO
  %BALLS!!!!!!! OTHERWISE, USE abs(springsMatrix)
    acels              = acels + abs(springsMatrix)*dampSprings;
  %   acels(springEnds(:,1),:) = acels(springEnds(:,1),:) + dampSprings;
  %   acels(springEnds(:,2),:) = acels(springEnds(:,2),:) + dampSprings;
% % The above is incorrect if a different anisoviscosity "au" is used for each spring
% %   acels(springEnds(:,1),1) = acels(springEnds(:,1),1) + dampSprings(:,1).*au(springEnds(:,1));
% %   acels(springEnds(:,1),2) = acels(springEnds(:,1),2) + dampSprings(:,2).*au(springEnds(:,1));
% %   acels(springEnds(:,2),1) = acels(springEnds(:,2),1) + dampSprings(:,1).*au(springEnds(:,2));
% %   acels(springEnds(:,2),2) = acels(springEnds(:,2),2) + dampSprings(:,2).*au(springEnds(:,2));
%   addAnisoDamp(acels, au, dampSprings, uint32(springEnds));
%%%%ANISOTROPIC DAMP%%%%%
end

if ss.useSAP
  rad  = ss.stick.rad;
  allr = ss.stick.allr;
  penk = ss.stick.penk;
  [prox ss.sap] = doSweepAndPrune(pos, allr, ss.sap);
  relative    = pos(prox(:,2),:)-pos(prox(:,1),:);
  sqrdD       = relative.*relative;
  if d3
  sqrdD       = realsqrt(sqrdD(:,1)+sqrdD(:,2)+sqrdD(:,3));
  else
  sqrdD       = realsqrt(sqrdD(:,1)+sqrdD(:,2));
  end

  minD        = rad(prox(:,1)) + rad(prox(:,2));
  stickD      = allr(prox(:,1)) + allr(prox(:,2));

  toConsider  = (sqrdD<=stickD);
  units       = relative(toConsider,:);
  sqrdD       = sqrdD(toConsider);
  prox        = prox(toConsider,:);
  if ~isempty(prox)
    dists       = sqrdD-minD(toConsider);
    forces      = min(penk(prox(:,1)), penk(prox(:,2))).*dists./sqrdD;
    % forces      = dists./sqrdD;
    % selected    = dists<0;
    % forces(selected,:) = forces(selected,:).*min( ink(prox(selected,1)),  ink(prox(selected,2)));
    % selected    = ~selected;
    % forces(selected,:) = forces(selected,:).*min(outk(prox(selected,1)), outk(prox(selected,2)));
    units(:,1) = units(:,1).*forces;
    units(:,2) = units(:,2).*forces;
    if d3
    units(:,3) = units(:,3).*forces;
    end
  %   for k=1:size(units, 1)
  %     acels(prox(k,1),:) = acels(prox(k,1),:)+units(k,:);
  %     acels(prox(k,2),:) = acels(prox(k,2),:)-units(k,:);
  %   end
    addIndexesToAcels(acels, units, prox);
  end
end
  
  acels(:,1) = acels(:,1)./m;
  acels(:,2) = acels(:,2)./m;
  if d3
  acels(:,3) = acels(:,3)./m;
  end
  
end
