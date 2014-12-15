function [relaxation deformation] = tryPerturbation(ss, seedRnd, pertParams, points)

sowSeed(seedRnd);

ss = prepareSimulationStructure(ss);

if ~exist('points', 'var')
  points   = get3RefPoints(ss);
end

dps      = distanceMatrix(ss.pos(points,:));

idxps    = sub2ind(size(ss.pos), reshape(repmat(points(:)', 3, 1), 1, []), [1 2 3 1 2 3 1 2 3]);

forceSum = pertParams.forceSum;
pertT    = pertParams.pertT;
relaxT   = pertParams.relaxT;
odefun   = pertParams.odefun;
odeOps   = odeset;
odeOps.Vectorized  = 'on';
odeOps.JPattern    = getJPattern(ss.pos, ss.springEnds);
odeOps.InitialStep = pertT/1000;
odeOps.AbsTol      = 1e-3;
%odeOps.Stats       = 'on';

nsteps   = 1000;
tspan    = linspace(0, pertT, nsteps);

forceDirections = rand(size(ss.pos))-0.5;
forceDirections = bsxfun(@rdivide, forceDirections, realsqrt(sum(realpow(forceDirections, 2), 2)));
forceMags       = rand(size(ss.pos,1), 1);
forceMags       = forceMags*(forceSum/sum(forceMags));
forces          = bsxfun(@times, forceDirections, forceMags);

useRec = false;

if useRec
  rec = makeNewRecorder(ss, 'MemRecorder', 10000, 10000, 100, 100);
  rec = recordAllState(rec, ss);
end

[T Y] = odefun(@innerLoop1, tspan, ss.pos(:), odeOps, forces(:), ss.r, ss.springEnds, ss.springsMatrix, ss.npoints);

last = Y(end,:)';

data = Y(:,idxps)';
d12  = realsqrt(sum(realpow([data(1,:)-data(4,:); data(2,:)-data(5,:); data(3,:)-data(6,:)], 2), 1))/dps(1,2)-1; 
d13  = realsqrt(sum(realpow([data(1,:)-data(7,:); data(2,:)-data(8,:); data(3,:)-data(9,:)], 2), 1))/dps(1,3)-1; 
d23  = realsqrt(sum(realpow([data(4,:)-data(7,:); data(5,:)-data(8,:); data(6,:)-data(9,:)], 2), 1))/dps(2,3)-1;

deformation = struct('t', T, 'd12', d12, 'd13', d13, 'd23', d23);

if useRec
  for k=1:size(Y,1);
    rec = recordDynState(rec, ss, T(k), struct('pos', reshape(Y(k,:), [], 3)));
  end
end

clear T Y;

nsteps   = 1000;
tspan    = linspace(0, relaxT, nsteps);

[T Y] = odefun(@innerLoop2, tspan, last, odeOps, ss.r, ss.springEnds, ss.springsMatrix, ss.npoints);

data = Y(:,idxps)';
d12  = realsqrt(sum(realpow([data(1,:)-data(4,:); data(2,:)-data(5,:); data(3,:)-data(6,:)], 2), 1))/dps(1,2)-1; 
d13  = realsqrt(sum(realpow([data(1,:)-data(7,:); data(2,:)-data(8,:); data(3,:)-data(9,:)], 2), 1))/dps(1,3)-1; 
d23  = realsqrt(sum(realpow([data(4,:)-data(7,:); data(5,:)-data(8,:); data(6,:)-data(9,:)], 2), 1))/dps(2,3)-1;

relaxation = struct('t', T, 'd12', d12, 'd13', d13, 'd23', d23);

if useRec
  T=T+max(deformation.t);
  for k=1:size(Y,1);
    rec = recordDynState(rec, ss, T(k), struct('pos', reshape(Y(k,:), [], 3)));
  end

  clear T Y;
  plot3(deformation.d12, deformation.d13, deformation.d23, 'b', relaxation.d12, relaxation.d13, relaxation.d23, 'r');
  grid on;
  axis equal;
end
function jpat = getJPattern(pos, springEnds)
%this is only valid if all points are dynamic!!!!
sps           = springEnds;
np            = size(pos,1);
%a==accelerations depends on positions of springs of points, plus self-positions and self-velocities, plus position and velocity of ball  
a             = double(logical(sparse([sps(:,1); sps(:,2); (1:np)'], [sps(:,2); sps(:,1); (1:np)'], 1, np, np)));
% %b==velocities and accelerations depend on self-positions and self-velocities 
% b             = speye(np);
jpat          = [a a a; a a a; a a a];

function output = innerLoop2(t, input, r, springEnds, springsMatrix, np)
  ns = size(input,2);
  pos = reshape(input, np, 3, []);%  [reshape(input(1:np,:), [], 1), reshape(input((np+1):(2*np),:), [], 1), reshape(input((2*np+1):end,:), [], 1)];
%   mind=(1:size(springEnds,1))';
%   springEnds = springEnds(mind(:,ones(1,size(input,1))),[1;2]); %springEnds = repmat(springEnds, size(input,2), 1);
  springVectors = pos(springEnds(:,2),:,:)-pos(springEnds(:,1),:,:);
  %lengths of spring vectors
  springLengths = springVectors.*springVectors;
  springLengths = realsqrt(springLengths(:,1,:)+springLengths(:,2,:)+springLengths(:,3,:));
  %springs' displacements from rest length, signed
     %springDisplacements = r-springLengths;
  %hookean forces exerted by springs. They are for springs' start
  %points. Springs' end points are the same with switched sign
  springUnitVectors        = zeros(size(springVectors));
  springUnitVectors(:,1,:) = springVectors(:,1,:)./springLengths;
  springUnitVectors(:,2,:) = springVectors(:,2,:)./springLengths;
  springUnitVectors(:,3,:) = springVectors(:,3,:)./springLengths;

  %vector force for each spring: hookean + damp forces
  hookeanForces = -(r(:,1,ones(1,ns))-springLengths);
  springVectorForces = springUnitVectors;
  springVectorForces(:,1,:) = springVectorForces(:,1,:).*hookeanForces;
  springVectorForces(:,2,:) = springVectorForces(:,2,:).*hookeanForces;
  springVectorForces(:,3,:) = springVectorForces(:,3,:).*hookeanForces;
  %compute hookean forces for each dynamic point:
  %spring forces for first ends minus spring forces for second ends
  acels = zeros(size(pos));
  for k=1:ns
    acels(:,:,k) = springsMatrix*springVectorForces(:,:,k);
  end
  output = reshape(acels, size(input));%[reshape(acels(:,1), np, []); reshape(acels(:,2), np, []); reshape(acels(:,3), np, [])];

function output = innerLoop1(t, input, forces, r, springEnds, springsMatrix, np)
output = innerLoop2(t, input, r, springEnds, springsMatrix, np) + forces(:,ones(1, size(input,2)));