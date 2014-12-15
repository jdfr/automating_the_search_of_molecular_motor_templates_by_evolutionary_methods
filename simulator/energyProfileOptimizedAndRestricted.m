%this is an optimized version of energyProfile. It only works if ts is such
%that
%   -ts(1) is a timestep with a fully recorded state (the start of a "dyn"
%    series in rec.dyn
%   -ts(2:end) are "dyn" timesteps, all of them relative to the fully
%    recorded state in ts(1)
%   -all timesteps in ts are consecutive
%   -all points are "dynamical_p"
%   -no point is "dynamical_m"
%   -no spring is "dynamical_k" or "dynamical_r"
%   -number of dimensions is 2
%   -rec is of class MemRecorder
%
% This version is used because the safer one is too slow for very long
% sequences
function [e, ts] = energyProfileOptimizedAndRestricted(rec, ts, verbose)

if ~exist('verbose', 'var')
  verbose = false;
end

e      = zeros(numel(ts),1);
baseSS = [];
index  = 1;
if verbose
  fprintf('num iters %d\n', numel(ts));
end
%only for the first, "notdyn" timestep
playSimulation(rec, ts(1), @energyRecorder);

ss      = baseSS;
svel1   = ss.selector.vel(1,1);
svel2   = ss.selector.vel(end,2);
%now, extract posX and posY
selPX = ss.selector.pos(1,1):ss.selector.pos(1,2);
selPX = reshape(selPX, numel(selPX)/2, 2);
selPY = selPX(:,2);
selPX = selPX(:,1);
%make them relative to spring ends
selPX1 = selPX(ss.springEnds(:,1));
selPX2 = selPX(ss.springEnds(:,2));
selPY1 = selPY(ss.springEnds(:,1));
selPY2 = selPY(ss.springEnds(:,2));
m       = ss.m;
k       = ss.k;
r       = ss.r;
initDyn = ts(2);
initI   = find(rec.dynTs==initDyn);
lastI   = initI+numel(ts)-1-1; %-1 to adjust. another -1 for the first one, not included
if any(rec.dynTs(initI:lastI)~=ts(2:end))
  error('This optimized code needs the timesteps to be consecutive and all of them to be derived of the same "notdyn" state');
end
zzz=svel1:svel2;
cellfun(@derivedEnergyRecorder, rec.dyn(initI:lastI));

e(1) = e(1)*2;
e    = e/2;

  function energyRecorder(t, ss)
    baseSS = ss;
    e(index) = calculateEnergy(ss);
    index = index+1;
  end

  function derivedEnergyRecorder(dyn)
    %add them to find each point's kinetic energy by multiplying them by the
    %point mass, then sum it to find all kinetic energies
    e_kinetic = sum(sum(reshape(realpow(dyn(zzz), 2), [], 2), 2).*m);
    %calculate potential energy for each timestep
    e_potential = sum(k.*realpow(r-realsqrt(realpow(dyn(selPX2)-dyn(selPX1), 2) + realpow(dyn(selPY2)-dyn(selPY1), 2)),2));

    e(index) = e_kinetic+e_potential;
    index = index+1;
    if verbose && (mod(index,1000)==0)
      fprintf('iters %d\n', index);
    end
  end
end


% e      = zeros(numel(ts),1);
% baseSS = [];
% index  = 1;
% 
% if verbose
%   fprintf('num iters %d\n', numel(ts));
% end
% %only for the first, "notdyn" timestep
% playSimulation(rec, ts(1), @energyRecorder);
% 
% ss      = baseSS;
% npoints = size(ss.pos, 1);
% initDyn = ts(2);
% initI   = find(rec.dynTs==initDyn);
% lastI   = initI+numel(ts)-1-1; %-1 to adjust. another -1 for the first one, not included
% if any(rec.dynTs(initI:lastI)~=ts(2:end))
%   error('This optimized code needs the timesteps to be consecutive and all of them to be derived of the same "notdyn" state');
% end
% cellfun
% c       = rec.dyn(initI:lastI); %we want it to be a row vector
% %make a giant matrix!
% states  = cat(2,c{1,:});
% 
% %%%extract kinetic energy%%%
% %power the speed components
% vp = realpow(states(ss.selector.vel(1,1):ss.selector.vel(end,2), :), 2);
% %add them to find each point's kinetic energy by multiplying them by the
% %point mass, then sum it to find all kinetic energies
% e_kinetic = sum(bsxfun(@times, vp(1:npoints,:)+vp((npoints+1):end,:), ss.m));
% %now, extract posX and posY
% selPX = ss.selector.pos(1,1):ss.selector.pos(1,2);
% selPX = reshape(selPX, numel(selPX)/2, 2);
% selPY = selPX(:,2);
% selPX = selPX(:,1);
% %make them relative to spring ends
% selPX1 = selPX(ss.springEnds(:,1));
% selPX2 = selPX(ss.springEnds(:,2));
% selPY1 = selPY(ss.springEnds(:,1));
% selPY2 = selPY(ss.springEnds(:,2));
% %find spring lengths and displacements
% springLengths = realsqrt(realpow(states(selPX2,:)-states(selPX1,:), 2) + realpow(states(selPY2,:)-states(selPY1,:), 2));
% springDisplacements = bsxfun(@minus, ss.r, springLengths);
% %calculate potential energy for each timestep
% e_potential = sum(bsxfun(@times, ss.k, realpow(springDisplacements,2)));
% 
% e_total = (e_kinetic+e_potential)/2;
% 
% e(2:end) = e_total;
% 
% %   springVectors = ss.pos(ss.springEnds(:,2),:)-ss.pos(ss.springEnds(:,1),:);
% %   springLengths = realsqrt(sum(realpow(springVectors,2), 2));
% %   springDisplacements = ss.r-springLengths;
% %   e_pot = sum( ss.k.*realpow(springDisplacements,2) )/2;
% %   e_vel = sum(sum(realpow(ss.vel,2),2).*ss.m)/2;
% %   e = e_pot+e_vel;
% 
%   function energyRecorder(t, ss)
%     baseSS = ss;
%     eInit(index) = calculateEnergy(ss);
%     index = index+1;
%     if verbose && (mod(index,1000)==0)
%       fprintf('iters %d\n', index);
%     end
%   end
% end