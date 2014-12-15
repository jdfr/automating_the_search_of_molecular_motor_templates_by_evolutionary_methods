%find energy profile
function [e, ts] = energyProfile(rec, limitA, limitB, verbose)

if ismember('TYPE', fieldnames(rec))
  type = rec.TYPE;
else
  type = 'MemRecorder';
end

switch type %rec.TYPE
  case 'MemRecorder'
    if ~exist('verbose', 'var'); verbose = false;             end
    ts      = allTimeSteps(rec);
    if exist('limitB', 'var') && (limitB~=inf)
      ts = ts((ts>=limitA)&(ts<=limitB));
    else
      ts = ts((ts>=limitA));
    end
    [e, ts] = energyProfileOptimizedAndRestricted(rec, ts, verbose);
  case 'EnergyRecorder'
    ts = rec.ts(1:(rec.index-1));
    if exist('limitB', 'var') && (limitB~=inf)
      toGET = (ts>=limitA)&(ts<=limitB);
    else
      toGET = (ts>=limitA);
    end
    ts   = ts(toGET);
    e    = rec.e(toGET);
  otherwise
    error('unsupported Recorder type %s\n', rec.TYPE);
end

% tic;
% e     = zeros(numel(ts),1);
% index = 1;
% if verbose
%   fprintf('num iters %d\n', numel(ts));
% end
% playSimulation(rec, ts, @energyRecorder);
% toc
% 
% tic;
% [e2, ts2] = energyProfileOptimizedAndRestricted(rec, ts, verbose);
% toc
% 
% if any(e~=e2) || any(ts~=ts2)
%   error('joder');
% else
%   fprintf('ok\n');
% end
% return
% 
% 
%   function energyRecorder(t, ss)
%     e(index) = calculateEnergy(ss);
%     index = index+1;
%     if verbose && (mod(index,1000)==0)
%       fprintf('iters %d\n', index);
%     end
%   end
end