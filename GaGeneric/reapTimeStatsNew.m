function info = reapTimeStatsNew(runstr,show)
%reap the time statistics for a given simulation. Runstr is the directory
%of the simulation

ranges = 1;
gens = 1;

if ~exist('show', 'var')
  show = true;
end

stfields = {'minl', 'ming', 'maxl', 'maxg', 'meanl', 'meang', 'medianl', 'mediang', 'numStrings', ...
            'timeSpawned', 'cputimeSpawned', 'timeReaped', 'cputimeReaped', 'totalClusterSpawned', 'totalClusterReaped', ...
            'bestFitness', ...
            };
statdim  = cell(gens, ranges);
args = reshape([stfields; repmat({statdim}, size(stfields))], 1, []);
info = struct(args{:});
clear statdim args;

nomdir = [runstr filesep];

% fpob = fopen([nomdir,'poblacion.txt'],'r');
fres = fopen([nomdir,'resultado.txt'],'r');
ftst = fopen([nomdir,'timestats.txt'],'r');

ntok = '([0-9\.eE+]+)';

AUTOMATED_SPAWNED = {...
  {'numStrings',          'STRINGS: ', true}, ...
  {'timeSpawned',         'ELAPSED=',  false}, ...
  {'cputimeSpawned',      'CPU: ',     false}, ...
  {'totalClusterSpawned', 'CLUSTER=',  true}, ...
  {'minl',                'MIN=',      false}, ...
  {'maxl',                'MAX=',      false}, ...
  {'meanl',               'MEAN=',     false}, ...
  {'medianl',             'MEDIAN=',   false}, ... %true}, ...
...%  {'ming',                'MIN=',      false}, ...
...%  {'maxg',                'MAX=',      false}, ...
...%  {'meang',               'MEAN=',     false}, ...
...%  {'mediang',             'MEDIAN=',   false}, ...
  };

maxgen = -inf;

while ~feof(ftst)
  line = fgetl(ftst);
  if ~isempty(strmatch('SPAWNED', line))
    tok   = regexp(line, ['GEN=' ntok], 'tokens', 'once');
    gen   = str2double(tok{1});
    tok   = regexp(line, ['RANGE=' ntok], 'tokens', 'once');
    if ~isempty(tok)
      range = str2double(tok{1});
    else
      range = 1;
    end
%     if (gen>gens) || (range>ranges)
%       for k=1:numel(AUTOMATED_SPAWNED)
%         if AUTOMATED_SPAWNED{k}{3}
%           fgetl(ftst);
%         end
%       end
%       continue
%     end
    for k=1:numel(AUTOMATED_SPAWNED)
      FIELD   = AUTOMATED_SPAWNED{k}{1};
      LABEL   = AUTOMATED_SPAWNED{k}{2};
      READNEW = AUTOMATED_SPAWNED{k}{3};
      tok   = regexp(line, [LABEL ntok], 'tokens', 'once');
      if ~isempty(tok)
        val   = str2double(tok{1});
        info(gen,range).(FIELD) = val;
      end
      if READNEW
        line = fgetl(ftst);
      end
    end
  elseif (~isempty(strmatch('FINISHED', line))) || (~isempty(strmatch('REAPED', line))) || (~isempty(strmatch('AT LAST, REAPED', line)))
    if (~isempty(strmatch('FINISHED', line)))
      tok   = regexp(line, ['GENERATION ' ntok], 'tokens');
      ges   = cellfun(@(x)str2double(x{1}), tok);
      ras   = 1;
    else
      tok   = regexp(line, ['G=' ntok], 'tokens');
      ges   = cellfun(@(x)str2double(x{1}), tok);
      tok   = regexp(line, ['R=' ntok], 'tokens');
      ras   = cellfun(@(x)str2double(x{1}), tok);
    end
    line  = fgetl(ftst);
    tok   = regexp(line, ['ELAPSED: ' ntok], 'tokens', 'once');
    elaps = str2double(tok{1});
    tok   = regexp(line, ['CPU: ' ntok], 'tokens', 'once');
    if ~isempty(tok)
      cput  = str2double(tok{1});
    else
      cput = nan;
    end
    tok   = regexp(line, ['CLUSTER: ' ntok], 'tokens', 'once');
    clust = str2double(tok{1});
    maxgen = max(maxgen, max(ges));
    for z=1:numel(ges)
%       if (ges(z)>gens) || (ras(z)>ranges)
%         continue;
%       end
      info(ges(z),ras(z)).timeReaped         = elaps;
      info(ges(z),ras(z)).cputimeReaped      = cput;
      info(ges(z),ras(z)).totalClusterReaped = clust;
    end
%   elseif ~isempty(strmatch('   POPULATION', line))
%     %do nothing
%   elseif ~isempty(strmatch('--BYTES', line))
%     %do nothing
%   elseif (~isempty(strmatch('JOB ', line)))           || ...
%          (~isempty(strmatch('  Message:', line)))     || ...
%          (~isempty(strmatch('  Identifier:', line)))  || ...
%          (~isempty(strmatch('  Stack:', line)))       || ...
%          (~isempty(strmatch('  Name:', line)))        || ...
%          (~isempty(strmatch('  Line:', line)))        || ...
%          (~isempty(strmatch('    File:', line)))      
%     %do nothing
%   elseif ~isempty(strmatch('THEY WILL BE ASSIGNED', line)) || ...
%          ~isempty(strmatch('evolucionGeneric:FUSS:', line))
%     %do nothing
%   elseif (~isempty(strmatch('BEFORE FINAL REAPING', line))) || ...
%          (~isempty(strmatch('BEFORE FINAL DRAWING', line))) || ...
%          (~isempty(strmatch('AFTER FINAL DRAWING', line)))
%     fgetl(ftst);
%     %do nothing
%   else
%     fprintf('I CANNOT UNDERSTAND THIS <%s>\n', line);
  end
end

if isempty(info(end).timeReaped)
  info = info(1:end-1);
end

%bestFitness = zeros(1,numel(info));
for k=1:numel(info)
  if feof(fres)
    break;
  end
  line = fgetl(fres);
  %info(k,1).bestFitness = sscanf(line, '%f');
  info(k,1).bestFitness = sscanf(line, '%*d %*d %*d %*s %f');
end

fclose(ftst);
fclose(fres);

if show
  bestFitness = [info(:,1).bestFitness];
  r1ts    = diff([0, [info(:,1).timeSpawned]]);
  r1tr    = diff([0, [info(:,1).timeReaped]]);
  r1tsc   = diff([0, [info(:,1).totalClusterSpawned]]);
  r1trc   = diff([0, [info(:,1).totalClusterReaped]]);
  numstrs = [info(:,1).numStrings];
  meanl   = [info(:,1).meanl];
  maxl    = [info(:,1).maxl];
  totlen  = meanl.*numstrs;

  tims = 1:maxgen;
  plot(tims, r1tsc, 'm', tims, meanl, 'r', tims, maxl, 'b', tims, totlen, 'k', tims, bestFitness, 'g'); grid on;
  legend({'t. cluster', 'mean ind. len.', 'max ind. len.', 'tot. ind. len.', 'best fitness'});
end