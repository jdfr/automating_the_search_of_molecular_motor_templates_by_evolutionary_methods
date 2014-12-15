function [rasterDevs rasterEvs varargout] = evaluarMolMots(P,indexesP,rndSeedDev,rndSeedEv,returnIndividuals,evparams,taskinfo) %#ok<INUSL,INUSD>


n           = numel(P);
rasterDevs  = cell(n,1);
rasterEvs   = cell(n,1);

nout        = size(evparams.statNames,2);

for k=1:nout
  varargout{k} = repmat(evparams.statNames{2,k}, n, 1);
end
fns = evparams.statNames(1,:);

retInds     = returnIndividuals>0;
useKeyboard = false;
outputrec   = false;
lggng = ~isempty(taskinfo.dir); 
if lggng;
  zz = [taskinfo.dir sprintf('zlog_%d_%d_%d.txt', taskinfo.ngen, taskinfo.rid, taskinfo.ntask)];
  [f msg]=fopen(zz, 'w');
  if ~isempty(msg)
    error('We have not been able to open <%s>\nCause: <%s>\n', zz, msg);
  end
  fprintf(f, 'NUMBER OF GENOMES TO EVALUATE: %d\n', n);
  fclose(f);
  toCopyStr = makeToCopy(evparams.ss, 'str');
end
if any(retInds)
  toCopyDev = makeToCopy(evparams.ss, 'dev');
  toCopyEv  = makeToCopy(evparams.ss, 'ev');
end
for i=1:n
  if lggng;
    f=fopen(zz, 'a');
    if isstruct(P{i})
      g = any2str(getSS(P{i}, toCopyStr));
    else
      g = taskinfo.g2s(P{i});
    end
    fprintf(f, 'RNDSEED DEV %s, EV %s\nGENOME: %s\nbefore computing %d', num2hex(rndSeedDev(i)), num2hex(rndSeedEv(i)), g,i);
    g = []; %#ok<NASGU>
    fclose(f);
  end
  structs = evaluarMolMot(P{i}, evparams, rndSeedDev(i), rndSeedEv(i), useKeyboard, outputrec);
  thereIsError = isfield(structs, 'ssError') || isempty(structs.stats);
  if (~thereIsError)
    stats = structs.stats;
    for k=1:nout
      varargout{k}(i) = stats.(fns{k});
    end
  end
  if lggng; f=fopen(zz, 'a'); fprintf(f, 'after computing %d\n', i); fclose(f); end
  
  if retInds(i) && not(thereIsError)
    ssFin = structs.ssAfter;
%     if evparams.relax.doRelax
%       ssFin = structs.ssRelaxed;
%     else
%       ssFin = structs.ssAfter;
%     end
    rasterEvs{i}  = getSS(ssFin,          toCopyEv);
    rasterDevs{i} = getSS(structs.ssDev,  toCopyDev);
  elseif thereIsError
    rasterDevs{i} = struct('isError', {true}, 'genome', {P{i}}, 'evparams', evparams, 'rndSeedDev', {rndSeedDev(i)}, 'rndSeedEv', {rndSeedEv(i)}, 'state', {structs});
  end
  if lggng; f=fopen(zz, 'a'); fprintf(f, 'after computing raster for %d\n', i); fclose(f); end
end

  if lggng; delete(zz); end; %f=fopen(zz, 'a'); fprintf(f, 'about to delete this file\n'); fclose(f); delete(zz); end


function toCopy = makeToCopy(ss, mode)
switch mode
  case 'dev'
    %toCopy = [ss.pointVars, ss.springVars, ss.atomVars, ss.pShiftingVars, 't', 'tick', 'ndims', 'nsprings', 'npoints', 'atomConstants'];%'organism'};
    %these vars are stored in the rasterDev structure because they are
    %either: (a) costly to be calculated (b) require an unreasonably high
    %amount of other variables to be calculated
    bufferedVars = {'ncellsTeo', 'overload1', 'developed'};
    toCopy = [ss.pointVars, ss.springVars, ss.atomVars, bufferedVars{:}, 't', 'tick', 'ndims', 'nsprings', 'npoints', 'atomConstants'];%'organism'};
  case 'ev'
    toCopy = [ss.pointVars, ss.springVars, 't', 'tick', 'ndims', 'nsprings', 'npoints'];%'organism'};
  case 'str'
    toCopy = ['t', 'u', 'tick', 'pos', 'vel', 'm', 'springEnds', 'k', 'r'];%'organism'};
end
for k=1:numel(toCopy)
  if iscell(toCopy{k})
    toCopy{k} = toCopy{k}{1};
  end
end
toCopy = unique(toCopy);

function res = getSS(ss, toCopy)

aux  = [toCopy; cellfun(@(x){ss.(x)}, toCopy, 'uniformoutput', false)];
res = struct(aux{:});
% for k=1:numel(toCopy)
%   res.(toCopy{k}) = ss.(toCopy{k});
% end

  
  
