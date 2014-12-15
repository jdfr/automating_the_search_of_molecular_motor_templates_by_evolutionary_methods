function [rasters varargout] = evaluarMolMots3D(P,indexesP,rndSeed,returnIndividuals,evparams,taskinfo) %#ok<INUSL,INUSD>

n           = numel(P);
rasters     = cell(n,1);

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
end
if any(retInds)
  toCopy  = makeToCopy(evparams.ss);
end
for i=1:n
  if lggng;
    f=fopen(zz, 'a');
    g = taskinfo.g2s(P{i});
    fprintf(f, 'RNDSEED %s, \nGENOME: %s\nbefore computing %d', num2hex(rndSeed(i)), g,i);
    g = []; %#ok<NASGU>
    fclose(f);
  end
  structs = evaluarMolMot(P{i}, evparams, rndSeed(i), [], useKeyboard, outputrec);
  thereIsError = isfield(structs, 'ssError') || isempty(structs.stats);
  if (~thereIsError)
    stats = structs.stats;
    for k=1:nout
      varargout{k}(i) = stats.(fns{k});
    end
  end
  if lggng; f=fopen(zz, 'a'); fprintf(f, 'after computing %d\n', i); fclose(f); end
  
  if retInds(i) && not(thereIsError)
    ssFin         = structs.ssAfter;
    rasters{i}    = getSS(ssFin, toCopy);
  elseif thereIsError
    rasters{i} = struct('isError', {true}, 'genome', {P{i}}, 'evparams', evparams, 'rndSeed', {rndSeed(i)}, 'state', {structs});
  end
  if lggng; f=fopen(zz, 'a'); fprintf(f, 'after computing raster for %d\n', i); fclose(f); end
end

  if lggng; delete(zz); end; 


function toCopy = makeToCopy(ss)
toCopy = [ss.pointVars, ss.springVars, 'u', 't', 'tick', 'ndims', 'nsprings', 'npoints'];%'organism'};
for k=1:numel(toCopy)
  if iscell(toCopy{k})
    toCopy{k} = toCopy{k}{1};
  end
end
toCopy = unique(toCopy);

function res = getSS(ss, toCopy)

aux  = [toCopy; cellfun(@(x){ss.(x)}, toCopy, 'uniformoutput', false)];
res = struct(aux{:});

