function poph = postProcessPoph(poph, verbose, processByGenerationFun, resType)

mng = min(poph.generation);
mxg = max(poph.generation);

notFirstRasters = false;
if ~exist('resType', 'var')
  resType       = [];
end

if isfield(poph, 'rangeid')
  fname = @(k) [poph.basedir filesep sprintf('lastFrames%04d_%d.mat', k, 1)];
else
  fname = @(k) [poph.basedir filesep sprintf('lastFrames%04d_1.mat', k)];
end

if isempty(resType)
  load(fname(0));
  if isfield(res, 'rasters');        resType = 1;
  elseif isfield(res, 'rasterDevs'); resType = 2;
  else                               error('Not expected!!!');
  end
end

if isnumeric(resType)
  switch resType
    case 1
      resType = struct('unit',       {{[]}}, ...
                       'getRasters', {@(res)res.rasters});
    case 2
      resType = struct('unit',       {struct('dev', [], 'ev', [])}, ...
                       'getRasters', {@(res)struct('dev', res.rasterDevs, 'ev', res.rasterEvs)});
    otherwise
      error('resType not understood!!!!');
  end
elseif ~isstruct(resType)
  error('resType not understood!!!!');
end

rewindnans = isfield(poph.params, 'rewindNans') && poph.params.rewindNans;

if verbose
  nc = fprintf('Initializing...');
end
poph = processByGenerationFun(poph, []);

for k=mng:mxg
  load(fname(k), 'res');
  thisGen                          = find(poph.generation==k);
  notAreNew                        = true(size(thisGen));
  notAreNew(res.indexNewG)         = false;
  notAreNew                        = find(notAreNew);
  rasters                          = repmat(resType.unit, size(thisGen));
  if notFirstRasters
    indexesA                       = poph.tree.iA(poph.tree.gD==k);
    rasters(notAreNew)             = oldRasters(indexesA(notAreNew));
  end
  if numel(res.indexNewG)<numel(res.rasters)
    rasters                        = resType.getRasters(res);
  else
    rasters(res.indexNewG)         = resType.getRasters(res);
  end
  if verbose
    fprintf(repmat('\b', 1, nc));
    nc                             = fprintf('Before processing generation %d...', k);
  end
  poph                             = processByGenerationFun(poph, k, thisGen, res.indexNewG, rasters);
  if notFirstRasters && rewindnans
    gonenan                        = isnan(poph.fitness(thisGen)) & (indexesA~=0);
    rasters(poph.individual(thisGen(gonenan)))  = oldRasters(indexesA(gonenan));
  end
  oldRasters                       = rasters;
  notFirstRasters                  = true;
end

if verbose
  fprintf([repmat('\b', 1, nc) '\n']);
end
