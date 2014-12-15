function saveCompiledExamples(conf, name, selecs, freq)

if ~exist('freq', 'var')
  %freq = 20;
  error('frequency must be provided!!!!');
end

%ndata = numel(compileExamples('data'));

for z=1:numel(selecs);%ndata
  k = selecs(z);
  results = compileExamples('sims', conf, k, false, true, freq); 
  path    = results.path; 
  sps     = find(path==filesep); 
  nam     = path(sps(end-1)+1:sps(end)-1); 
  pth     = sprintf(name, nam);
  fprintf('saving simulation <%s> in <%s> ...\n', nam, pth); 
  save(pth, 'results');
  clear results;
end