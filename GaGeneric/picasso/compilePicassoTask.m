function compilePicassoTask(basedir, fname, toDeploy)

mx = ['.' mexext];

mxf = getFiles(basedir, mx);

mxf = cellfun(@(x) ['-a ' x ' '], mxf, 'uniformoutput', false);
mxf = [mxf{:}];

%comm = ['mcc -v -R -nojvm -m ' fname ' -o ..' filesep 'deployed' filesep fname '_jaja ' mxf];
comm = ['mcc -v -R -nojvm -m ' fname ' -d ' toDeploy ' '  mxf];

eval(comm);



function files = getFiles(basedir, mx)
d=dir(basedir);
files = {};
for k=1:numel(d)
  nm = d(k).name;
  if d(k).isdir
    if not(any(strcmp(nm, {'.', '..'}))) && not(nm(1)=='+')
      files = [files; getFiles([basedir filesep nm], mx)]; %#ok<AGROW>
    end
  else
    if (numel(nm)>numel(mx)) && all(nm(end-numel(mx)+1:end)==mx)
      files = [files; nm]; %#ok<AGROW>
    end
  end
end