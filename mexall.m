function mexall(folder, changeAtEnd)

if nargin<2
  changeAtEnd = true;
end

if changeAtEnd
  old = pwd;
end
d = dir(folder);
nams = {d.name};
for k=1:numel(d)
  nam = nams{k};
  if d(k).isdir
    if all(nam(1)~='+@.')
      mexall([folder filesep nam], false);
    end
  else
    if ( (numel(nam)>2) && all(lower(nam(end-1:end))=='.c') ) || ...
       ( (numel(nam)>4) && all(lower(nam(end-3:end))=='.cpp') )
      lastp       = find(nam=='.',1,'last');
      mexname     = [nam(1:lastp) mexext];
      mexk        = find(strcmp(mexname, nams));
      switch numel(mexk)
        case 0
          doit = true;
        case 1
          if d(mexk).isdir %#ok<FNDSB>
            error('The mex file is a directory!!!!!');
          end
          doit = d(k).datenum > d(mexk).datenum;
          %fprintf('datenum  .c: %s\ndatenum mex: %s\n', mat2str(d(k).datenum), mat2str(d(mexk).datenum));
        otherwise
          error('This is plainly weird!!!');
      end
      if doit
        cd(folder);
        fprintf('Compiling %s...\n', [folder filesep nam]);
        try
          mex(nam);
        catch ME
          fprintf('The file was not compiled: %s\n', showError(ME));
          %break;
        end
      else
        fprintf('Found %s, but compilation is not necessary\n', [folder filesep nam]);
      end
    end
  end
end

if changeAtEnd
  cd(old);
end



