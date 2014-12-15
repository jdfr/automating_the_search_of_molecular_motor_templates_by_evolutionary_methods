function [files] = my_filenamedirrec(pattern)
%the same as filenamedirrec, but does not include files from directories
%starting with '+', because we do not want them to be sent to the cluster
    if(nargin == 0)
        pattern=['.' filesep '*'];
        posslash=0;
    else
        posslash=findstr(filesep,pattern);
        if(~isempty(posslash))
            posslash=posslash(end);
        else
            posslash=0;
        end
    end
    
    if(pattern(end) == filesep)
        pattern = [pattern '*'];
    end
    
    
    files={};
    
    dirr = dir(pattern);
    
    ndirr = size(dirr, 1);
    for i=1:ndirr
        file=dirr(i);
        if(file.isdir == 0)
            file.path = pattern(1:posslash);
            files{end+1} = [file.path file.name];
        end
    end
    
    dirr = dir(pattern(1:posslash));
    
    xx=cell(numel(dirr),1);
    for k=1:numel(dirr)
      xx{k}=dirr(k).name;
    end
    
    ndirr = size(dirr, 1);
    for i=1:ndirr
        file=dirr(i);
        if(file.isdir && ~strcmp(file.name,'.') && ~strcmp(file.name,'..') && (file.name(1)~='+') && (file.name(1)~='_'))
            newfiles = my_filenamedirrec([pattern(1:posslash) file.name filesep pattern(posslash+1:end)]);
            files = {files{:}, newfiles{:}};
        end
    end
    