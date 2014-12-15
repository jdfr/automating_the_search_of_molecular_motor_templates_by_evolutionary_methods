function [files] = filenamedirrec(pattern)
    if(nargin == 0)
        pattern='./*';
        posslash=0;
    else
        posslash=findstr('/',pattern);
        if(length(posslash > 0))
            posslash=posslash(end);
        else
            posslash=0;
        end
    end
    
    if(pattern(end) == '/')
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
    
    ndirr = size(dirr, 1);
    for i=1:ndirr
        file=dirr(i);
        if(file.isdir && ~strcmp(file.name,'.') && ~strcmp(file.name,'..'))
            newfiles = filenamedirrec([pattern(1:posslash) file.name '/' pattern(posslash+1:end)]);
            files = {files{:}, newfiles{:}};
        end
    end
    