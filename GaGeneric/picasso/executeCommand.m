function [channel2, result, reserr]  =  executeCommand(channel2,command)
% modification of sshfrommatlabissue to also recover stderr from the remotely invoked command

  import java.io.BufferedReader;
  import java.io.IOException;
  import java.io.InputStream;
  import java.io.InputStreamReader;
  import ch.ethz.ssh2.Connection;
  import ch.ethz.ssh2.Session;
  import ch.ethz.ssh2.StreamGobbler;
%
% Invocation checking
%
  if(nargin  ~=  2)
    error('Error: two input arguments are required...');
  end
  if(~isa(channel2,'ch.ethz.ssh2.Connection'))
    error('Error: input argument CHANNEL2 is not a Java Connection object...');
  end
  if(~ischar(command))
    error('Error: input argument COMMAND is not a string...');
  end
% 
% Send the commands
%
  result  =  {''};
  reserr  =  {''};
  channel  =  channel2.openSession();
  channel.execCommand(command);
%
% Report the result to screen and to the string result...
%  
  stdout = StreamGobbler(channel.getStdout());
  stderr = StreamGobbler(channel.getStderr());
  br = BufferedReader(InputStreamReader(stdout));
  while(true)
    line = br.readLine();
	  if(isempty(line))
	    break
    else
      if(isempty(result{1}))
        result{1}  =  char(line);
      else
        result{end+1}  =  char(line);
      end
	  %fprintf(1,'\n%s',char(line));
    end
  end
  
  br = BufferedReader(InputStreamReader(stderr));
  while(true)
    line = br.readLine();
	  if(isempty(line))
	    break
    else
      if(isempty(result{1}))
        reserr{1}  =  char(line);
      else
        reserr{end+1}  =  char(line);
      end
	  %fprintf(1,'\n%s',char(line));
    end
  end
  
  channel.close();
  result  =  result';
  reserr  =  reserr';
  
  if ~isempty(result)
    ln = sprintf('\n');
    result = result(:)';
    result(2,:) = {ln};
    result = horzcat(result{:});
  end
  
  if ~isempty(reserr)
    ln = sprintf('\n');
    reserr = reserr(:)';
    reserr(2,:) = {ln};
    reserr = horzcat(reserr{:});
  end
%   fprintf('\nMIRA: <<%s>>, <<%s>>, <<%s>>\n', command, result, result2);

