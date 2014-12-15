function state = getPicassoJobState(jobID, terclus, params)

    % Actually, for running/pending calculation, we don't care about the state
    % of the subjobs, so treat everything the same.
    [anyRunning, anyPending, FAILED] = iGetRunningPending( jobID, terclus, params );
    
    %fprintf('Mira getPicassoJobState for jobID %s: %s %s %s\n', jobID, mat2str(anyRunning), mat2str(anyPending), mat2str(FAILED));

    if FAILED
        % Already warned in iGetRunningPending
        return
    end
    
    % Now deal with the logic surrounding these
    % Any running indicates that the job is running
    if anyRunning
        state = 'running';
        return
    end

    % We know numRunning == 0 so if there are some still pending then the
    % job must be queued again, even if there are some finished
    if anyPending
        state = 'queued';
        return
    end

    % Ensure that all tasks have the right state
    if ~anyRunning && ~anyPending
        state = 'finished';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iGetRunningPending - return the number of jobs/tasks that are running or
% pending. Uses "qstat -f" on each job ID in turn. The cunning use of
% qselect doesn't work because you need multiple calls, and the state of a
% job can change in between those calls, so this method is more
% conservative, and ensures that it retrieves the state of each PBS job
% precisely once.
function [anyRun, anyPend, FAILED] = iGetRunningPending( jobID, terclus, params )

anyRun  = false;
anyPend = false;

id = jobID;
areNums = ((jobID>='0') & (jobID<='9')) | (jobID=='[') | (jobID==']');
first   = find(~areNums, 1, 'first');
if ~isempty(first)
  jobID   = jobID(1:(first(1)-1));
end


    cmdLine = sprintf( 'qstat -f "%s"', jobID );
%fprintf('querying for job %s (%s) with line <<%s>>...\n', id, jobID, cmdLine);
    
    if terclus.remoteSSH
      
      channel = createSSHChannel(params, terclus);

      [channel, out1, out2]  =  executeCommand(channel,cmdLine);
      
      out = [out1 out2];
      clear out1 out2;
      
      FAILED = isempty(out) || all(isspace(out));

      channel.close;
      
    else
      
      [FAILED, out] = system( cmdLine );
      
    end
%fprintf('output to the query: <<%s>>...\n', out);
    
    if ~isempty( strfind( out, 'Unknown Job Id' ) )
        FAILED = 0;
        % Get here for a job that PBS has no knowledge of. 
        return;
    end
    
    if FAILED
        % Some other problem with qstat -f
        warning( 'distcomp:pbsscheduler:UnableToQueryState', ...
                 'Error executing the PBS command ''%s''. The reason given is \n %s', ...
                 cmdLine, out );
        return
    end
    stateLine = regexp( out, 'job_state = [A-Z]', 'match', 'once' );
    if isempty( stateLine )
        warning( 'distcomp:pbsscheduler:FailedToParseQstat', ...
                 ['Couldn''t interpret job state from the output of qstat -f.\n', ...
                  'qstat returned:\n%s'], out );
    else
        stateLetter = stateLine(end);
        
        % Use the scheduler-specific state indicators to see if this particular
        % letter indicates running or pending.
        if ~isempty( strfind( 'BRE', stateLetter ) )
            anyRun  = true;
        elseif ~isempty( strfind( 'HQSTUW', stateLetter ) )
            anyPend = true;
        end
    end
    
