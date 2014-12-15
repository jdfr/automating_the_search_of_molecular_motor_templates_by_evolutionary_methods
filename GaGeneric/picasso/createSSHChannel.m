function channel = createSSHChannel(params, terclus)
  strs = javaArray('java.lang.String', 1);
  strs(1) = java.lang.String(terclus.ssh.password); x.responses = strs;
  x = jd.PreloadedInteractive;
  x.responses = strs;
  x.showThings = false;
  done  = false;
  if isfield(terclus, 'pauseTimeAfterError')
    pt = terclus.pauseTimeAfterError;
  else
    pt = 10;
  end
  if isfield(terclus, 'submitTimes')
    st = terclus.submitTimes;
  else
    st = 3;
  end
  channel = [];
  for z=1:st
    try
      channel = ch.ethz.ssh2.Connection(terclus.ssh.nameurl);
      info = channel.connect;
      done = true;
    catch ME
      mprintf(params.ftst, 'Try %d. There has been an error while trying to open a ssh channel:\n%s\n', z, showError(ME));
      try
        if ~isempty(channel)
          channel.close;
        end
      catch ME2
        mprintf(params.ftst, 'Error while trying to close a bad channel:\n%s\n', showError(ME2));
      end
      if z<st
        mprintf(params.ftst, 'Waiting %s seconds... ', mat2str(pt));
        pause(pt);
        mprintf(params.ftst, 'Now, trying again...\n');
      else
        mprintf(params.ftst, 'After %d tries spaced %s seconds between them, we have not been able to open a channel. Simulation will abort...\n', z, mat2str(pt));
        rethrow(ME);
      end
    end
    if done
      break
    elseif z==st
      mprintf(params.ftst, ':( After %d tries spaced %s seconds between them, we have not been able to open a channel. Simulation will abort...\n', z, mat2str(pt));
    end
  end
  ok = channel.authenticateWithKeyboardInteractive(terclus.ssh.username, x);
  if ~ok
    channel.close;
    error('Authentication was not possible!!!');
  end
