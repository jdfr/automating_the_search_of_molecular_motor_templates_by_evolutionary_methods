function recordingCallback(rec, varargin) %, ss, Ts, states)
  if ~isempty(rec.callback) && ishandle(rec.callback)
    handle = rec.callback;
    handle(varargin);
  end
end
