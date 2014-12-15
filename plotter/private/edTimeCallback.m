    %callback for specifying a TIME to see 
    function edTimeCallback(p, rec, Ts, objs)
      [edTime, edIndex, slider] = objs{:};
      t = str2double(get(edTime, 'String')); %parse and check TIME
      if ~isnan(t) && ~isinf(t) && (t>=min(Ts)) && (t<=max(Ts))
        [nevermind idx] = min(abs(Ts-t)); %get the INDEX
        idx = idx(1); %keep it safe==vayamos a pollas
        set(edTime,  'String', num2str(Ts(idx)));
        set(edIndex, 'String', num2str(idx));
        set(slider,  'Value',  idx);
        playSimulation(rec, Ts(idx), @(t, st) frameWrapper(p, t, st));
%         updateControls(p, rec, Ts, idx, edTime, edIndex, slider);
      end
    end

