    %callback for specifying an INDEX to see 
    function edIndexCallback(p, rec, Ts, objs)
      [edTime, edIndex, slider] = objs{:};
      idx = str2double(get(edIndex, 'String')); %parse and check INDEX
      if ~isnan(idx) && ~isinf(idx) && (idx>=get(slider, 'Min')) && (idx<=get(slider, 'Max'))
        idx = round(idx);
        set(edTime,  'String', num2str(Ts(idx)));
        set(edIndex, 'String', num2str(idx));
        set(slider,  'Value',  idx);
        playSimulation(rec, Ts(idx), @(t, st) frameWrapper(p, t, st));
%         updateControls(p, rec, Ts, idx, edTime, edIndex, slider);
      end
    end
