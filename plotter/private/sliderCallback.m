    %callback for sliding through INDEXES 
    function sliderCallback(p, rec, Ts, objs)
      [edTime, edIndex, slider] = objs{:};
      idx = round(get(slider, 'Value')); %get INDEX
        set(edTime,  'String', num2str(Ts(idx)));
        set(edIndex, 'String', num2str(idx));
        set(slider,  'Value',  idx);
        playSimulation(rec, Ts(idx), @(t, st) frameWrapper(p, t, st));
%       updateControls(p, rec, Ts, idx, edTime, edIndex, slider);
    end