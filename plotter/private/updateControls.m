    %update controls in navigation mode
    function updateControls(p, rec, Ts, idx, edTime, edIndex, slider)
      set(edTime,  'String', num2str(Ts(idx)));
      set(edIndex, 'String', num2str(idx));
      set(slider,  'Value',  idx);
%       set(slider,  'TooltipString', sprintf('SLIDER control. T:%g, i:%g', Ts(idx), idx));
      playSimulation(rec, Ts(idx), @(t, st) frameWrapper(p, t, st));
    end
