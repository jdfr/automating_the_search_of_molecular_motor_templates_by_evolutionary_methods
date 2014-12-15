%this function draws the required timesteps Ts
function p = drawSimulation(p, rec, Ts)
  p = initPlotting(p);
  playSimulation(rec, Ts, @wrapper);
  p = endPlotting(p);
  function wrapper(t, st)
    p = initFrame(p, st);
    drawFrame(p, t, st);
    p = endFrame(p);
  end
end


