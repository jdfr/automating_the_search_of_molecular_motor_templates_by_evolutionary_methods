%utility function for drawing a frame
function frameWrapper(p, t, st)
  p = initFrame(p, st);
  drawFrame(p, t, st);
  p = endFrame(p);
end
