%finalize frame drawing
function p = endFrame(p)
  if ~isempty(p.AVIfilename) %record AVI frame
      if p.visible
          F = getframe(get(p.figura,'CurrentAxes'));
          p.fileAVI = addframe(p.fileAVI,F);
      else
          p.fileAVI = addframe(p.fileAVI,get(p.figura,'CurrentAxes'));
      end
  else
      %dump new graphics to screen
      drawnow;
      %do pause, if needed
      if isempty(p.pauseTime)
        pause;
      elseif p.pauseTime>0
        pause(p.pauseTime);
      end
  end;
end
    
