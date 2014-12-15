%finalize plotting round
function p = endPlotting(p)
  if ~isempty(p.AVIfilename) %close AVI
      p.fileAVI = close(p.fileAVI);
      close all hidden;
  end;
end

