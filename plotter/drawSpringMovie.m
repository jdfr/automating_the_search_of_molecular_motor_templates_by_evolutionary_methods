function drawSpringMovie(AVIfilename, codec, fpsAVI, rec, Ts, imageLimits, imageDimensions)

lastNumChars     = 0;
numFramesPlotted = 0;

fileAVI = avifile(AVIfilename, 'compression',codec, 'fps', fpsAVI, 'keyframe',1, 'quality',100);

try 
  playSimulation(rec, Ts, @drawFrame);
catch ME 
  fileAVI = close(fileAVI);
  rethrow(ME);
end
fileAVI = close(fileAVI);

  function drawFrame(t, st)

    mappingHue = (updateMod(t, st.rmod, st.rtb, st.rtm)-1)*2;
    colors = squeeze(hsv2rgb(mappingHue*(2/3), ones(size(st.rmod)), ones(size(st.rmod))));
    colors(colors<0) = 0;
    colors(colors>1) = 1;
    image = drawAsImage(st.pos, st.springEnds, imageLimits, imageDimensions, colors);
    frame = im2frame(double(image));
    fileAVI = addframe(fileAVI,frame);

    numFramesPlotted = numFramesPlotted+1;

    if mod(numFramesPlotted, 10)==0%50)==0
      if lastNumChars>0
        fprintf(repmat('\b', 1, lastNumChars));
      end
      lastNumChars = fprintf('frames plotted: %10d\n', numFramesPlotted);
    end

  end

end



