%set up a plotting round
function p = initPlotting(p)
  if p.createWindow %it is sane to set this always to true
    if isempty(p.AVIfilename) %if we are not recording AVI
        if p.useOpenGL 
            p.figura = figure('Renderer', 'OpenGL', 'RendererMode', 'manual');
        else
            p.figura = figure;
        end
        p.visible = true;
        %maximize drawn area
        set(p.figura, 'Units', 'normalized');
        set(p.figura, 'Position', p.figurePos);
    else
        %while recording AVI, the window is not visible
        p.figura = figure('Visible', 'off');
        p.visible = false;
        set(p.figura, 'Units', 'normalized');
        set(p.figura, 'Position', p.figurePos);%[0 0 2.5 2.5] [1 1 round(1280*2) round(1024*2)]);
        %prepare AVI file
        if exist(p.AVIfilename, 'file')
            delete(p.AVIfilename);
        end
        p.fileAVI = avifile(p.AVIfilename, 'compression',p.codecAVI, 'fps',p.fpsAVI, 'keyframe',1, 'quality',100);
    end;
    %white background
    set(p.figura, 'Color', 'w');
  end
  p.justStarted = true;
  if p.showLineInfo
    set(p.figura,'Toolbar','figure');
    p.textLineInfo = uicontrol('Style', 'text', 'Units', 'normalized', ...
    'FontName', 'FixedWidth', 'Position', [0.1 0 0.9 0.07]);
  end
end
      
