classdef BasicPlotter < hgsetget
  
  properties
    figura;      %figure's handle
    figurePos;   %figure's position in normalized units
    figureUnits; %units for expressing figure's position
    
    AVIfilename; %name of AVI file to record
    AVImethod;   %Two options: 'MATLAB' (use MATLAB infrastructure, mostly
                 %for easy use in Windows), 'VIDEOIO' (use videoIO library,
                 %mostly for use in UNIX, please take care that ffmpeg DEV
                 %packages are installed (in DEBIAN: libavcodec-dev
                 %libpostproc-dev libavformat-dev libswscale-dev), and that
                 %videoIO is also installed, compiled (it might be a
                 %pain in the arse to get it compiled) and in the MATLAB
                 %path)
    createWindow;%this should be always true, unless you really know what you're doing 
    useOpenGL;   %whether to use openGl. Enfatically recommended when not rendering movies
    visible;     %info flag for figure's visibility. It is automatically set (AVI=not visible) 
    AVIWidthHeight; %only used for AVImethod=='VIDEOIO': if not empty, it must be a 2-vector specifying the movie frame's width and height
    
    fileAVI;     %AVI file's handle
    codecAVI;    %AVI file's codec
    fpsAVI;      %AVI file's fps
    
    fixedAxisFrame; %whether the axes lengths are fixed. It is implied by 'axisWindow' 
    axisWindow;  %if []=>not fixed frames. if num=>squared axes'semiradius. if vector=>axes
    hasGrid;     %whether to use a grid
    hasAxis;     %whether to display coordinate axes
    axisEqual;   %whether to use 'axis equal'
    axisSquare;  %whether to use 'axis square'
    hasTitle;    %whether to put a title
    pauseTime;   %after each frame: []=>do pause. 0=>no pause. n=>pause n seconds 
    viewPoint;   %[]=> standard viewpoint. T=>transformation matrix for 'view' 
    axesShown;   %axes to be shown
    axesLengths; %lengths of axes
    axesOrigins; %origin of axes
    axesColors;  %colors of axes
    axesArgs;    %arguments for axes' lines
    camProps;    %camera properties
    camMode;     %mode of camera setting
    handleImage; %do whatever you want with the image
  end
  
  properties (GetAccess = protected, SetAccess = protected)
    justStarted; %flag to tell if the frame is the first to be drawn
    lastNumChars;     %to be used when making an AVI
    numFramesPlotted; %to be used when making an AVI
    defaultCodecs;
  end
  
  methods
    function p = BasicPlotter(varargin)
      %set defaults
      p.figureUnits = 'normalized';
      p.figurePos   = [0.05 0.05 0.8 0.8];
      
      p.AVIfilename     = '';
      p.AVImethod       = 'MATLAB'; %'VIDEOIO'
      p.fileAVI         = [];
      p.createWindow    = true;
      p.useOpenGL       = true;
      p.visible         = true;
      p.AVIWidthHeight  = [];
      
      p.defaultCodecs = {'CinePak', 'mpeg4'};
      p.codecAVI      = '';
      p.fpsAVI        = 15;
      
      p.fixedAxisFrame = true;
      p.axisWindow = 2;
      p.hasGrid    = true;     
      p.hasAxis    = true;     
      p.axisEqual  = true;   
      p.axisSquare = false;  
      p.hasTitle   = true;
      p.pauseTime  = 0;
      p.viewPoint  = [];
      
      p.justStarted = [];
      p.lastNumChars = 0;
      p.numFramesPlotted = 0;
      p.axesShown   = '';
      p.axesLengths = [];
      p.axesColors  = '';
      p.axesOrigins  = [];
      p.axesArgs     = {};
      
      p.handleImage  = [];
      
      p.camProps     = {};
      p.camMode      = 'none'; %'none'; %'before'; %'after';
      %set customized options
      if numel(varargin)>0; set(p, varargin{:}); end
    end

    %changing axisWindow implies changing fixedAxisFrame
    function p = set.axisWindow(p, v)
      if isempty(v)
        p.fixedAxisFrame = false;
        p.axisWindow = v;
      else
        p.fixedAxisFrame = true;
        p.axisWindow = v;
      end
    end
    
    %set up a plotting round
    function initPlotting(p)
      if p.createWindow %it is sane to set this always to true
        if isempty(p.AVIfilename) %if we are not recording AVI
            if p.useOpenGL 
                p.figura = figure('Renderer', 'OpenGL', 'RendererMode', 'manual');
            else
                p.figura = figure;
            end
            p.visible = true;
            %maximize drawn area
            set(p.figura, 'Units',    p.figureUnits, ...
                          'Position', p.figurePos);
        else
            %while recording AVI, the window is not visible
            p.lastNumChars = 0;
            p.numFramesPlotted = 0;
            p.visible = false;
            p.figura = figure('Visible',  'off', ...
                              'Units',    p.figureUnits, ...
                              'Position', p.figurePos);%[0 0 2.5 2.5] [1 1 round(1280*2) round(1024*2)]);
            %prepare AVI file
            if exist(p.AVIfilename, 'file')
                delete(p.AVIfilename);
            end
            p.AVImethod = upper(p.AVImethod);
            switch p.AVImethod
              case 'MATLAB'
                if isempty(p.codecAVI)
                  codec = p.defaultCodecs{1};
                else
                  codec = p.codecAVI;
                end
%                 puold = get(p.figura, 'PaperUnits');
%                 set(p.figura, 'PaperUnits', p.figureUnits, 'PaperPosition', p.figurePos);
%                 set(p.figura, 'PaperUnits', puold);
                p.fileAVI = avifile(p.AVIfilename, 'compression',codec, 'fps',p.fpsAVI, 'keyframe',1, 'quality',100);
              case 'VIDEOIO'
                if isempty(p.AVIWidthHeight) %frame's width and height are chosen to be the same as figure's ones
                  p.figureUnits = 'pixels';
                  set(p.figura, 'Units',    p.figureUnits);
                  pos = round(get(p.figura, 'Position'));
                  set(p.figura, 'Position', pos);
                  widthheight = pos(3:4);
                else
                  widthheight = p.AVIWidthHeight;
                end
                if isempty(p.codecAVI)
                  codec = p.defaultCodecs{2};
                else
                  codec = p.codecAVI;
                end
                p.fileAVI = videoWriter(p.AVIfilename, 'codec', codec, 'width',widthheight(1), 'height',widthheight(2), 'fps',p.fpsAVI);
              case 'IMAGES'
              otherwise
                error('AVI method not recognized: <%s>', p.AVImethod);
            end
        end;
        %white background
        set(p.figura, 'Color', 'w');
      end
      p.justStarted = true;
    end
      
    %set up a frame for drawing it
    function initFrame(p, st)
      %erase axes' contents
      cla(get(p.figura, 'CurrentAxes'));
      if ~isempty(p.AVIfilename)
        set(p.figura, 'Position', p.figurePos);
      end
      hold on;                               
      %set up axes
      if p.fixedAxisFrame
        if numel(p.axisWindow)>1
          axis(p.axisWindow);
        else
          switch size(st.pos,2)
            case 2; axis([-p.axisWindow p.axisWindow -p.axisWindow p.axisWindow]);
            case 3; axis([-p.axisWindow p.axisWindow -p.axisWindow p.axisWindow -p.axisWindow p.axisWindow]);
          end
        end
      end;
      %the first time we set the viewpoint
      if p.justStarted 
        if ~isempty(p.viewPoint)
          view(p.viewPoint);
        else
          switch size(st.pos,2)
            case 2; view(2);
            case 3; view(3);
          end
        end
        p.justStarted = false;
      end
      if (~isempty(p.camProps)) && strcmp(p.camMode, 'before')
        cameraProps(gca, p.camProps);
      end
% units = get(gca,'Units'); set(gca,'Units','Pixels')
% a = get(gca,'Position'); set(gca,'Units',units)
% set(gca,'DataAspectRatio',[1 1 1]);
% dx = diff(get(gca,'xlim')); dy = diff(get(gca,'ylim'));
% dz = diff(get(gca,'zlim'));
% error('Mira esto: dx %s - dy %s - dz %s - a %s', mat2str(dx), mat2str(dy), mat2str(dz), mat2str(a));
% error('Mira esto: %s - %s', mat2str(get(gca, 'PlotBoxAspectRatio')), mat2str(get(p.figura, 'Position')));
      if p.axisEqual;  axis equal;  end
      if p.axisSquare; axis square; end
      if p.hasAxis; axis on; else axis off; end;
      if p.hasGrid; grid on; else grid off; end;
    end
    
    %finalize frame drawing
    function endFrame(p)
      if (~isempty(p.camProps)) && strcmp(p.camMode, 'after')
        cameraProps(gca, p.camProps);
      end
      if ~isempty(p.AVIfilename) %record AVI frame
        switch p.AVImethod
          case 'MATLAB'
            if p.visible
                F = getframe(get(p.figura,'CurrentAxes'));
                p.fileAVI = addframe(p.fileAVI,F);
            else
                p.fileAVI = addframe(p.fileAVI,get(p.figura,'CurrentAxes'));
            end
          case 'VIDEOIO'
            addframe(p.fileAVI, p.figura);
          case 'IMAGES'
            img = imcapture(get(p.figura, 'CurrentAxes'), 'img', p.AVIWidthHeight);
            if ~isempty(p.handleImage);
              p.handleImage(p.AVIfilename, img, p.numFramesPlotted);
            end
        end
        p.numFramesPlotted = p.numFramesPlotted+1;
        if mod(p.numFramesPlotted, 1)==0%50)==0
          if p.lastNumChars>0
            fprintf(repmat('\b', 1, p.lastNumChars));
          end
          p.lastNumChars = fprintf('frames plotted: %10d', p.numFramesPlotted);
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
    
    %finalize plotting round
    function endPlotting(p)
      if ~isempty(p.AVIfilename) && any(strcmp(p.AVImethod, {'MATLAB', 'VIDEOIO'})) %close AVI
        %FOR BOTH METHODS 'MATLAB' AND 'VIDEOIO'
          p.fileAVI = close(p.fileAVI);
          close(p.figura);
          p.fileAVI = [];
      end;
    end

    %subclasses must redefine this
    function drawFrame(p, t, st)
      if ~isempty(p.axesShown)
        axisWindow           = p.axisWindow;
        axesShown            = p.axesShown;
        axesOrigins          = p.axesOrigins;
        axesLengths          = p.axesLengths;
        axesColors           = p.axesColors;
        axesArgs             = p.axesArgs;
        if ~iscell(axesShown)
          axesShown = {axesShown};
        end
        if ~iscell(axesOrigins)
          axesOrigins = {axesOrigins};
        end
        if ~iscell(axesColors)
          axesColors = repmat({axesColors}, size(axesShown));
        end
        if isempty(axesArgs)
          axesArgs = repmat({{{},{},{}}}, size(axesShown));
        elseif (~iscell(axesArgs{1})) || (~iscell(axesArgs{1}{1}))
          axesArgs = repmat({axesArgs}, size(axesShown));
        end
        if ~iscell(axesLengths)
          axesLengths = repmat({axesLengths}, size(axesShown));
        end
        for k=1:numel(axesShown)
          shown    = axesShown{k};
          lengths  = axesLengths{k};
          colors   = axesColors{k};
          if isempty(colors)
            colors = ('kkk')';
          end
          if ischar(colors)
            colors = colors(:);
          end
          args     = axesArgs{k};
          if ~isempty(args) && ~iscell(args{1})
            args   = repmat({args}, size(shown));
          end
          origin   = axesOrigins{k};
          if isempty(origin)
            origin = [0 0 0];
          end
          if isa(origin, 'function_handle')
            origin = origin(st);
          end
          origin   = arrayfun(@(x)[x x], origin, 'uniformoutput', false);
          if isempty(lengths)
            if numel(axisWindow)>1
              lengths      = reshape(axisWindow, 2, [])';
            else
              lengths      = repmat(axisWindow, 3, 2);
              lengths(:,1) = -lengths(:,1);
            end
          elseif numel(lengths)==numel(shown)
            lengths = [-lengths(:), lengths(:)];
          end
          for m=1:numel(shown)
            z       = shown(m);
            lpos    = origin;
            lpos{z} = lpos{z}+lengths(z,:);
            line(lpos{:}, 'Color', colors(z,:), args{z}{:});
          end
        end
      end
      if p.hasTitle
        title(sprintf('time: %-10g; e:%-10g', t, calculatePotentialEnergy(st) ));% calculateEnergy(st)));
      end
    end
    
    %utility function for drawing a frame
    function frameWrapper(p, t, st)
      initFrame(p, st);
      drawFrame(p, t, st);
      endFrame(p);
    end
    
    function drawSS(p, st)
      if isempty(st)
        fprintf('I cannot draw an empty thing!!!');
      else
        initPlotting(p);
        frameWrapper(p, st.t, st);
        endPlotting(p);
      end
    end
    
    %this function draws the required timesteps Ts
    function drawSimulation(p, rec, Ts)
      initPlotting(p);
      try 
        playSimulation(rec, Ts, @(t, st) frameWrapper(p, t, st));
      catch ME
        endPlotting(p);
        rethrow(ME);
      end
      endPlotting(p);
    end

    %this function is a helper to draw the states as they are passed down 
    %to the Recorder (the args {ss, Ts, states} are passed to the Recorder
    %callback, see Recorder.recordingCallback).
    function recordFrames(p, ss, Ts, states)
      if nargin<3 %full state change
        frameWrapper(p, ss.t, ss);
      else %minor change
        for k=1:numel(Ts);
          ssOtro = dumpEvolvedState(ss, states(:,k));
          ssOtro.t = Ts(k);
          frameWrapper(p, ssOtro.t, ssOtro);
        end
      end
    end    
    
    %this function navigates through the required timesteps
    %returns info for automating navigation
    %please DO NOT USE THIS while recording AVI files! ERRORS ARE GUARANTEED TO HAPPEN! 
    function navigatorInfo = navigateSimulation(p, rec, Ts)
      if isempty(rec)
        fprintf('I cannot show an empty simulation record!!!');
        navigatorInfo = [];
        return
      end
      if nargin<3
        Ts = allTimeSteps(rec);
      end
      initPlotting(p);
      set(p.figura,'Toolbar','figure');
      commonProps = {'Units', 'normalized', 'Interruptible', 'off', 'BusyAction', 'cancel'};
      %edit control: sets the time
      edTime = uicontrol('Style', 'edit', commonProps{:}, ...
        'Position', [0 0.93 0.1 0.07], 'String', sprintf('%g', Ts(1)), ...
        'TooltipString', sprintf('TIME control. minT:%g, maxT:%g', min(Ts), max(Ts)));
      %edit control: sets the index
      edIndex = uicontrol('Style', 'edit', commonProps{:}, ...
        'Position', [0 0 0.1 0.07], 'String', '1', ...
        'TooltipString', sprintf('INDEX control. min:1, max:%g', numel(Ts)));
      %slider control: slides over indexes
      slider = uicontrol('Style', 'slider', commonProps{:}, ...
        'Min', 1, 'Max', numel(Ts), 'SliderStep', [1, 10]/numel(Ts), ...
        'Position', [0 0.07 0.02 0.86], 'Value', 1);
      %inexplicably, our callbacks cannot use more than 4 args. Thus, we pack them 
      lastarg = {edTime, edIndex, slider};
      %set the callbacks. If they are not static, error pop up everywhere
      set(slider,  'Callback', @(obj, event)sliderCallback(p, rec, Ts, lastarg));
      set(edTime,  'Callback', @(obj, event)edTimeCallback(p, rec, Ts, lastarg));
      set(edIndex, 'Callback', @(obj, event)edIndexCallback(p, rec, Ts, lastarg));
      %start it!
      sliderCallback(p, rec, Ts, lastarg);
      %store results
      if nargout==1
        navigatorInfo = struct('rec', rec, 'Ts', Ts, 'edTime', edTime, 'edIndex', edIndex, 'slider', slider, 'p', p);
        set(p.figura, 'UserData', navigatorInfo);
      end
    end
    
    function navigatorInfo = reshowNavigatorInfo(p, navigatorInfo)
      navigatorInfo = navigateSimulation(p, navigatorInfo.rec, navigatorInfo.Ts);
    end
    
    %given that the navigator whose navigatorInfo is provided is still
    %alive, this function plays the desired timesteps
    function playNavigator(p, navigatorInfo, Ts)
      %extract info
      rec     = navigatorInfo.rec;
      edTime  = navigatorInfo.edTime;
      edIndex = navigatorInfo.edIndex;
      slider  = navigatorInfo.slider;
      mainTs  = navigatorInfo.Ts;
      %determine the indexes
      [esmiembro, locs] = ismember(Ts, mainTs);
      %for each timestep
      for k=1:numel(Ts)
        if esmiembro(k)
          %force the uicontrols
          updateControls(p, rec, mainTs, locs(k), edTime, edIndex, slider);
        end
      end
    end

    %given that the navigator whose navigatorInfo is provided is still
    %alive, this function glides through the desired timesteps
    function glideNavigator(p, navigatorInfo, hop, numhops)
      %extract info
      rec     = navigatorInfo.rec;
      edTime  = navigatorInfo.edTime;
      edIndex = navigatorInfo.edIndex;
      slider  = navigatorInfo.slider;
      mainTs  = navigatorInfo.Ts;
      idx = round(get(slider, 'Value'));
      n   = 1;
      while ((idx+hop)<=numel(mainTs)) && ((idx+hop)>=1) && (n<=numhops)
        idx = idx+hop; n = n+1;
        set(edTime,  'String', num2str(mainTs(idx)));
        set(edIndex, 'String', num2str(idx));
        set(slider,  'Value',  idx);
        playSimulation(rec, mainTs(idx), @(t, st) frameWrapper(p, t, st));
%         updateControls(p, rec, mainTs, idx, edTime, edIndex, slider);
      end
    end
    
    function setCameraProps(p, modo, props)
      if nargin<3
        props = cameraProps(get(p.figura, 'CurrentAxes'));
      end
      p.camProps = props;
      p.camMode  = modo;
    end
  end
  
end