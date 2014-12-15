function p = SpringPlotter(varargin)

fields = {...
    ... % BasicPlotter public
    'figura';...      %figure's handle
    'figurePos' ;...   %figure's position in normalized units
    'AVIfilename';... %name of AVI file to record
    'createWindow';...%this should be always true, unless you really know what you're doing 
    'useOpenGL';...   %whether to use openGl. Enfatically recommended
    'visible';...     %info flag for figure's visibility. It is automatically set (AVI=not visible) 
    'fileAVI';...     %AVI file's handle
    'codecAVI';...    %AVI file's codec
    'fpsAVI';...      %AVI file's fps
    'fixedAxisFrame';... %whether the axes lengths are fixed. It is implied by 'axisWindow' 
    'axisWindow';...  %if []=>not fixed frames. if num=>squared axes'semiradius. if vector=>axes
    'hasGrid';...     %whether to use a grid
    'hasAxis';...     %whether to display coordinate axes
    'axisEqual';...   %whether to use 'axis equal'
    'axisSquare';...  %whether to use 'axis square'
    'hasTitle';...    %whether to put a title
    'pauseTime';...   %after each frame: []=>do pause. 0=>no pause. n=>pause n seconds 
    'viewPoint';...   %[]=> standard viewpoint. T=>transformation matrix for 'view' 
    ... % BasicPlotter protected
    'justStarted';... %flag to tell if the frame is the first to be drawn
    ... % SpringPlotter public
    'showLineInfo';...            %whether to show springs and points info
    'selectedSpringByDefault';... %if showLineInfo==true, this var tells what spring is selected by default
    'acelFactor';...              %if 0, does nothing. Otherwise, shows acelerations for points, scaled by a factor 'acelFactor'
    'showSpringAcels';...         %if acelFactor>0, shows spring acels (compressed and stretched springs) 
    'showOtherAcels';...          %if acelFactor>0, shows other acels (drag force, perturbations, enforcements and customized ones) 
    'showTotalAcels';...          %if acelFactor>0, shows total acels (sum of the previous ones) 
    'velFactor';...               %if 0, does nothing, Otherwise, shows scaled vels for points
    'shownSprings';...            %indexes of drawn springs. If empty, show all of them
    'shownPoints';...             %indexes of drawn points.  If empty, show all of them
    'colorsCodeForce';...         %if true, colors code force. If false, colors code relative length and spring stiffness.
    'maxForce';...                %when colorsCodeForce=true: the range of spring forces to be mapped in [blue...green...red] will be [maxForce, -maxForce]
    'springMaxRel';...            %when colorsCodeForce=false: the range of spring stresses to be mapped in a map [blue...green...red] will be [-springMaxRel, springMaxRel]
    'rangeK';...                  %when colorsCodeForce=false: rangeK=[min max] is a range for springs' K values: springs with K=max or higher will be displayed fully saturated; springs with K=min or lower will be displayed minimally saturated
    'minSat';...                  %when colorsCodeForce=false: between 0 and 1: minimum saturation, for springs with K=rangeK(1) (the maximum saturation is 1)
    'springMaxForce';...
    ... % SpringPlotter protected
    'textLineInfo';... %text uicontrol to show info
    'end1';... %selected spring's point 1
    'end2';... %selected spring's point 2
    'buttonDownFcn';... %callback for selected springs
};

aux = reshape([fields repmat({[]}, size(fields))]', [], 1);

p = struct(aux{:});

p.figurePos = [0.05 0.05 0.8 0.8];

p.AVIfilename = '';
p.createWindow = true;
p.useOpenGL = true;
p.visible   = true;

p.codecAVI = 'CinePak';
p.fpsAVI   = 15;

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

p.showLineInfo = true;
p.acelFactor = 0;
p.showSpringAcels = true;
p.showOtherAcels  = false;
p.showTotalAcels  = false;
p.velFactor       = 0;
p.selectedSpringByDefault = 1;
p.shownSprings = [];
p.shownPoints  = [];
p.colorsCodeForce = true;
p.maxForce     = 0.5;
p.springMaxRel = 0.03;
p.rangeK       = [0.5 1];
p.minSat       = 1/3;

p = setSPFields(p, varargin{:});

p.buttonDownFcn = @(obj, event) linePressed(p, obj);
