function y = sigmoidP(xDest, slope, phaseX, rangeOrig, rangeDest, shiftX, amplitudeY, shiftY)
%parameterized sigmoid
%xDest:  x coords to be evaluated
%slope:  sigmoid's slope
%phaseX: logarithmic parameter controlling the phase (similar to shiftX),
%        (default to 1)
%rangeOrig, rangeDest: xDest values will be interpreted in the scale
%                      [rangeDest(1)...rangeDest(2)], but will be
%                      transformed to values in scale
%                      [rangeOrig(1)...rangeOrig(2)]
%shiftX:     in rangeDest coords, shift to be applied to x coords 
%            (default to 0)
%amplitudeY, shiftY:   y values will go be in the range
%                      [-shiftY..amplitudeY-shiftY], default values to
%                      amplitudeY=1, shiftY=0

if nargin<2
  slope=1;
end
if nargin<3
  phaseX=1;
end
if nargin<4 %slope
  rangeOrig=[-10 10];
end
if nargin<5 %slope
  rangeDest=[-10 10];
end
if nargin<6 %slope
  shiftX = 0;
end
if nargin<7 %slope
  amplitudeY=1;
end
if nargin<8 %slope
  shiftY=0;
end

xOrig = ((xDest-shiftX)-rangeDest(2))*((rangeOrig(2)-rangeOrig(1))/(rangeDest(2)-rangeDest(1))) + rangeOrig(2);

y = amplitudeY./(1+phaseX*exp(-slope*xOrig)) - shiftY;

