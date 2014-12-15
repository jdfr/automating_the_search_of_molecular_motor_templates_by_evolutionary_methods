function xOrig = affinTransform(xDest, rangeOrig, rangeDest)
%Range transformation. Parameters:
%xDest:  x coords to be evaluated
%rangeOrig, rangeDest: xDest values will be interpreted in the scale
%                      [rangeDest(1)...rangeDest(2)], but will be
%                      transformed to values in scale
%                      [rangeOrig(1)...rangeOrig(2)]

xOrig = ((xDest-0)-rangeDest(2))*((rangeOrig(2)-rangeOrig(1))/(rangeDest(2)-rangeDest(1))) + rangeOrig(2);



