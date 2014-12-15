function image = drawAsImage(pos, springEnds, imageLimits, imageDimensions, springColors)

% corners = imageLimits([1 3; 1 4; 2 3; 2 4]);
% corners = [corners; corners+[1 1; 1 -1; -1 1; -1 -1]];
% pos     = [pos; corners];

if ~isa(springColors, 'single')
  springColors  = single(springColors)';
else
  springColors  = springColors';
end
% springEnds      = int32(springEnds);
pos             = int32([affinTransform(pos(:,1)', [0, imageDimensions(1)-1], imageLimits([1 2])); ...
                         affinTransform(pos(:,2)', [imageDimensions(2)-1, 0], imageLimits([3 4]))]);
% imageDimensions = int32(imageDimensions);

image = drawLinesAsImageTrueColor(pos, int32(springEnds)'-1, int32(imageDimensions), springColors);