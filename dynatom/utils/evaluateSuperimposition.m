function [counts ndatoms] = evaluateSuperimposition(evparams, ss)

pos = ss.pos;

imageDimensions = repmat(evparams.imageDimensions, 1, 2);
imageLimits     = [min(pos(:,1)) max(pos(:,1)) min(pos(:,2)) max(pos(:,2))];
lx = imageLimits(2)-imageLimits(1);
ly = imageLimits(4)-imageLimits(3);
toadd = abs(lx-ly)/2;
if lx<ly
  idxs = [1 2];
else
  idxs = [3 4];
end
imageLimits(idxs) = imageLimits(idxs) + [-toadd, +toadd];

change = @(x) [affinTransform(x(:,1)', [1, imageDimensions(1)], imageLimits([1 2])); ...
               affinTransform(x(:,2)', [imageDimensions(2), 1], imageLimits([3 4]))];
% change = @(x) [affinTransform(x(:,1)', [0, imageDimensions(1)-1], imageLimits([1 2])); ...
%                affinTransform(x(:,2)', [imageDimensions(2)-1, 0], imageLimits([3 4]))];

atoms   = int32(unique(sort(ss.atomPoints, 2), 'rows'));
ndatoms = size(atoms, 1);

chulls  = cell(ndatoms,1);
centers = zeros(ndatoms,2);
for k=1:ndatoms
  ap = atoms(k,:);
  chulls{k}    = int32(ap(ConvHull2D(pos(ap,1),pos(ap,2))));
  centers(k,:) = mean(pos(chulls{k}(1:end-1),:));
end

pos     = int32(change(pos));
centers = round(change(centers))';

image = zeros(imageDimensions([2 1]), 'uint8');

imageDimensions = int32(imageDimensions);

for k=1:ndatoms
  img = drawLinesAsBoolImage(pos, chulls{k}([(1:end-1); (2:end)])-1, imageDimensions);
  img2 = imfill(img, centers(k,[2 1]));
  %this check is necessary for degenerate cases: if an atom is just a line,
  %imfill will flood all the image. To prevent incorrect calculations, in
  %these cases we just add the lines, not the flooded img
  if all(img2(:))
    image = image + uint8(img);
  else
    image = image + uint8(img2);
  end
%   imshow(img, [0 0 0; 1 1 1; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0]);
%   pause;
%   image = image + uint8(img);
%   imshow(image, [0 0 0; 1 1 1; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0]);
%   pause;
end

counts = hist(single(image(:)), single(0:max(image(:))));

