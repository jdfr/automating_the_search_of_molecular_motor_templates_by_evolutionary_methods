%points: 1x4: indexes of atom's points. Zero means new point
%pos: for new points, their positions
%springsR, springsK: parameters for atom's springs
function [ss atomIndex] = addAtomComplex(ss, points, pos, springsK, springsR, varargin) %varargin={atomGenome, atomState}

atomIndex = size(ss.atomPoints,1)+1;

fixedPoints = points(points~=0);
newps = find(points==0);
if numel(newps)~=size(pos,1)
  error('the parameter pos represents the position of new points, ergo its row number must match with the parameter point''s amount of zeroes');
end

if numel(newps)>0
  [ss, indexesNewPoints] = addPoints(ss, pos);
  points(newps) = indexesNewPoints;
else
  indexesNewPoints = [];
end

%update/initialize the pointAtoms entry for each point
ss.pointAtoms(fixedPoints)      = cellfun(@(atoms) [atoms; atomIndex], ss.pointAtoms(fixedPoints),      'UniformOutput', false);
for k=1:numel(indexesNewPoints)
  ss.pointAtoms{indexesNewPoints(k)} = atomIndex;
end
%ss.pointAtoms(indexesNewPoints) = repmat({atomIndex}, cellfun(@(x)      atomIndex,         ss.pointAtoms(indexesNewPoints), 'UniformOutput', false);

localSprings = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

[ss springs] = addSprings(ss, points(localSprings), springsK, springsR);
ss.springAtom(springs) = atomIndex;

ss = auxAssignment(ss, 1, 1, {points, springs', varargin{:}}, ss.atomDefaultVals, ss.atomVars, 'atom');
