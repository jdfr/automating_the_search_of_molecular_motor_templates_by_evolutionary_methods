function [signo idxSpring points] = decodeBest(structs)

if (~isfield(structs, 'best')) || ( isempty(structs.best) || (structs.best==0) )
  signo = 'none';
  idxSpring = 0;
  points = [];
  return
end

if structs.best>size(structs.selecteds,1)
  signo = 'stressed';
  best = structs.best-size(structs.selecteds,1);
else
  signo = 'compressed';
  best = structs.best;
end

s = sort(structs.ss.springEnds, 2);

points = structs.selecteds(best,:);

idxSpring = find( (s(:,1)==min(points)) & (s(:,2)==max(points)) );
