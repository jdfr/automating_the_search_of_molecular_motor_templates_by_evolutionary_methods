function ss = getProtein(fname)

f = fopen(fname, 'r');

lines = textscan(f, '%[^\n]');

if (numel(lines)==1) && iscell(lines{1})
  lines = lines{1};
end

lines = lines(cellfun(@(x)(numel(x)>15)  && all(x([1 2 3 4 14 15 22])=='ATOMCAA'), lines));
pos   = cellfun(@(x)str2num(['[' x(32:55) ']']), lines, 'uniformoutput', false); %#ok<ST2NM>
pos   = vertcat(pos{:});
idxs  = cellfun(@(x)str2double(x(24:26)), lines);

L0    = 10;

[sp1 sp2] = find(triu(distanceMatrixSquared(pos)<(L0*L0), 1));

r = realsqrt(sum(realpow(pos(sp1,:)-pos(sp2,:), 2), 2));

ss = BasicSystem(3, 0, 1);

ss = addPoints(ss, pos);

ss = addSprings(ss, [sp1 sp2], 1, r);

ss.idxs = idxs;
