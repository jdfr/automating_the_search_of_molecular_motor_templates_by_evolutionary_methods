function [relaxations deformations] = tryPerturbations(ss, seedRnd, pertParams, points)

it = cputime;
ic = clock;

if ~exist('points', 'var')
  points   = get3RefPoints(ss);
end

ntimes = pertParams.ntimes;

nb = 0;

for k=ntimes:-1:1
  if nb>0; fprintf(repmat('\b', 1, nb)); end;
  nb = fprintf('Doing try %05d; seconds in total: %s; cpu seconds in total: %s', k, mat2str(etime(clock, ic)), mat2str(cputime-it));
  [relaxations(k) deformations(k)] = tryPerturbation(ss, seedRnd, pertParams, points); %#ok<AGROW>
end

fprintf('\n Showing it...');

figure; hold on; grid on;
mrk=@(c,s){'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'Marker', 'o', 'MarkerSize', s};
m1=mrk('b',6);
m2=mrk('g',4);
for k=1:ntimes
  line(relaxations(k).d12, relaxations(k).d13, relaxations(k).d23, 'Color', 'r');
end
for k=1:ntimes
  line(relaxations(k).d12(end), relaxations(k).d13(end), relaxations(k).d23(end), m2{:});
end
line(0,0,0, m1{:});
view(3);
axis equal;


%figure; hold on; grid on; mrk=@(c,s){'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'Marker', 'o', 'MarkerSize', s}; m1=mrk('b',6); m2=mrk('g',4); for k=1:ntimes;   line(zrelax(k).d12, zrelax(k).d13, zrelax(k).d23, 'Color', 'r'); end; for k=1:ntimes;   line(zrelax(k).d12(end), zrelax(k).d13(end), zrelax(k).d23(end), m2{:}); end; line(0,0,0, m1{:}); axis equal;

