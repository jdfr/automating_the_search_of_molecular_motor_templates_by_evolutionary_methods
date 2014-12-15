function showCAlpha(ss, mut)%, mode)

if nargin>1
  doIt(ss, @(ss)showMut(ss, mut));
else

%   doIt(ss, @showMol);

  if ss.legs(1).atp>ss.legs(2).atp
    a = 1;
    b = 2;
  else
    a = 2;
    b = 1;
  end

  ss = removeSprings(ss, ss.legs(a).atpb(:));
  ss = removeSprings(ss, ss.legs(b).atpb(:));
  ss =  removePoints(ss, ss.legs(a).atp);
  ss =  removePoints(ss, ss.legs(b).atp);
  ss = rmfield(ss, 'legs');

%   doIt(ss, @showMol);
  doIt(ss, @showThirdEigenVector);
end

function doIt(ss, fun)

figure;

hold on;

fun(ss);

%light('Position',[0 1000 0],'Style','infinite');

axis equal;
view(180,0);
set(gca, 'CameraUpVector',  get(gca, 'CameraUpVector')  * (-1) );
%set(gca, 'CameraViewAngle', get(gca, 'CameraViewAngle') * 0.4  );
set(gcf, 'Color', [1 1 1]);
axis off;
cameratoolbar(gcf, 'Show');


function showMut(ss, mut)
[vertices, links] = showMol(ss);

c=0.6;
c=[c c c];

for k=1:numel(vertices)
  set(vertices{k}, 'FaceColor', c);
end
for k=1:numel(links)
  set(links{k}, 'FaceColor', c);
end

mut = str2num(mut);
mut = [sort(mut(:,1:2), 2), mut(:,3)>=1, mut(:,3)];
s = sort(ss.springEnds, 2);

    [zX zY] = meshgrid(0:0.1:10*2*pi); 
    
    img1=((sin(zX/3)+1)/2);
    img1 = repmat(img1, [1, 1, 3]);
    img2=((sin(zY)+1)/2);
    img2 = repmat(img2, [1, 1, 3]);

for z=1:size(mut,1)
  km = z;
  ln = find((s(:,1)==mut(km,1)) & (s(:,2)==mut(km,2)) );
  if mut(km,3)
    img = img1;
  else
    img = img2;
  end
  set(links{ln}, 'FaceColor', 'texturemap', 'CData', img);
  %set(links{ln}, 'FaceColor', 0.3*ones(1,3));
end


function [vertices, links] = showMol(ss)
% switch mode
%   case 'myosin'
    ns = 40;
    rs = 3.8/2;
    cs = 0.3;

    nc = 10;
    rc = 0.3;
    cc = 0.4;
%   case 'struct1'
%     ns = 20;
%     rs = 3.8/2;
%     cs = 0.3;
% 
%     nc = 10;
%     rc = 0.3;
%     cc = 0.4;
%   otherwise
%     error('jarllll!!!');
% end


md = false;true;false;

[X,Y,Z] = sphere(ns);

if md
  X=-X;
  Y=-Y;
  Z=-Z;
end

[F,V]   = surf2patch(X,Y,Z);

st = struct('vertices', [], 'faces', F);

args = {...
  'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
  'SpecularStrength', 0, 'DiffuseStrength', 1, ...
  'AmbientStrength', 1, 
  };

vertices = cell(size(ss.pos,1),1);

for k=1:size(ss.pos,1)
  st.vertices = bsxfun(@plus, V*rs, ss.pos(k,:));
  s = patch(st);
  set(s, 'FaceColor', cs*ones(3,1), args{:});
  vertices{k} = s;
end

links = cell(size(ss.springEnds,1),1);

for k=1:size(ss.springEnds)
  es = ss.springEnds(k,:);
  [X,Y,Z]=cylinder2P(rc,nc,ss.pos(es(1),:), ss.pos(es(2),:));
  c = surf(X,Y,Z);
  set(c, 'FaceColor', cc*ones(3,1), args{:});
  links{k} = c;
end

if isfield(ss, 'legs') && ~isempty([ss.legs.atp])
  toe  = cs*2;
  atpV = 0;
  atpE = 0;
  for k=1:numel(ss.legs)
    for z=1:numel(ss.legs(k).toes)
      t = ss.legs(k).toes(z);
      set(vertices{t}, 'FaceColor', toe*ones(3,1));
    end
    set(vertices{ss.legs(k).atp}, 'FaceColor', atpV*ones(3,1));
    for z=1:numel(ss.legs(k).atpb)
      set(links{ss.legs(k).atpb(z)}, 'FaceColor', atpE*ones(3,1));
    end
  end
end

function showThirdEigenVector(ss)

    ns  = 40;
    rs  = 3.8/2;
    cs1 = 0.3;
    cs2 = 0.6;

    nc = 10;
    rc = 2;0.05;1;
    cc = 0.4;

%numver of vertices
N = size(ss.pos, 1);
%connectivity matrix
conns = sparse(ss.springEnds(:), reshape(ss.springEnds(:,[2,1]), [], 1), true, N, N);
%kirchoff/laplacian matrix
kirchoff = sparse(1:N, 1:N, sum(conns), N,N) - conns;
[vecs vals] = eig(full(kirchoff));
[idx3 idx3] = sort(diag(vals));
vec3 = vecs(:,idx3(3));
vec3 = vec3/max(abs(vec3))*1.5;

[X,Y,Z] = sphere(ns);
md = false;
if md
  X=-X;
  Y=-Y;
  Z=-Z;
end

[F,V]   = surf2patch(X,Y,Z);

st = struct('vertices', [], 'faces', F);

args = {...
  'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
  'SpecularStrength', 0, 'DiffuseStrength', 1, ...
  'AmbientStrength', 1, 
  };

vertices = cell(size(ss.pos,1),1);

for k=1:size(ss.pos,1)
  st.vertices = bsxfun(@plus, V*rs*abs(vec3(k)), ss.pos(k,:));
  if sign(vec3(k))>=0
    csU = cs1;
  else
    csU = cs2;
  end
  s = patch(st);
  set(s, 'FaceColor', csU*ones(3,1), args{:});
  vertices{k} = s;
end

links = cell(size(ss.springEnds,1),1);

for k=1:size(ss.springEnds)
  es  = ss.springEnds(k,:);
  ps1 = ss.pos(es(1),:);
  ps2 = ss.pos(es(2),:);
  ps  = [ps1; ps2];
  c = line(ps(:,1), ps(:,2), ps(:,3), 'Color', cc*ones(3,1), 'LineWidth', rc);
%   [X,Y,Z]=cylinder2P(rc,nc,ps1, ps2);
%   c = surf(X,Y,Z);
%   set(c, 'FaceColor', cc*ones(3,1), args{:});
  links{k} = c;
end
