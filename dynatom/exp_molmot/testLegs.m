function [rec ss] = testLegs(ss, evparams)
walkerParams.minAmountOfToes = @(conns, vs) min(8, ceil(numel(vs)/2));
walkerParams.cutoff = 1.5;
walkerParams.rewiringMode = {'geotopological', 'georecurtopo'};'geotopological';'topological';'geotopological';'georecurtopo';
walkerParams.rewireLength = [40 10];(1.5*10*1.3)*1.5;
walkerParams.initialDepth = 1;
walkerParams.ABCMode = 'heurconcave'; 'meanpos'; 'heurconcave';
walkerParams.tooSmallNeighbour = 2/3;
ss.legs = structure2walker(ss, walkerParams);
ss = rewireLegs(ss, walkerParams);

% testANM(ss, true);
% testANM(ss2, true);


ss.t                 = 0;

ss.brownianFractor   = 1;1;0.1;

ss.k(1:end)          = 100;

ss.stick.rad(1:end)  = 3;
ss.stick.allr(1:end) = 4;
ss.row.ballAllr      = ss.stick.allr(1);
ss.row.ballSep       = 8;12;
ss.row.ballPenk      = 500;ss.stick.penk(1);
ss.row.typeSpec      = [5 5];
ss.row.ballRad       = 3;
ss.row.active        = true;
ss.row.allPoints     = true;false;
ss.row.locked        = false(size(ss.m));
ss.row.lockednball   = zeros(size(ss.m));

ss.legs(1).state     = 'sticky';
ss.legs(1).atpt      = nan;
ss.legs(2).state     = 'relaxation';
ss.legs(2).atpt      = ss.t;
ss.row.activePoints  = [ss.legs(1).toes; ss.legs(2).toes];
ss.row.attractive    = false(size(ss.m));
ss.row.attractive(ss.legs(1).toes) = true;

for k=1:numel(ss.legs)
  ss.legs(k).dynLengthFactor = 0;0.1;0.5;0.1;
  bnds          = ss.legs(k).ABCbinding;
  ABCw          = ss.legs(k).ABCweight;
  ABCb          = ss.pos(bnds,:);
  atpp          = sum([ABCb(:,1).*ABCw, ABCb(:,2).*ABCw])/sum(ABCw);
  ss.legs(k).atpLengths = realsqrt(sum(realpow([ABCb(:,1)-atpp(1), ABCb(:,2)-atpp(2)], 2), 2)) * ss.legs(k).dynLengthFactor;
  ss.legs(k).dynLengths = false;
  ss.legs(k).timeAtp = 50;
  ss.legs(k).timeRelax = 50;
  ss.legs(k).atpK = 200;
  ss.legs(k).atp  = [];
  ss.legs(k).atpb = [];
  ss.legs(k).atpCount = 0;
  ss.legs(k).relaxCount = 0;
end

ss.chgStep           = @walkerMachine;
ss.postProcess       = @collisionsWithRow;

ss.allStateVars{end+1}  = 'row';
ss.allStateVars{end+1}  = 'legs';

ss = rotateLegs(ss, walkerParams);

evparams.ss             = ss;
evparams.fixedTimeSpent = 3000;
evparams.maxCPUTimeEval = 300;
evparams.devEThreshold  = -inf;

useKeyboard = true;
rndSeed     = [];
genome      = [];

[ss rec] = DADevelopGenome(genome, evparams, rndSeed, useKeyboard, []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = rewireLegs(ss, params)

pos = ss.pos;

N = size(pos,1);

sps = ss.springEnds;


conns = sparse(sps(:), [sps(:,2); sps(:,1)], true, N,N);
conns2 = false(size(conns));

rewireModes   = params.rewiringMode;
rewireLengths = params.rewireLength;

if ~iscell(rewireModes)
  rewireModes = {rewireModes};
end

for q=1:numel(rewireModes)
  c2    = full(conns);
  switch rewireModes{q}
    case 'topological'
      for k=1:N
        nghs       = conns(:,k)';
        neighs     = find(nghs);
        for z=1:numel(neighs)
          nz       = neighs(z);
          c2(nz,:) = c2(nz,:) | nghs;
        end
      end
    case 'geotopological'
      distsOK = distanceMatrix(pos)<=rewireLengths(q);
      distsOK(1:(N+1):(N*N)) = false;
      for k=1:N
        nghs         = conns(:,k)';
        neighs       = find(nghs);
        for z=1:numel(neighs)
          nz         = neighs(z);
          c2(nz,:)   = c2(nz,:) | (nghs & distsOK(nz,:));
        end
      end
    case 'georecurtopo' %georecursivetopological
      distsOK = distanceMatrix(pos)<=rewireLengths(q);
      distsOK(1:(N+1):(N*N)) = false;
      tempconns        = conns;
      changed = true;
      while changed
        changed        = false;
        for k=1:N
          nghs         = tempconns(:,k)';
          neighs       = find(nghs);
          for z=1:numel(neighs)
            nz         = neighs(z);
            newneighs  = (nghs & distsOK(nz,:));
            changed    = changed || any(newneighs & (~c2(nz,:)));
            c2(nz,:)   = c2(nz,:) | newneighs;
          end
        end
        tempconns = c2;
      end
    otherwise
      error('evparams.walker.rewiring %s not understood!!!', any2str(rewireModes{q}));
  end
  conns2 = conns2 | c2;
end

conns2 = tril(conns2, -1);
[a b]  = find(conns2);
sps    = [a b];
k      = ss.k(1);
r      = realsqrt(sum(realpow(pos(sps(:,1),:)-pos(sps(:,2),:), 2), 2));
ss     = removeSprings(ss, (1:numel(ss.k))');
ss     = addSprings(ss, sps, k, r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = rotateLegs(ss, params)

%center of mass of toes of leg 1
t1 = ss.legs(1).toes;
w1 = ss.m(t1)/sum(ss.m(t1));
cm1 = [sum(ss.pos(t1,1).*w1) sum(ss.pos(t1,2).*w1)];

%center of mass of toes of leg 2
t2 = ss.legs(2).toes;
w2 = ss.m(t2)/sum(ss.m(t2));
cm2 = [sum(ss.pos(t2,1).*w2) sum(ss.pos(t2,2).*w2)];

ss.pos = bsxfun(@minus, ss.pos, cm1);   %displace cm1 to the origin
t12 = cm2-cm1;                          %get the vector from cm1 to cm2
t12d = realsqrt(sum(realpow(t12, 2)));  %get distance
t12u = t12/t12d;                        %get unit vector

%rotate the structure so the line from cm1 to cm2 is horizontal
%(t12u=[cos(alpha), sin(alpha)], and we want to rotate the structure by
%the amount -alpha)
ss.pos = ss.pos*[t12u(1), -t12u(2); t12u(2), t12u(1)];

%if most of the structure is below, flip it
if sum(ss.pos(:,2))<0
  ss.pos(:,2) = -ss.pos(:,2);
end

%get the lowest toe in each leg
[nevermind mt1] = min(ss.pos(t1,2));
mt1 = ss.pos(t1(mt1),:);
[nevermind mt2] = min(ss.pos(t2,2));
mt2 = ss.pos(t2(mt2),:);

ss.pos = bsxfun(@minus, ss.pos, mt1);   %displace mt1 to the origin
t12 = mt2-mt1;                          %get the vector from mt1 to mt2
t12d = realsqrt(sum(realpow(t12, 2)));  %get distance
t12u = t12/t12d;                        %get unit vector

%rotate the structure so the line from mt1 to mt2 is horizontal
%(t12u=[cos(alpha), sin(alpha)], and we want to rotate the structure by
%the amount -alpha). After this, both toes are at line y=0
ss.pos = ss.pos*[t12u(1), -t12u(2); t12u(2), t12u(1)];

%now, offset the structure so the lowest toes are just in contact with the
%row
ss.pos(:,2) = ss.pos(:,2)+ss.row.ballAllr+ss.stick.allr(1)-params.initialDepth;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function axiswindow = getAxisWindow(pos, margen)
mins = min(pos);
maxs = max(pos);
intervs = maxs-mins;
meds = (mins+maxs)/2;
maxi = max(intervs)*(1+margen)/2;
axiswindow = [meds(1)+[-1 1]*maxi, meds(2)+[-1 1]*maxi];
