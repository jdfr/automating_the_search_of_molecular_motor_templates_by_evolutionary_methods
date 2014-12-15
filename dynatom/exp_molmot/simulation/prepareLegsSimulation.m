%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evparams legerror ss] = prepareLegsSimulation(ss, evparams)

if isfield(ss, 'atomSprings') && (size(ss.atomSprings, 1)>0)
  ss        = removeAtoms(ss, (1:size(ss.atomSprings, 1))', false, false);
end
ss          = fuseRepeatedSprings(ss, false);
ss.k(1:end) = evparams.walker.springK;

legerror = '';

switch evparams.walker.rewiringPlace
  case 'none'
    ss = structure2walker(ss, evparams);
    if (~isfield(evparams.walker.correct, 'rewiring')) || ~(evparams.walker.correct.rewiring)
      ss = structure2walker(ss, evparams);
    end
    if ischar(ss.legs);
      legerror = ss.legs;
    end
  case 'before'
    ss = rewireLegs(ss, evparams);
    ss = structure2walker(ss, evparams);
    if ischar(ss.legs);
      legerror = ss.legs;
    end
  case 'after'
    ss = structure2walker(ss, evparams);
    if ischar(ss.legs);
      legerror = ss.legs;
    else
      ss = rewireLegs(ss, evparams);
    end
  otherwise
    error('param evparams.walker.rewiring==%s, not understood!!!', any2str(evparams.walker.rewiringPlace));
end

if isempty(legerror)
  ss.t                  = 0;

  ss.brownianFractor    = evparams.walker.brownianFactor;

  ss.k(1:end)           = evparams.walker.springK;
  ss.r                  = realsqrt(sum(realpow(ss.pos(ss.springEnds(:,1),:)-ss.pos(ss.springEnds(:,2),:), 2), 2));

  %although strictly speaking this is a bug correction, the bug had some
  %interesting consequences: the legs came interlocked wherever they met.
  if isfield(evparams.walker.correct, 'rowRad') && evparams.walker.correct.rowRad
    ss.stick.rad(1:end) = evparams.walker.pointAllr;
  else
    ss.stick.rad(1:end) = evparams.walker.pointRad;
  end
  ss.stick.allr(1:end)  = evparams.walker.pointAllr;
  ss.row                = evparams.walker.row;
  if ~isfield(ss.row, 'attachMode')
    ss.row.attachMode   = 0;
  end
  if ~isfield(ss.row, 'lockFreeze')
    ss.row.lockFreeze   = false;%false; true;
    %ss.row.freezeK      = 0.01;
  end
  ss.row.correct.attach = isfield(evparams.walker.correct, 'attach') && evparams.walker.correct.attach;
  ss.row.locked         = false(size(ss.m));
  ss.row.lockednball    = zeros(size(ss.m));
  ss.row.activePoints   = [ss.legs(1).toes(:); ss.legs(2).toes(:)];
  ss.row.attractive     = false(size(ss.m));
  ss.row.attractive(ss.legs(1).toes) = true;
  ss.row.rad            = evparams.walker.pointRad(ones(size(ss.m)), 1);
  if ~any(cellfun(@(x)iscell(x) && strcmp(x{1}, 'row'), ss.pointVars))
    ss.pointVars(end+1:end+4)        = {{'row', 'locked'}, {'row', 'lockednball'}, {'row', 'attractive'}, {'row', 'rad'}};
    ss.pointDefaultVals(end+1:end+4) = {false,             0,                      false,                 0};
  end

  for k=1:numel(ss.legs)
    if (~isfield(ss.legs(k), 'dynLengthFactor')) || isempty(ss.legs(k).dynLengthFactor)
      ss.legs(k).dynLengthFactor = evparams.walker.leg.dynLengthFactor;
    end
    bnds          = ss.legs(k).ABCbinding;
    ABCw          = ss.legs(k).ABCweight;
    ABCb          = ss.pos(bnds,:);
    atpp          = sum(bsxfun(@times, ABCb, ABCw))/sum(ABCw);
    ss.legs(k).atpLengths = realsqrt(sum(realpow(bsxfun(@minus, ABCb, atpp), 2), 2)) .* ss.legs(k).dynLengthFactor;
    ss.legs(k).dynLengths = evparams.walker.leg.dynLengths;
    ss.legs(k).timeAtp    = evparams.walker.leg.timeAtp;
    ss.legs(k).timeRelax  = evparams.walker.leg.timeRelax;
    ss.legs(k).atpK       = evparams.walker.leg.atpK;
    ss.legs(k).atp        = [];
    ss.legs(k).atpb       = [];
    ss.legs(k).atpCount   = 0;
    ss.legs(k).relaxCount = 0;
    ss.legs(k).state      = evparams.walker.leg.stateInit{k};
    ss.legs(k).atpt       = ss.t+evparams.walker.leg.stateInitTime(k);
  end

  if isfield(evparams.walker, 'stateMachine')
    ss.chgStep         = evparams.walker.stateMachine;
  else
    ss.chgStep         = @walkerMachine;
  end
  ss.postProcess       = @collisionsWithRow;

  ss.allStateVars{end+1}  = 'row';
  ss.allStateVars{end+1}  = 'legs';

  if isfield(evparams.walker, 'd3') && strcmp(evparams.walker.d3.rotationMode, 'handoverhand')
    hoh                                     = ss.hoh;
    ss                                      = rmfield(ss, 'hoh'); %remove these data calculated in rotateLegs
    ss.updateSAP                            = doCreateSAPFun(evparams, hoh);
    ss.sap                                  = ss.updateSAP(ss);
    ss.useSAP                               = true;
  end
end

evparams.ss             = ss;
evparams.fixedTimeSpent = evparams.walker.fixedTimeSpent;
evparams.maxCPUTimeEval = evparams.walker.maxCPUTimeEval;
evparams.devEThreshold  = evparams.walker.devEThreshold;

function fun = doCreateSAPFun(evparams, hoh)
fun = @(ss)createSAPFun(ss, evparams, hoh);

function sap = createSAPFun(ss, evparams, hoh)
  %set up the sweep and prune classical parameters
  ss.row.rad(1:end)                       = evparams.walker.pointRad; 
  ss.stick.rad(1:end)                     = evparams.walker.pointAllr; 
  ss.stick.allr(1:end)                    = evparams.walker.pointAllr;
  ss.stick.penk(1:end)                    = evparams.walker.row.ballPenk; %set up the penetration penalties
  %set up the sweep and prune filter: we do not want hinge vertices 
  %to suffer interpenetration forces
  sapFilter                               = zeros(size(ss.pos,1), size(ss.pos,1), 'int8'); 
  pivot                                   = hoh.pivot;
  hinge                                   = unique(ss.springEnds(any(ss.springEnds==pivot,2),:));
  sapFilter(hinge, :)                     = -10;
  sapFilter(:,     hinge)                 = -10;
  %set up sweep and prune mode: there will be two sets (two legs), and
  %interpenetrations will happen between members of different sets
  sets                                    = zeros(size(ss.m), 'uint8');
  sets(hoh.idxoffset+1:2*hoh.idxoffset-1) = uint8(1);
  for k=1:numel(ss.legs)
    atp = ss.legs(k).atp;
    if ~isempty(atp)
      sets(atp) = k-1;
    end
  end
  mode                                    = false;
  sap                                     = createSweepAndPrune(ss.pos, ss.stick.allr, sets, mode, sapFilter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = rewireLegs(ss, evparams)

pos = ss.pos;

N = size(pos,1);

sps = ss.springEnds;


conns = sparse(sps(:), [sps(:,2); sps(:,1)], true, N,N);
conns2 = false(size(conns));

rewireModes   = evparams.walker.rewiringMode;
rewireLengths = evparams.walker.rewireLength;

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
