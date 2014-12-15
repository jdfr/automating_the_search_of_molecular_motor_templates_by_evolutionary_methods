function [ss rec] = simulateSystem(ss, rec, tSteps, recordInitialState, clearOldModifications, recFreq, varargin)
try
  %defalt values for optional arguments
  if nargin<5; clearOldModifications = true; end
  
  doRec = ~isempty(rec);
  %if needed, record the given state
  if doRec
    if (nargin<6) || isempty(recFreq); recFreq = 20; end
    if recordInitialState
      rec = recordAllState(rec, ss);
    end
  end
  stepsByTick = ss.stepsByTick;
  
  step           = ss.tick/stepsByTick;
  nsteps         = (tSteps(end)-ss.t)/step;
  steps          = linspace(ss.t, tSteps(end), nsteps+1);
  step05         = step/2;
  
  chgStep        = ss.chgStep;
  if isempty(chgStep)
    chgStep = @localUG;
  end
  useBrownian   = ss.brownianFractor~=0;
  if useBrownian
    brownType   = ss.brownianFractor>0;
    brownFactor = abs(ss.brownianFractor);
  end
  postProcess   = ss.postProcess;
  pP            = ~isempty(postProcess);
  d3            = size(ss.pos,2)>2;
  checkVecs     = ss.checkVecs;
  FO            = ss.firstOrder;
  if FO
    %in first order mode, the nodes have no momentum, and acelerations and
    %velocities are one and the same
    v0          = zeros(size(ss.pos));
  end
  
  if ss.rdyn.rchg; [rdyn r] = modifyR(ss.rdyn, ss.r, ss.t); else; r=ss.r; rdyn = ss.rdyn; end; %#ok<NOSEM>
  for k=1:nsteps
    if checkVecs && ( any(isnan(ss.pos(:))) || any(isnan(ss.vel(:))) || any(abs(ss.pos(:))>1e5) || any(abs(ss.vel(:))>1e5)  )
      ss.checkFailed = true;
      return
    end
    
    [ss r rdyn veryChanged]   = chgStep(ss, r, rdyn, ss.t);
    midt  = ss.t+step05;
    nextt = steps(k+1);
    
    if FO
      if numel(v0)~=numel(ss.pos)
        v0 = zeros(size(ss.pos));
      end
      [k1v ss] = calculateAcels(ss, ss.pos, v0, r, d3);
      if pP; [nevermind k1v] = postProcess(ss, v0, k1v, r, ss.t); end
      pk1 = ss.pos+step05*k1v;

      if rdyn.rchg; [rdyn r] = modifyR(rdyn, r, midt); end
      [k2v ss] = calculateAcels(ss, pk1, v0, r, d3);
      if pP; [nevermind k2v] = postProcess(ss, v0, k2v, r, midt); end
      pk2 = pk1+step05*k2v;

      [k3v ss] = calculateAcels(ss, pk2, v0, r, d3);
      if pP; [nevermind k3v] = postProcess(ss, v0, k3v, r, midt); end
      pk3 = pk2+step*k3v;

      if rdyn.rchg; [rdyn r] = modifyR(rdyn, r, nextt); end
      [k4v ss] = calculateAcels(ss, pk3, v0, r, d3);
      if pP; [nevermind k4v] = postProcess(ss, v0, k4v, r, nextt); end

                          ss.vel       = step*(k1v+2*k2v+2*k3v+k4v)/6;
                          ss.pos       = ss.pos + ss.vel;
    else
      k1p = ss.vel;
      [k1v ss] = calculateAcels(ss, ss.pos, ss.vel, r, d3);
      if pP; [k1p k1v] = postProcess(ss, k1p, k1v, r, ss.t); end
      pk1 = ss.pos+step05*k1p;
      vk1 = ss.vel+step05*k1v;

      if rdyn.rchg; [rdyn r] = modifyR(rdyn, r, midt); end
      k2p = vk1;
      [k2v ss] = calculateAcels(ss, pk1, vk1, r, d3);
      if pP; [k2p k2v] = postProcess(ss, k2p, k2v, r, midt); end
      pk2 = pk1+step05*k2p;
      vk2 = vk1+step05*k2v;

      k3p = vk2;
      [k3v ss] = calculateAcels(ss, pk2, vk2, r, d3);
      if pP; [k3p k3v] = postProcess(ss, k3p, k3v, r, midt); end
      pk3 = pk2+step*k3p;
      vk3 = vk2+step*k3v;

      if rdyn.rchg; [rdyn r] = modifyR(rdyn, r, nextt); end
      k4p = vk3;
      [k4v ss] = calculateAcels(ss, pk3, vk3, r, d3);
      if pP; [k4p k4v] = postProcess(ss, k4p, k4v, r, nextt); end

                          ss.pos       = ss.pos    + step*(k1p+2*k2p+2*k3p+k4p)/6;
                          ss.vel       = ss.vel    + step*(k1v+2*k2v+2*k3v+k4v)/6;
    end
%     if ~isempty(r0chs); ss.r(chs)    = ss.r(chs) + step*(k1r+2*k2r+2*k3r+k4r)/6; end
%     if ~isempty(r0direct); ss.r(direct) = r(direct); end
    if useBrownian
      perturbation = (rand(size(ss.pos))-0.5)*brownFactor;%1e-3;
      rad = ss.stick.rad;
      perturbation(:,1) = perturbation(:,1).*rad;
      perturbation(:,2) = perturbation(:,2).*rad;
      if d3
      perturbation(:,3) = perturbation(:,3).*rad;
      end
      if brownType
        %random acceleration added
        ss.vel = ss.vel + perturbation;
      else
        %random velocity added
        ss.pos = ss.pos + perturbation;
      end
    end

%     if cP; [ss veryChanged] = postConstrs(ss, veryChanged); end
    
    ss.t = nextt;
    
    %recording
    if doRec %&& mod(k, stepsByTick)==0
      if veryChanged %changed && doRec
        rec = recordAllState(rec, ss);
      elseif mod(k, recFreq)==0
        rec = recordDynState(rec, ss, ss.t, makeEvolvedState(ss));
      end
    end
    
  end
  
  %prune old modifications
  if clearOldModifications
    %ss = pruneOldModificationsToDerivatives(ss, false, false);
    [rdyn r] = clearOldPerturbations(ss, rdyn, r, ss.t);
  end
  
  ss.r = r;
  ss.rdyn = rdyn;

catch ME
  if exist('r', 'var')
    ss.r = r;
    ss.rdyn = rdyn;
  end
  ss.simError = showError(ME);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rdyn r] = modifyR(rdyn, r, time)
ra = rdyn.radapt;
r(ra) = rdyn.r0(ra)+(exp((rdyn.rtb(ra)-time)./rdyn.rt(ra)).*(rdyn.rspan(ra)));


% towipe = rdyn.ltime>=time;
% if any(towipe)
%   r(rdyn.rdirect(towipe)) = rdyn.rfinal;
%   rdyn.rfinal(towipe) = [];
%   rdyn.lineA(towipe) = [];
%   rdyn.lineB(towipe) = [];
%   rdyn.ltime(towipe) = [];
%   rdyn.rdirect(towipe) = [];
%   rdyn.rchg  = ~isempty(rdyn.rdirect);
% end
% r(rdyn.rdirect) = rdyn.lineA*time+rdyn.lineB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rmod = makeRMOD(ss, time, rc, re) %#ok<DEFNU>


rdyn                = ss.rdyn;
rmod                = rdyn.rmod;

if any(rc)
rmod(rc)            = 1+rdyn.ramp(rc).*sin(rdyn.rfreq(rc)*time+rdyn.rphase(rc));
end
if any(re)
rmod(re)            = (1+exp(-(time-rdyn.rtb(re))./rdyn.rtm(re)).*(rmod(re)-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ss, rmod] = updateEM(ss, rmod, time, prox)%#ok<DEFNU> %, inContact)

rdyn                = ss.rdyn;

excitants           = (rmod<=rdyn.em.excitedRange(1)) & (rmod>=rdyn.em.excitedRange(2));
if any(excitants)
  
  %THIS IS IF POINTS WITH MANY CONNECTIONS ARE ALLOWED
  % neighs              = find(any(ss.springsNeighs(:,excitants), 2));
  % neighs              = neighs(rdyn.rexcitable(neighs));
  % neighs              = neighs(rmod(neighs)<=rdyn.em.restMax);
  % changed             = ~isempty(neighs);
  % if changed
  %   ss.rdyn.rmod(neighs)   = rdyn.em.excitedVal;
  %   ss.rdyn.rtb(neighs)    = time;
  %   rmod(neighs)           = rdyn.em.excitedVal;
  % end
  %THIS IS ONLY POINTS ATTACHED TO JUST ONE SPRING ARE ALLOWED
  relative    = ss.pos(prox(:,2),:)-ss.pos(prox(:,1),:);
  sqrdD       = realsqrt(sum(relative.*relative, 2));
  inContact = prox(sqrdD<=(ss.allr(prox(:,1)) + ss.allr(prox(:,2))),:);
  %translate inContact (pairs of points which are in contact) from points to
  %springs, given that each point is attached to just one spring
  inContact           = ss.pSpring(inContact);
  if size(inContact,2)==1
    inContact = inContact';
  end
  springsNeighs       = logical(sparse(inContact(:,1), inContact(:,2), 1, ss.nsprings, ss.nsprings));
  neighs              = find(any(springsNeighs(:,excitants), 2));
  neighs              = neighs(rdyn.rexcitable(neighs));
  neighs              = neighs(rmod(neighs)<=rdyn.em.restMax);
  changed             = ~isempty(neighs);
  if changed
    ss.rdyn.rmod(neighs)   = rdyn.em.excitedVal;
    ss.rdyn.rtb(neighs)    = time;
    rmod(neighs)           = rdyn.em.excitedVal;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drchs = calculateR(r, chs, r0chs, rmodchs, rtchs) %#ok<DEFNU>
if ~isempty(r0chs)
  drchs  = (((r0chs.*rmodchs))-r(chs))./rtchs;
else
  drchs = zeros(size(r0chs));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ss, r, rdyn, modif] = localUG(ss, r, rdyn, t) %#ok<INUSD>
modif= false;

