%subclasses must redefine this
function drawFrame(p, t, st)
  if p.hasTitle
    title(sprintf('time: %-10g; u:%-10g, e:%-10g', t, st.u, calculateEnergy(st)));
  end
  if (p.acelFactor~=0)
    [springDisplacements, compressedSpringAcels, stretechedSpringAcels, otherAcels, dynamicVel, changedDynamicVel] = getKinematicVectors(t, st);
    if p.showSpringAcels || p.showTotalAcels
      compressedSpringAcels = compressedSpringAcels*p.acelFactor;
      stretechedSpringAcels = stretechedSpringAcels*p.acelFactor;
    end
    if p.showOtherAcels || p.showTotalAcels; otherAcels = otherAcels*p.acelFactor; end
    if p.velFactor~=0; dynamicVel = dynamicVel*p.velFactor; end
    plotSprings(p, t, st, springDisplacements);
    plotPointKinematics(p, t, st, compressedSpringAcels, stretechedSpringAcels, otherAcels, dynamicVel);
  elseif (p.velFactor~=0)
    plotSprings(p, t, st);
    dynamicVel = st.vel(st.dynamical_p,:);
    plotPointKinematics(p, t, st, [], [], [], dynamicVel);
  else
    plotSprings(p, t, st);
  end
end

    function plotPointKinematics(p, t, st, compressedSpringAcels, stretechedSpringAcels, otherAcels, dynamicVel)
      showPoints      = p.shownPoints;
      showSpringAcels = p.showSpringAcels && (p.acelFactor~=0);
      showOtherAcels  = p.showOtherAcels && (p.acelFactor~=0);
      showTotalAcels  = p.showTotalAcels && (p.acelFactor~=0);
      velFactor       = p.velFactor;
      dynPoints  = find(st.dynamical_p);
      dynpos     = st.pos(st.dynamical_p,:);
      if ~isempty(showPoints) %filter out shown points which are not dynamical
        showPoints(~st.dynamical_p(showPoints)) = [];
        dynPointsNotShown = ~ismember(dynPoints, showPoints);
        %dynPoints(dynPointsNotShown) = [];
        dynpos(dynPointsNotShown,:)  = [];
        if showSpringAcels || showTotalAcels
          compressedSpringAcels(dynPointsNotShown,:) = [];
          stretechedSpringAcels(dynPointsNotShown,:) = [];
        end
        if showOtherAcels || showTotalAcels; otherAcels(dynPointsNotShown,:) = []; end
        if velFactor~=0; dynamicVel(dynPointsNotShown,:) = []; end
      end
      if showSpringAcels;
        endcsa = dynpos+compressedSpringAcels;
        endssa = dynpos+stretechedSpringAcels;
      end
      if showOtherAcels; endoa = dynpos+otherAcels; end
      if showTotalAcels; enda  = dynpos+compressedSpringAcels+stretechedSpringAcels+otherAcels; end
      if velFactor~=0;   endv  = dynpos+dynamicVel; end
      for k=1:size(dynpos,1)
        pos = dynpos(k,:);
        if showSpringAcels
          csa = [pos;endcsa(k,:)]; CSA = {'Color', 'r', 'LineWidth', 2};
          ssa = [pos;endssa(k,:)]; SSA = {'Color', 'b', 'LineWidth', 2};
        end
        if showOtherAcels; oa  = [pos;endoa(k,:)];  OA  = {'Color', 'c', 'LineWidth', 2};  end
        if showTotalAcels; ta   = [pos;enda(k,:)];  TA   = {'Color', 'k', 'LineWidth', 2}; end
        if velFactor~=0;   tv   = [pos;endv(k,:)];  TV   = {'Color', 'm', 'LineWidth', 2}; end
        switch st.ndims
          case 2
            if showSpringAcels
              line(csa(:,1), csa(:,2), CSA{:});
              line(ssa(:,1), ssa(:,2), SSA{:});
            end
            if showOtherAcels; line(oa(:,1),  oa(:,2),  OA{:}); end
            if showTotalAcels; line(ta(:,1),  ta(:,2),  TA{:});  end
            if velFactor~=0;   line(tv(:,1),  tv(:,2),  TV{:});  end
          case 3
            if showSpringAcels
              line(csa(:,1), csa(:,2), csa(:,3), CSA{:});
              line(ssa(:,1), ssa(:,2), ssa(:,3), SSA{:});
            end
            if showOtherAcels; line(oa(:,1),  oa(:,2),  oa(:,3),  OA{:}); end
            if showTotalAcels; line(ta(:,1),  ta(:,2),  ta(:,3),  TA{:}); end
            if velFactor~=0;   line(tv(:,1),  tv(:,2),  tv(:,3),  TV{:}); end
        end
      end
    end

    %plots springs
    function plotSprings(p, t, st, springDisplacements)
      showSprings = p.shownSprings;
      mxrel = p.springMaxRel;
      ranK = p.rangeK;
      mnSat = p.minSat;
      mxForce = p.maxForce;
      if isempty(showSprings)
        showSprings = 1:size(st.springEnds,1);
      end
      if ~exist('springDisplacements', 'var')
        springVectors = st.pos(st.springEnds(showSprings,2),:)-st.pos(st.springEnds(showSprings,1),:);
        springLengths = realsqrt(sum(springVectors.*springVectors,2));
        springDisplacements = st.r(showSprings)-springLengths;
      else
        springDisplacements = springDisplacements(showSprings,:);
      end
      if p.colorsCodeForce
        pseudoHookeanForces = st.k(showSprings).*springDisplacements;
        pseudoHookeanForces(pseudoHookeanForces<-mxForce) = -mxForce;
        pseudoHookeanForces(pseudoHookeanForces>mxForce)  = mxForce;
        mappingHue = (1-((pseudoHookeanForces-(-mxForce))/(mxForce-(-mxForce))));
        colors = hsv2rgb(mappingHue*(2/3), ones(size(pseudoHookeanForces)), ones(size(pseudoHookeanForces)));
      else
        springRelDisplacements = springDisplacements./st.r(showSprings);
        springRelDisplacements(springRelDisplacements<-mxrel)=-mxrel;
        springRelDisplacements(springRelDisplacements>mxrel)=mxrel;
        mappingHue = (1-((springRelDisplacements-(-mxrel))/(mxrel-(-mxrel))));
        mappingSat = (mnSat-1)*(st.k(showSprings)-ranK(2))/(ranK(1)-ranK(2))+1;
        mappingSat(mappingSat<mnSat) = mnSat;
        mappingSat(mappingSat>1) = 1;
        colors = hsv2rgb(mappingHue*(2/3), mappingSat, ones(size(springRelDisplacements)));
      end
      colors(colors<0) = 0;
      colors(colors>1) = 1;
      %if we must show info, reset points
      if p.showLineInfo
        p.end1 = [];
        p.end2 = [];
      end
      %draw springs
      for i=1:numel(showSprings)
          indexes = st.springEnds(showSprings(i),:);
          if p.showLineInfo
            infoData = {'UserData', {showSprings(i), st}, 'ButtonDownFcn', p.buttonDownFcn};
          else
            infoData = {};
          end
          %write line. This is valid for 2 and 3 dimensions, but
          %not for any other amount of them, but it is easily generalizable
          %if we need to do so
          switch st.ndims
            case 2
              h = line(st.pos(indexes,1),st.pos(indexes,2), 'Color',colors(i,:), infoData{:});
            case 3
              h = line(st.pos(indexes,1),st.pos(indexes,2),st.pos(indexes,3), 'Color',colors(i,:), infoData{:});
          end
          
          %if we must show info and this is the selected spring by default,
          %let's select it
          if p.showLineInfo && ~isempty(p.selectedSpringByDefault) && p.selectedSpringByDefault==showSprings(i)
            p.buttonDownFcn(h,[]);
          end
      end
    end

    %decompose system's dynamics
    function [springDisplacements, compressedSpringAcels, stretechedSpringAcels, otherAcels, dynamicVel, changedDynamicVel] = getKinematicVectors(t, st)
      st = prepareSimulationStructure(st);
      dynamicVel = st.vel(st.dynamical_p,:);
      masses = st.m(st.dynamical_p);
      [compressiveForces, stretechedForces, springDisplacements] = calculateDetailedSpringForces(st);
      compressedSpringAcels = bsxfun(@rdivide, compressiveForces, masses);
      stretechedSpringAcels = bsxfun(@rdivide, stretechedForces, masses);
      dragAcels = bsxfun(@rdivide, -dynamicVel*st.u, masses);
      [changedDynamicVel AfterAcel] = nonNaturalDynamics(st, t, dynamicVel, compressedSpringAcels+stretechedSpringAcels+dragAcels, [], [], [], [], []);
      otherAcels = AfterAcel-compressedSpringAcels-stretechedSpringAcels;
    end

    %this is a copy of calculatedSpringForces, with a twist: calculate them
    %separatedly for compressed and stretched springs
    function [compressiveForces, stretechedForces, springDisplacements] = calculateDetailedSpringForces(ss)
      %first, deal some local variables to cut access time to properties
      springEnds    = ss.springEnds;
      pos           = ss.pos;
      vel           = ss.vel;
      k             = ss.k;
      r             = ss.r;
      c             = ss.c;
      springsMatrix = ss.springsMatrix;
      %spring vectors
      springVectors = pos(springEnds(:,2),:)-pos(springEnds(:,1),:);
      %lengths of spring vectors
      springLengths = realsqrt(sum(springVectors.*springVectors, 2));
      %springs' displacements from rest length, signed
      springDisplacements = r-springLengths;
      %hookean forces exerted by springs. They are for springs' start
      %points. Springs' end points are the same with switched sign
      hookeanForces = -k.*springDisplacements;
      %spring unit vectors
      %  this may be needed in some cases
      %  springUnitVectors = zeros(size(springVectors));
      %  notzero = springLengths~=0;
      %  springUnitVectors(notzero,:) = bsxfun(@rdivide, springVectors(notzero,:), springLengths(notzero));
      springUnitVectors = bsxfun(@rdivide, springVectors, springLengths);
      %relative speeds of springs' ends
      springSpeeds = vel(springEnds(:,2),:)-vel(springEnds(:,1),:);
      %relative speeds of springs' ends along springVector direction:
      %projection of springSpeeds into springVectors: dot product with unit
      %vector
      springDirectSpeeds = sum(springSpeeds.*springUnitVectors,2);
      %damping forces for each spring. Note that they aren't negated
      %as it may seem sound, because of the chosen reference frame
      dampForces = c .* springDirectSpeeds;
      %vector force for each spring: hookean + damp forces
      springVectorForces = bsxfun(@times, springUnitVectors, hookeanForces+dampForces);
      %compute hookean forces for each dynamic point:
      %spring forces for first ends minus spring forces for second ends
      compressedSprings = springDisplacements>0;
      stretechedSprings = ~compressedSprings;
      compressiveForces = springsMatrix(:,compressedSprings)*springVectorForces(compressedSprings,:);
      stretechedForces  = springsMatrix(:,stretechedSprings)*springVectorForces(stretechedSprings,:);
    %   %profilaxis, it pays to locate errors inside this function, before
    %   %they spread over the code
    %   if any(isnan(forces(:)))
    %       error('Some point''s force is NaN');
    %   end;
    %   if ~isreal(forces)
    %       error('Some point''s force is complex');
    %   end
    end


