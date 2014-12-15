classdef SpringPlotter < BasicPlotter
  properties (GetAccess = public, SetAccess = public)
    showLineInfo;            %whether to show springs and points info
    selectedSpringByDefault; %if showLineInfo==true, this var tells what spring is selected by default
    acelFactor;              %if 0, does nothing. Otherwise, shows acelerations for points, scaled by a factor 'acelFactor'
    showSpringAcels;         %if acelFactor>0, shows spring acels (compressed and stretched springs) 
    showOtherAcels;          %if acelFactor>0, shows other acels (drag force, perturbations, enforcements and customized ones) 
    showTotalAcels;          %if acelFactor>0, shows total acels (sum of the previous ones) 
    velFactor;               %if 0, does nothing, Otherwise, shows scaled vels for points
    shownSprings;            %indexes of drawn springs. If empty, show all of them
    shownPoints;             %indexes of drawn points.  If empty, show all of them
    highlightChangedSprings; %highlight springs that are to change
    colorsCodeForce;         %if true, colors code force. If false, colors code relative length and spring stiffness.
    maxForce;                %when colorsCodeForce=true: the range of spring forces to be mapped in [blue...green...red] will be [maxForce, -maxForce]
    springMaxRel;            %when colorsCodeForce=false: the range of spring stresses to be mapped in a map [blue...green...red] will be [-springMaxRel, springMaxRel]
    rangeK;                  %when colorsCodeForce=false: rangeK=[min max] is a range for springs' K values: springs with K=max or higher will be displayed fully saturated; springs with K=min or lower will be displayed minimally saturated
    minSat;                  %when colorsCodeForce=false: between 0 and 1: minimum saturation, for springs with K=rangeK(1) (the maximum saturation is 1)
    springMaxForce;
    showMode;        %'analogic' 'digital'
    digitalDimensions;
  end

  properties (GetAccess = protected, SetAccess = protected)
    textLineInfo; %text uicontrol to show info
    end1; %selected spring's point 1
    end2; %selected spring's point 2
    buttonDownFcn; %callback for selected springs
  end
  
  methods
    %constructor
    function p = SpringPlotter(varargin)
      p = p@BasicPlotter;
      p.showLineInfo = true;
      p.acelFactor = 0;
      p.showSpringAcels = true;
      p.showOtherAcels  = false;
      p.showTotalAcels  = false;
      p.velFactor       = 0;
      p.selectedSpringByDefault = [];
      p.shownSprings = [];
      p.shownPoints  = [];
      p.colorsCodeForce = false;%true;
      p.maxForce     = 0.5;
      p.springMaxRel = 0.03;
      p.rangeK       = [0.5 1];
      p.minSat       = 1/3;
      p.highlightChangedSprings = false;
      p.showMode     = 'analogic';
      p.digitalDimensions = [500 500];
      if numel(varargin)>0; set(p, varargin{:}); end
      p.buttonDownFcn = @(obj, event) linePressed(p, obj);
    end
    
    %selectedSpringByDefault may imply showLineInfo
    function p = set.selectedSpringByDefault(p,v)
      p.selectedSpringByDefault = v;
      if ~isempty(p.selectedSpringByDefault)
        p.showLineInfo = true;
      end
    end
    
    function drawSSAndText(p, st, txt)
      drawSS(p, st);
      if p.showLineInfo
        set(p.textLineInfo, 'String', txt);
      end
    end
    
    function drawTextInFigure(p, txt)
      if p.showLineInfo && (~isempty(p.textLineInfo))
        set(p.textLineInfo, 'String', txt);
      end
    end
    
    %redefine parent's drawFrame to plot springs
    function drawFrame(p, t, st)
      drawFrame@BasicPlotter(p, t, st);
      if (p.acelFactor~=0)
        [springDisplacements, compressedSpringAcels, stretechedSpringAcels, otherAcels, dynamicVel, changedDynamicVel] = SpringPlotter.getKinematicVectors(t, st);
        dynamicVel = changedDynamicVel;
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
        dynamicVel = dynamicVel*p.velFactor;
        plotPointKinematics(p, t, st, [], [], [], dynamicVel);
      else
        plotSprings(p, t, st);
      end
    end
    
    %redefine parent's initPlotting to create, if necessary uicontrols for
    %showing info
    function initPlotting(p)
      initPlotting@BasicPlotter(p);
      switch p.showMode
        case {'analog', 'analogic'}
          if p.showLineInfo
            set(p.figura,'Toolbar','figure');
            p.textLineInfo = uicontrol('Style', 'text', 'Units', 'normalized', ...
            'FontName', 'FixedWidth', 'Position', [0.1 0 0.9 0.07]);
          end
        case 'digital'
        otherwise
          error('showMode <<%s>> not understood!!!!', mat2str(p.showMode));
      end
    end
    
    function plotPointKinematics(p, t, st, compressedSpringAcels, stretechedSpringAcels, otherAcels, dynamicVel) %#ok<INUSL>
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
        colors = squeeze(hsv2rgb(mappingHue*(2/3), ones(size(pseudoHookeanForces)), ones(size(pseudoHookeanForces))));
      else
        springRelDisplacements = springDisplacements./st.r(showSprings);
        springRelDisplacements(springRelDisplacements<-mxrel)=-mxrel;
        springRelDisplacements(springRelDisplacements>mxrel)=mxrel;
        mappingHue = (1-((springRelDisplacements-(-mxrel))/(mxrel-(-mxrel))));
        mappingSat = (mnSat-1)*(st.k(showSprings)-ranK(2))/(ranK(1)-ranK(2))+1;
        mappingSat(mappingSat<mnSat) = mnSat;
        mappingSat(mappingSat>1) = 1;
        colors = squeeze(hsv2rgb(mappingHue*(2/3), mappingSat, ones(size(springRelDisplacements))));
      end
%       mappingHue = (updateMod(t, st.rmod(showSprings), st.rtb(showSprings), st.rtm(showSprings))-1)*2;
%       colors = squeeze(hsv2rgb(mappingHue*(2/3), ones(size(springDisplacements)), ones(size(springDisplacements))));
%       colors(colors<0) = 0;
%       colors(colors>1) = 1;
      if size(colors,2)==1
        colors = colors';
      end
      
      switch p.showMode(1)
        case 'd'
          aw = p.axisWindow;
          if numel(aw)==1
            aw = [-aw aw -aw aw];
          end
          imgsiz = p.digitalDimensions;
          if numel(imgsiz)==1
            imgsiz = [imgsiz imgsiz];
          end
          img = drawAsImage(st.pos, st.springEnds, aw, imgsiz, colors);
          axis([0 imgsiz(1)-1 0 imgsiz(2)-1]); axis ij; imshow(img);
        case 'a'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %if we must show info, reset points
          if p.showLineInfo
            p.end1 = [];
            p.end2 = [];
          end
          %draw springs
          st2 = st; st2.t = t;
          if islogical(showSprings)
            showSprings = find(showSprings);
          end
          for i=1:numel(showSprings)
              indexes = st.springEnds(showSprings(i),:);
              if p.showLineInfo
                infoData = {'UserData', {showSprings(i), st2}, 'ButtonDownFcn', p.buttonDownFcn};
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
          %black & striped for redesignFixK
          %gray  & striped for redesignFixR
          %black & solid   for redesignFactorK
          %gray  & slid    for redesignFactorR
          if p.highlightChangedSprings
            colorA = 'k'; colorB = [0.8 0.8 0.8]; linestyleA = '-'; linestyleB = '-.';
            for k=1:numel(st.atomState)
              stateAS                        = st.atomState{k}.afterStabilized;
              %TODO: hack to avoid displaying the whole-organism modification
              if (numel(stateAS.redesignListTypes)>0) && (stateAS.redesignListTypes(1)==st.geneTypes.redesignFactorR) && (abs(stateAS.redesignListArgs{1}(1)-2.1449451076582799)<1e-6) && (stateAS.redesignListArgs{1}(2)==2)
                start = 2;
              else
                start = 1;
              end
              for m=start:numel(stateAS.redesignListTypes)
                %fprintf('cosa generations %d\n', stateAS.generations(m));
    %             if stateAS.generations(m)==0
    %               continue;
    %             end
                args = stateAS.redesignListArgs{m};
                switch stateAS.redesignListTypes(m)
                  case st.geneTypes.redesignFixK
                    springids = args(3:end);
                    %fprintf('redesignFixK atom %03d gene %03d spring %s\n', k, m, mat2str(springids));
                    for s=1:numel(springids);
                      spring = st.atomSprings(k,springids(s));
                      if spring ==0; continue; end
                      indexes = st.springEnds(spring,:);
                      switch st.ndims
                        case 2
                          h = line(st.pos(indexes,1),st.pos(indexes,2), 'Color',colorA, 'LineWidth', 2, 'LineStyle', linestyleB); %#ok<NASGU>
                      end
                    end
                  case st.geneTypes.redesignFixR
                    springids = args(2:end);
                    %fprintf('redesignFixR atom %03d gene %03d spring %s\n', k, m, mat2str(springids));
                    for s=1:numel(springids);
                      spring = st.atomSprings(k,springids(s));
                      if spring ==0; continue; end
                      indexes = st.springEnds(spring,:);
                      switch st.ndims
                        case 2
                          h = line(st.pos(indexes,1),st.pos(indexes,2), 'Color',colorB, 'LineWidth', 2, 'LineStyle', linestyleB); %#ok<NASGU>
                      end
                    end
                  case st.geneTypes.redesignFactorK
                    springids = args(2:end);
                    %fprintf('redesignFactorK atom %03d gene %03d spring %s\n', k, m, mat2str(springids));
                    for s=1:numel(springids);
                      spring = st.atomSprings(k,springids(s));
                      if spring ==0; continue; end
                      indexes = st.springEnds(spring,:);
                      switch st.ndims
                        case 2
                          h = line(st.pos(indexes,1),st.pos(indexes,2), 'Color',colorA, 'LineWidth', 2, 'LineStyle', linestyleA); %#ok<NASGU>
                      end
                    end
                  case st.geneTypes.redesignFactorR
                    springids = args(2:end);
                    %fprintf('redesignFactorR atom %03d gene %03d spring %s\n', k, m, mat2str(springids));
                    for s=1:numel(springids);
                      spring = st.atomSprings(k,springids(s));
                      if spring ==0; continue; end
                      indexes = st.springEnds(spring,:);
                      switch st.ndims
                        case 2
                          h = line(st.pos(indexes,1),st.pos(indexes,2), 'Color',colorB, 'LineWidth', 2, 'LineStyle', linestyleA); %#ok<NASGU>
                      end
                    end
    %               case st.geneTypes.redesignEraseSpring
    %                 springids = args(1:end);
    %                 %fprintf('redesignEraseSpring atom %03d gene %03d spring %s\n', k, m, mat2str(springids));
    %                 for s=1:numel(springids);
    %                   spring = st.atomSprings(k,springids(s));
    %                   if spring ==0; continue; end
    %                   indexes = st.springEnds(spring,:);
    %                   switch st.ndims
    %                     case 2
    %                       h = line(st.pos(indexes,1),st.pos(indexes,2), 'Color','m', 'LineWidth', 2, 'LineStyle', ':');
    %                   end
    %                 end
                end
              end
            end
          end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        otherwise
          error('showMode <<%s>> not understood!!!!', mat2str(p.showMode));
      end
    end
    
    %format information for a spring
    function txt = makeSpringText(p, idx, st, p1, p2) %#ok<INUSL>
      springVector = st.pos(p2,:)-st.pos(p1,:);
      springLength = realsqrt(sum(realpow(springVector,2), 2));
      %springDisplacement = ss.r(idx)-springLength;
      %spring length relative to its rest length
      relative = st.r(idx)./springLength;
      hookeanForce = -st.k(idx).*(st.r(idx)-springLength);
      springUnitVector = springVector./springLength;
      springSpeed = st.vel(p2,:)-st.vel(p1,:);
      springDirectSpeed = sum(springSpeed.*springUnitVector,2);
      dampForce = st.c(idx) .* springDirectSpeed;
      force = hookeanForce+dampForce;
      pr = @(x) sprintf('%+.3e ', x);
      txt = sprintf('Line %-4g len=%.3e r=%+.3e (rel=%-1.3f) k=%+.3e c=%+.3e F=%s\nblack:p1 %-4g pos=%svel=%sm=%+.3e\n gray:p2 %-4g pos=%svel=%sm=%+.3e', ...
        idx, springLength, st.r(idx), relative, st.k(idx), st.c(idx), pr(force), ...
        p1, pr(st.pos(p1,:)), pr(st.vel(p1,:)), st.m(p1), ...
        p2, pr(st.pos(p2,:)), pr(st.vel(p2,:)), st.m(p2));
    end

    %retrieve info for clicked spring
    function linePressed(p, obj)
      %delete previous point markers
      if ~isempty(p.end1); delete(p.end1); end;
      if ~isempty(p.end2); delete(p.end2); end;
      %retrive index and structure
      if iscell(obj)
        data = obj;
      else
        data = get(obj, 'UserData');
      end
      idx = data{1};
      st  = data{2};
      writeText = (numel(data)<3) || (islogical(data{3}) && data{3});
      if (~writeText) && (numel(idx)==2)
        p1 = idx(1);
        p2 = idx(2);
      else
        p1 = st.springEnds(idx,1);
        p2 = st.springEnds(idx,2);
      end
      
      %change selected spring by default
      p.selectedSpringByDefault = idx;

      %plot marker for point 1
      pos = mat2cell(st.pos(p1,:), 1, ones(1,size(st.pos,2)));
      p.end1 = line(pos{:},'MarkerFaceColor','k','MarkerEdgeColor', 'k', 'MarkerSize', 8, 'Marker', 'o');
      
      %plot marker for point 2
      pos = mat2cell(st.pos(p2,:), 1, ones(1,size(st.pos,2)));
      p.end2 = line(pos{:},'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerSize', 8, 'Marker', 'o');
      
      %write info
      if writeText
        set(p.textLineInfo, 'String', makeSpringText(p, idx, st, p1, p2));%, 'FontName', 'FixedWidth', 'FontSize', 8); %, 'HorizontalAlignment', 'left');
      end
    end
  end
  methods (Static)
    %decompose system's dynamics
    function [springDisplacements, compressedSpringAcels, stretechedSpringAcels, otherAcels, dynamicVel, changedDynamicVel] = getKinematicVectors(t, st)
      st = prepareSimulationStructure(st);
      dynamicVel = st.vel(st.dynamical_p,:);
      masses = st.m(st.dynamical_p);
      [compressiveForces, stretechedForces, springDisplacements] = SpringPlotter.calculateDetailedSpringForces(st);
      compressedSpringAcels = bsxfun(@rdivide, compressiveForces, masses);
      stretechedSpringAcels = bsxfun(@rdivide, stretechedForces, masses);
      dragAcels = bsxfun(@rdivide, -bsxfun(@times, dynamicVel, st.u), masses);
%       [changedDynamicVel AfterAcel] = nonNaturalDynamics(st, t, dynamicVel, compressedSpringAcels+stretechedSpringAcels+dragAcels, [], [], [], [], []);
%       otherAcels = AfterAcel-compressedSpringAcels-stretechedSpringAcels;
      otherAcels = zeros(size(dragAcels));
      changedDynamicVel = dynamicVel;
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
    
  end
  
end
