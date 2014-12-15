classdef ElegantMolPlotter < BasicPlotter
  properties
    sphereDetail;
    sphereRad;
    sphereColor;
    cylinderDetail;
    cylinderRad;
    cylinderColor;
    atpColor;
    atpLinkColor;
    toeColor;
    jointIndex;
    jointColor;
    rowRad;
    rowSegment;
    rowColor;
    rowSpec;
    moveCM;
    lastCM;
    shiftCamera;
    light;
    toVaryPoint;
    showBalls; %to have compatbility with BallPlotter
  end

  properties (GetAccess = protected, SetAccess = protected)
  end
  
  methods
    %constructor
    function p = ElegantMolPlotter(varargin)
      p = p@BasicPlotter;
      p.sphereDetail        = 20;
      p.sphereRad           = 3.8/2;
      p.sphereColor         = 0.3;
      p.cylinderDetail      = 10;
      p.cylinderRad         = 0.3;
      p.cylinderColor       = 0.4;
      p.atpColor            = 0;
      p.atpLinkColor        = 0;
      p.jointIndex          = [];
      p.jointColor          = 0.4;
      p.toeColor            = 0.6;
      p.rowRad              = 3.8/2;
      p.rowSegment          = -30:30;
      p.rowColor            = 0;
      p.rowSpec             = 1;
      p.moveCM              = true;
      p.lastCM              = [];
      p.toVaryPoint         = false(1, 3);
      p.shiftCamera         = [0 0 0];
      if numel(varargin)>0; set(p, varargin{:}); end
    end
    
    function drawFrame(p, t, ss)
      drawFrame@BasicPlotter(p, t, ss);
      

          ns  = p.sphereDetail;
          rs  = p.sphereRad;
          cs  = p.sphereColor;

          nc  = p.cylinderDetail;
          rc  = p.cylinderRad;
          cc  = p.cylinderColor;
          
          ca  = p.atpColor;
          cl  = p.atpLinkColor;
          
          rr  = p.rowRad;
          sr  = p.rowSegment;
          if numel(sr)==2
            sr=sr(1):sr(2);
          end
          cr  = p.rowColor+0.2; 
          specr = p.rowSpec;


          toVaryPoint = p.toVaryPoint;
          %[70.43 70.44 120.45 122 122.87 142.87
      [X,Y,Z] = sphere(ns);
      if toVaryPoint(1) %#ok<PROP>
        X=-X;
      end
      if toVaryPoint(2) %#ok<PROP>
        Y=-Y;
      end
      if toVaryPoint(3) %#ok<PROP>
        Z=-Z;
      end
      [F,V]   = surf2patch(X,Y,Z);

      st = struct('vertices', [], 'faces', F);

      args = {...
        'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
        'SpecularStrength', 0, 'DiffuseStrength', 1, ...
        'AmbientStrength', 1, ...
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

      if isfield(ss, 'legs') 
        toe  = p.toeColor;
        for k=1:numel(ss.legs)
          for z=1:numel(ss.legs(k).toes)
            t = ss.legs(k).toes(z);
            set(vertices{t}, 'FaceColor', toe*ones(3,1));
          end
          for z=1:numel(ss.legs(k).atp)
            set(vertices{ss.legs(k).atp(z)}, 'FaceColor', ca*ones(3,1));
          end
          for z=1:numel(ss.legs(k).atpb)
            set(links{ss.legs(k).atpb(z)}, 'FaceColor', cl*ones(3,1));
          end
        end
      end
      
      ji = p.jointIndex;
      jc = p.jointColor;
      for k=1:numel(ji)
        set(vertices{ji(k)}, 'FaceColor', jc*ones(3,1));
      end
      
      [X,Y,Z] = sphere(ns);
      if toVaryPoint(1) %#ok<PROP>
        X=-X;
      end
      if toVaryPoint(2) %#ok<PROP>
        Y=-Y;
      end
      if toVaryPoint(3) %#ok<PROP>
        Z=-Z;
      end
      [F,V]   = surf2patch(X,Y,Z);
      rowc = [sr(:)*rr*2, zeros(numel(sr), 2)];
      row = cell(size(sr(:)));
      for k=1:numel(sr)
        st.vertices = bsxfun(@plus, V*rr, rowc(k,:));
        r = patch(st);
        set(r, 'FaceColor', cr*ones(3,1), args{:}, 'SpecularStrength', specr, 'DiffuseStrength', 1, 'AmbientStrength', 1);
        row{k} = r;
      end
      
      if p.moveCM
         w = ss.m/sum(ss.m);
         cm = sum(bsxfun(@times, ss.pos, w));
         camBefore = (~isempty(p.camProps)) && strcmp(p.camMode, 'before');
         ax      = get(p.figura, 'CurrentAxes');
         if isempty(p.lastCM)
           camtarget(get(p.figura, 'CurrentAxes'), cm);
           if camBefore
             p.lastCM=cm;
           end
         else
           veccmcm = cm-p.lastCM;
           camtarg = get(ax, 'CameraTarget');
           camp    = get(ax, 'CameraPosition');
           set(ax, 'CameraTarget',   camtarg+veccmcm);
           set(ax, 'CameraPosition',    camp+veccmcm);
         end
         camtarg = get(ax, 'CameraTarget');
         camp    = get(ax, 'CameraPosition');
         set(ax, 'CameraTarget',   camtarg+p.shiftCamera);
         set(ax, 'CameraPosition',    camp+p.shiftCamera);
         if ~camBefore
           p.lastCM=cm;
         end
      end
         
         
         p.light = camlight('headlight', 'infinite');
         axis off;
         
    end
    
    function initPlotting(p)
      initPlotting@BasicPlotter(p);
      p.lastCM = [];
      cameratoolbar(p.figura, 'Show');
    end
    
    function txt = makeSpringText(p, idx, st, p1, p2)
      txt = makeSpringText@SpringPlotter(p, idx, st, p1, p2);
%       if isfield(st, 'rdyn')
%         txt = [txt sprintf('\nr0: %-1.3f rMOD: %-1.3f tbaseMOD: %-1.3f decayMOD: %-1.3f realMOD: %-1.3f', st.rdyn.r0(idx), st.rdyn.rmod(idx), st.rdyn.rtb(idx), st.rdyn.rtm(idx), updateMod(st.t, st.rdyn.rmod(idx), st.rdyn.rtb(idx), st.rdyn.rtm(idx)))];
%       end
    end
    
    function fun = setTargetFun(p, cm)
      fun = @(varargin)set(get(p.figura, 'CurrentAxes'), 'CameraTarget', cm);
    end
    
  end
  
end
