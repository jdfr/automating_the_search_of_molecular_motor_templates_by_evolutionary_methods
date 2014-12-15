classdef EigenPlotter < SpringPlotter
  
  properties
    change;
  end
  
  methods
    function drawFrame(p, t, st)
      drawFrame@SpringPlotter(p, t, st);
      change = p.change;
      eivt = '';
      if isfield(st.stick, 'eigVector') && (~isempty(st.stick.eigVector))
        eiv = st.stick.eigVector;
        vec2    = realpow(eiv, 2);
        vec2s   = vec2/sum(vec2);
        ka      = 1/size(eiv, 1)*exp(-sum(vec2s.*log(vec2s)));
        eivt    = sprintf('\ndegree of collectivity: %s', mat2str(ka));
      else
        if isfield(st.stick, 'eigVectorss') && (~isempty(st.stick.eigVectorss))
          eiv     = st.stick.eigVectorss;
          vec2    = sum(realpow(eiv, 2), 2);
          vec2s   = vec2/sum(vec2);
          ka      = 1/size(eiv, 1)*exp(-sum(vec2s.*log(vec2s)));
          eivsp   = eiv(st.springEnds(:,1),:)-eiv(st.springEnds(:,2),:);
          vecsp2  = sum(realpow(eivsp, 2), 2);
          vecsp2s = vecsp2/sum(vecsp2);
          kb      = 1/size(eivsp, 1)*exp(-sum(vecsp2s.*log(vecsp2s)));
          eivt    = sprintf('\ndegree of collectivity A: %s\ndegree of collectivity B: %s', mat2str(ka), mat2str(kb));
        end
      end
      if (numel(change)==1) && islogical(change) && change
        if (size(st.pos,2)>2) && isfield(st.stick, 'eigVector') && (~isempty(st.stick.eigVector))
          [zx zy zz] = sphere(3);
          fvc=surf2patch(zx, zy, zz);
          for k=1:size(st.pos,1)
            fvc2 = fvc;
            mag = st.stick.eigVector(k);
            fvc2.vertices = bsxfun(@plus, fvc2.vertices*abs(mag)*10, st.pos(k,:));
            if mag>0
              col = 'm';
            else
              col = 'k';
            end 
            patch(fvc2, 'FaceColor', col);
          end
        else
          if ~isempty(st.stick.eigVec1)
            line(st.stick.eigVec1(:,1), st.stick.eigVec1(:,2), 'Color', 'm');
          end
          if ~isempty(st.stick.eigVec2)
            line(st.stick.eigVec2(:,1), st.stick.eigVec2(:,2), 'Color', 'k');
          end
        end
      else
        if ~isempty(st.stick.eigVec1)
          ps = zeros(3*size(st.pos,1), size(st.pos,2));
          ps(1:3:end,:) = st.pos;
          ps(2:3:end,:) = st.pos-st.stick.eigVec1;
          ps(3:3:end,:) = nan;
          args = {'Color', 'm', 'lineWidth', 1.5};
          if size(ps,2)<3
            line(ps(:,1), ps(:,2),          args{:});
          else
            line(ps(:,1), ps(:,2), ps(:,3), args{:});
          end
        end
      end
      drawTextInFigure(p, [mat2str(st.stick.eigTxt) eivt]);
    end
  end
end

