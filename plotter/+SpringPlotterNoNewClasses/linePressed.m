%retrieve info for clicked spring
function linePressed(p, obj)
  %delete previous point markers
  if ~isempty(p.end1); delete(p.end1); end;
  if ~isempty(p.end2); delete(p.end2); end;
  %retrive index and structure
  data = get(obj, 'UserData');
  idx = data{1};
  st  = data{2};

  %change selected spring by default
  p.selectedSpringByDefault = idx;

  %plot marker for point 1
  p1 = st.springEnds(idx,1);
  pos = mat2cell(st.pos(p1,:), 1, ones(1,size(st.pos,2)));
  p.end1 = line(pos{:},'MarkerFaceColor','k','MarkerEdgeColor', 'k', 'MarkerSize', 8, 'Marker', 'o');

  %plot marker for point 2
  p2 = st.springEnds(idx,2);
  pos = mat2cell(st.pos(p2,:), 1, ones(1,size(st.pos,2)));
  p.end2 = line(pos{:},'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerSize', 8, 'Marker', 'o');

  %write info
  set(p.textLineInfo, 'String', makeSpringText(p, idx, st, p1, p2));%, 'FontName', 'FixedWidth', 'FontSize', 8); %, 'HorizontalAlignment', 'left');
end


    %format information for a spring
    function txt = makeSpringText(p, idx, st, p1, p2)
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
