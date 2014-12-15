%set up a frame for drawing it
function p = initFrame(p, st)
  %erase axes' contents
  cla(get(p.figura, 'CurrentAxes'));
  hold on;                               
  %set up axes
  if p.fixedAxisFrame
    if numel(p.axisWindow)>1
      axis(p.axisWindow);
    else
      switch size(st.pos,2)
        case 2; axis([-p.axisWindow p.axisWindow -p.axisWindow p.axisWindow]);
        case 3; axis([-p.axisWindow p.axisWindow -p.axisWindow p.axisWindow -p.axisWindow p.axisWindow]);
      end
    end
  end;
  %the first time we set the viewpoint
  if p.justStarted 
    if ~isempty(p.viewPoint)
      view(p.viewPoint);
    else
      switch size(st.pos,2)
        case 2; view(2);
        case 3; view(3);
      end
    end
    p.justStarted = false;
  end
  if p.axisEqual;  axis equal;  end
  if p.axisSquare; axis square; end
  if p.hasAxis; axis on; else axis off; end;
  if p.hasGrid; grid on; else grid off; end;
end
    
