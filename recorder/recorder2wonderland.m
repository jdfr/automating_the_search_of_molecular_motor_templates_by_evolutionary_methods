function scene = recorder2wonderland(rec, indexes)

Ts = allTimeSteps(rec);
Ts = Ts(indexes);

empt     = cell(size(Ts));
index    = 1;

graphs   = struct('nodes', empt, 'edges', empt, 'rads', empt);%, 'edgesSpringConst', {empt});

switch size(rec.notdyn{1}.pos, 1)
  case 2
    fun = @translate2D;
  case 3
    fun = @translate3D;
end

playSimulation(rec, Ts, fun);

element  = struct('name', {'ELE1'}, 'edgetype', {'TYPE1'}, 'initTime', {2}, 'graphs', graphs);

scene    = struct('elements', {{element}}, 'framesByTimeUnit', {5});

  function translate2D(t, ss) %#ok<INUSL>
    graphs(index).nodes = [ss.pos*10, zeros(size(ss.pos,1),1)];
    graphs(index).rads  = (ss.rad+ss.stick.dist)*10;
    graphs(index).edges = ss.springEnds;
    index = index+1;
  end
  function translate3D(t, ss) %#ok<INUSL>
    graphs(index).nodes = ss.pos*10;
    graphs(index).rads  = ss.stick.allr*10;
    graphs(index).edges = ss.springEnds;
    index = index+1;
  end

end