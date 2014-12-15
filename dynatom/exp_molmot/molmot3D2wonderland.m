function scene = molmot3D2wonderland(params, wonderparams, ss, rec, indexes)
% class Element:
% 	name = 0
% 	edgetype = 0
% 	initframe = 0
% 	graphs = 0
% 	entities = 0
% 	textures = 0
% 
% class Entidad:
% 	entname = 0
% 	entkind = 0
% 	entidxs = 0
%         entlayers = 0
% 
% class Graph:
% 	nodes = 0
% 	rads  = 0
% 	edges = 0
% 	numedges = 0
% 	edgesTextures = 0
% 	edgesSpringConst = 0

Ts = allTimeSteps(rec);
if isempty(indexes)
  indexes = 1:numel(Ts);
end
Ts = Ts(indexes);

empt     = cell(size(Ts));
index    = 1;
factor   = 1;10;
factorRow = 2;
factorMol = 2;
factorSpr = 2;
bc       = 0;

springs = reshape(find(ss.r>=(params.evparams.walker.pointAllr*2)), 1, []);

%rowBalls = -20:10;
rowBalls = -50:15;

rowpos        = [factor*rowBalls(:)*params.evparams.walker.row.ballSep, zeros(numel(rowBalls), 2)];

molGraphs     = struct('nodes', empt,   'rads', empt, 'edges', empt);
rowGraph      = struct('nodes', rowpos, 'rads', [],   'edges', []);

npoints       = numel(ss.m);% params.genome.fold.numPoints;

rowEntity     = struct('entid', 'row balls',    'entname', 'R1', 'entkind', 0, 'entidxs', 1:numel(rowBalls),              'entlayers', 1,                       'entrad', factorRow*factor*params.evparams.walker.row.ballAllr(1, ones(1,numel(rowBalls))) );
molEntity     = struct('entid', 'mol balls',    'entname', 'M1', 'entkind', 0, 'entidxs', 1:npoints,                      'entlayers', ones(1, numel(indexes)), 'entrad', factorMol*factor*params.evparams.walker.pointAllr(1,    ones(1,npoints))         );
sprEntity     = struct('entid', 'mol springs',  'entname', 'M1', 'entkind', 2, 'entidxs', springs,                        'entlayers', ones(1, numel(indexes)), 'entrad', factorSpr*factor*params.evparams.walker.pointAllr );
atp1Entity    = struct('entid', 'atp1 ball',    'entname', 'A1', 'entkind', 1, 'entidxs', zeros(1, numel(indexes)),       'entlayers', ones(1, numel(indexes)), 'entrad', factorMol*factor*params.evparams.walker.pointAllr );
atp2Entity    = struct('entid', 'atp2 ball',    'entname', 'A2', 'entkind', 1, 'entidxs', zeros(1, numel(indexes)),       'entlayers', ones(1, numel(indexes)), 'entrad', factorMol*factor*params.evparams.walker.pointAllr );
s11Entity     = struct('entid', 'atp11 spring', 'entname', 'A1', 'entkind', 3, 'entidxs', zeros(1, numel(indexes)),       'entlayers', ones(1, numel(indexes)), 'entrad', factorSpr*factor*params.evparams.walker.pointAllr );
s12Entity     = struct('entid', 'atp12 spring', 'entname', 'A1', 'entkind', 3, 'entidxs', zeros(1, numel(indexes)),       'entlayers', ones(1, numel(indexes)), 'entrad', factorSpr*factor*params.evparams.walker.pointAllr );
s21Entity     = struct('entid', 'atp21 spring', 'entname', 'A2', 'entkind', 3, 'entidxs', zeros(1, numel(indexes)),       'entlayers', ones(1, numel(indexes)), 'entrad', factorSpr*factor*params.evparams.walker.pointAllr );
s22Entity     = struct('entid', 'atp22 spring', 'entname', 'A2', 'entkind', 3, 'entidxs', zeros(1, numel(indexes)),       'entlayers', ones(1, numel(indexes)), 'entrad', factorSpr*factor*params.evparams.walker.pointAllr );

playSimulation(rec, Ts, @translate3D);
fprintf('\n');

molElementSpr = struct('name', {'MOL1'}, 'edgetype', {'TYPE1'}, 'initTime', {1}, 'graphs', molGraphs, 'entities', [molEntity, atp1Entity, atp2Entity, sprEntity, s11Entity, s12Entity, s21Entity, s22Entity]);
molElement    = struct('name', {'MOL1'}, 'edgetype', {'TYPE1'}, 'initTime', {1}, 'graphs', molGraphs, 'entities', [molEntity, atp1Entity, atp2Entity]);
rowElement    = struct('name', {'ROW1'}, 'edgetype', {'TYPE2'}, 'initTime', {1}, 'graphs', rowGraph,  'entities', rowEntity);

if ~isempty(wonderparams) && isfield(wonderparams, 'springs') && wonderparams.springs
  scene       = struct('elements', {[rowElement, molElementSpr]},    'framesByTimeUnit', {5});
else
  scene       = struct('elements', {[rowElement, molElement]},       'framesByTimeUnit', {5});
end


  function translate3D(t, ss) %#ok<INUSL>
    if bc>0; fprintf(repmat('\b', 1, bc)); end
    bc = fprintf('Doing t=%s (%d of %d)...', mat2str(t), index, numel(indexes));
    molGraphs(index).nodes = ss.pos*factor;
    molGraphs(index).rads  = [];%ss.stick.allr*factor;
    molGraphs(index).edges = ss.springEnds;
    atp1 = ss.legs(1).atp;
    if isempty(atp1)
      atp1Entity.entlayers(index) = 2;
      atp1Entity.entidxs(index)   = 1;
      s11Entity.entlayers(index)  = 2;
      s11Entity.entidxs(index)    = 1;
      s12Entity.entlayers(index)  = 2;
      s12Entity.entidxs(index)    = 1;
    else
      atp1Entity.entlayers(index) = 1;
      atp1Entity.entidxs(index)   = atp1;
      s11Entity.entlayers(index)  = 1;
      s11Entity.entidxs(index)    = ss.legs(1).atpb(1);
      s12Entity.entlayers(index)  = 1;
      s12Entity.entidxs(index)    = ss.legs(1).atpb(2);
    end
    atp2 = ss.legs(2).atp;
    if isempty(atp2)
      atp2Entity.entlayers(index) = 2;
      atp2Entity.entidxs(index)   = 1;
      s21Entity.entlayers(index)  = 2;
      s21Entity.entidxs(index)    = 1;
      s22Entity.entlayers(index)  = 2;
      s22Entity.entidxs(index)    = 1;
    else
      atp2Entity.entlayers(index) = 1;
      atp2Entity.entidxs(index)   = atp2;
      s21Entity.entlayers(index)  = 1;
      s21Entity.entidxs(index)    = ss.legs(2).atpb(1);
      s22Entity.entlayers(index)  = 1;
      s22Entity.entidxs(index)    = ss.legs(2).atpb(2);
    end
    index = index+1;
  end

end