function ss = structure2walker(ss, evparams)

if ~isfield(evparams.walker, 'fixUnixWindows')
  evparams.walker.fixUnixWindows.doFix = false;
end

%numver of vertices
N = size(ss.pos, 1);

%connectivity matrix
conns = sparse(ss.springEnds(:), reshape(ss.springEnds(:,[2,1]), [], 1), true, N, N);

%kirchoff/laplacian matrix
kirchoff = sparse(1:N, 1:N, sum(conns), N,N) - conns;

if isfield(evparams.walker, 'eigCalc')
  eigCalc = evparams.walker.eigCalc;
else
  eigCalc = 'sparse';
end
switch eigCalc
  case 'sparse'
    %eigendecomposition
    [vecs vals flag] = eigs(kirchoff, 3, 'sa', struct('disp', 0));
    if flag~=0
      ss.legs = sprintf('Error while calculating legs: The eigendecomposition did not converge!!! (flag==%s)', mat2str(flag));
      return
    end
    %third mode; it should decompose the structure in 3 parts
    [idx3 idx3] = sort(diag(vals));
    vec3 = vecs(:,idx3(end));
  case 'full'
    [vecs vals] = eig(full(kirchoff));
    [idx3 idx3] = sort(diag(vals));
    vec3 = vecs(:,idx3(3));
  otherwise
    error('Parameter eigCalc=%s not understood!!!\n', any2str(eigCalc));
end

% %eigendecomposition
% [eivecs eivals] = eig(kirchoff);
% %third mode; it should decompose the structure in 3 parts
% vec3 = eivecs(:,3);

%find the 3 parts
[regions numregions] = findConnectedRegionsBySigns(vec3, conns, evparams);

if numregions ~= 3
  ss.legs = sprintf('Error while calculating legs: The shape of the graph is not lineal enough! There are not 3 regions in the 3rd mode, but %d!!!', numregions);
  return
end

%find the 2 interfaces between the 3 parts
interfaces = findInterfaces(regions, numregions, conns);

if numel(interfaces) ~= 2
  ss.legs = sprintf('Error while calculating legs: The shape of the graph is bizarre! In spite of having 3 regions, it does not have 2 interfaces but %d!!!!!', numel(interfaces));
  return
end

toes = calculateToes(vec3, conns, interfaces, regions, evparams);

%we rotate early, so the orientation of the structure can be used to decide
%some leg parameters in the 3d case.
[ss toes contactToes] = rotateLegs(ss, evparams, toes);

if ischar(contactToes)
  ss.legs = contactToes;
  return
end

%from the decomposition, calculate the active sites of the structure (ATP
%binding cassete, toes)
ss.legs = makeLegs(interfaces, regions, vec3, conns, toes, contactToes, ss, evparams);


function [labels label] = findConnectedRegionsBySigns(autovec, conns, evparams)
%give an autovector of the laplacian matrix, and the connectivity matrix,
%segment the graph in several connected subgraphs whose vertices share the
%same sign by the autovector

signos              = sign(autovec);
if evparams.walker.fixUnixWindows.doFix
  %new treatment. Tries to minimize sensibility to changes in signs and
  %epsilon values
  epsilon           = evparams.walker.fixUnixWindows.epsilon;
  maxilon           = evparams.walker.fixUnixWindows.maxilon;
  autoabs           = abs(autovec);
  epsilons          = autoabs<epsilon;
  maxilon           = find(autoabs>maxilon, 1, 'first');
  if isempty(maxilon)
    maxilon         = find(autoabs>epsilon, 1, 'first');
  end
  signos(epsilons)  = signos(maxilon);
else
  %old treatment. Though arbitrary, this would be OK if not were for
  %eig(A) in unix might yield different results in windows and in Unix.
  signos(signos==0) = 1;
end

labels          = nan(size(signos));

nolabeled       = true(size(signos));

label            = 0;
while any(nolabeled)
  v              = find(nolabeled, 1);
  signo          = signos(v);
  
  while ~isempty(v)
    labels(v)    = label;
    nolabeled(v) = false;
    v            = find((any(conns(v,:), 1)') & nolabeled & (signos==signo));
  end
  
  label          = label+1;
end

function interfaces = findInterfaces(regions, numregions, conns)
%find the interfaces between the segmented regions

c0 = cell(0,1);
interfaces = struct('region1', c0, 'region2', c0, 'points1', c0, 'points2', c0);

for a=1:numregions-1
  ra = a-1;
  for b=a+1:numregions
    rb = b-1;
    pa = regions==ra;
    pb = regions==rb;
    subconns = conns(pa, pb);
    if any(subconns(:))
      pa = find(pa);
      pb = find(pb);
      [r c] = find(subconns);
      interfaces(end+1,1) = struct('region1', ra, 'region2', rb, 'points1', pa(unique(r)), 'points2', pb(unique(c))); %#ok<AGROW>
    end
  end
end

function legs = makeLegs(interfaces, regions, autovec, conns, toes, contactToes, ss, evparams)
%this expects 3 regions with 2 interfaces. The distal regions will be the
%"legs" of the motor machine, each having an ABC (ATP binding cassete) and
%a docking region (capable of contacting the microtubule). The ABC will be
%defined in the interface between regions

autovecabs = abs(autovec);

c2 = cell(2,1);
legs = struct('ABCweight', c2, 'ABCbinding', c2, 'toes', toes);

pos = ss.pos;
spEnds = ss.springEnds;

d3 = size(pos,2)>2;

correct                   = evparams.walker.correct;
if isfield(correct, 'sign')
  correctSign      = correct.sign;
else
  correctSign      = false;
end
if isfield(correct, 'singleRow')
  correctSingleRow = correct.singleRow;
else
  correctSingleRow = false;
end
if isfield(correct, 'pairIndex')
  correctPairIndex = correct.pairIndex;
else
  correctPairIndex = false;
end
if isfield(evparams.walker, 'useHull')
  useHull                 = evparams.walker.useHull;
else
  useHull                 = true;
end
if isfield(evparams.walker, 'extendInterface')
  extendInterface         = evparams.walker.extendInterface;
else
  extendInterface         = true;
end

usehoh = isfield(evparams.walker, 'd3') && strcmp(evparams.walker.d3.rotationMode, 'handoverhand');

if usehoh
  hoh       = ss.hoh;
  if any(cellfun('prodofsize', toes)==0)
    error('some toe set is empty!!! toes=%s', any2str(toes));
  end
  %center of mass of toes of leg 1
  t1 = toes{1};
  w1 = ss.m(t1)/sum(ss.m(t1));
  cm1 = sum(bsxfun(@times, ss.pos(t1,:), w1), 1);
  %center of mass of toes of leg 2
  t2 = toes{2};
  w2 = ss.m(t2)/sum(ss.m(t2));
  cm2 = sum(bsxfun(@times, ss.pos(t2,:), w2), 1);
  %select the lowest toes as the leg of choice
  if abs(cm1(2))<abs(cm2(2))
    idxlegs = 1;
    cm      = cm1;
  else
    idxlegs = 2;
    cm      = cm2;
  end
else
  idxlegs = [1 2];
end

switch evparams.walker.ABCMode
  case 'heurconcave'
    %calculate ABC parameters
    for z=1:numel(idxlegs)
      k = idxlegs(z);
      %points in the interface
      p1           = interfaces(k).points1(:);
      p2           = interfaces(k).points2(:);
      %connections in the interface
      connsi       = conns(p1, p2);
      %extended interface
      extended     = false(size(pos, 1), 1);
      extended(p1) = true;
      extended(p2) = true;
      if extendInterface
        %loop over connections
        for q=1:size(connsi,1)
          p1i = p1(q);
          neighs = find(connsi(q,:));
          if ~isempty(neighs)
            %for each connection from p1i to its neighbours, find its absolute
            %autovec values
            p2n        = p2(neighs);
            valuesAbs  = [repmat(autovecabs(p1i), 1, numel(neighs)); reshape(autovecabs(p2n), 1, [])];
            %for each connection, normalize the autovec values relative to the
            %biggest of the two
            values     = bsxfun(@rdivide, valuesAbs, max(valuesAbs));
            %collect min values (readable as ratios to the biggest) and indexes
            [minv, mn] = min(values);
            %determine the connections where the imbalance is too great
            toosmall   = minv<evparams.walker.tooSmallNeighbour;
            if evparams.walker.fixUnixWindows.doFix
              %new treatment. Tries to minimize sensibility to epsilon values
              epsilon                = evparams.walker.fixUnixWindows.epsilon;
              epsilons               = valuesAbs<epsilon;
              bothepsilons           = epsilons(1,:) & epsilons(2,:);
              toosmall               = toosmall & (~bothepsilons);
            end
            if any(toosmall & (mn==1))
              %if any connection is imbalanced because of too small p1i, add
              %p1i's neighbours
              extended(conns(p1i,:))             = true;
            end
            mins2                                = toosmall & (mn==2);
            if any(mins2)
              %if any connection is imbalanced because of too small p2i, add
              %p2i's neighbours
              if correctSingleRow
                extended(any(conns(p2n(mins2),:), 1)) = true;
              else
                extended(any(conns(p2n(mins2),:))) = true;
              end
              connsi(:,mins2)                      = false; %to avoid to examine again these points
            end
          end
        end
      end
      %finally, the points to be considered are:
      ps     = pos(extended,:);
      idxps   = find(extended);
      %now, let's get the points of the convex hull
      if d3
        if useHull
          points = convhulln(ps);
        else
          points = 1:numel(idxps);
        end
      else
        points = ConvHull2D(ps(:,1),ps(:,2));
      end
      if correctSign
        signos = sign(autovec(idxps(points)));
      else
        signos = sign(autovec(points));
      end
      if d3
        %This procedure is highly arbitrary, and is sure to yield only a
        %very limited number of configurations, but is necessary to use it,
        %as imperfect as it is, in order to canalize the configuration
        %of the ABC to be more prone to yield useful movements.
        
        %get all reasonable pairs
        if useHull
          pairs  = [points(:,[1 2]); points(:,[2 3]); points(:,[1 3])];
          signos = [signos(:,[1 2]); signos(:,[2 3]); signos(:,[1 3])];
          %pairs having distinct sign
          pairs  = pairs(signos(:,1)~=signos(:,2),:);
          pairs  = unique([min(pairs,[],2), max(pairs,[],2)], 'rows');
        else
          signos1 = signos>0;
          ps1     = find(signos1(:));
          ps2     = find(~signos1(:));
          %pairs having distinct sign
          pairs   = [repmat(ps1, numel(ps2), 1), reshape(repmat(ps2, 1, numel(ps1))', [], 1)];
        end
        
        validPairs = true(size(pairs,1),1);
        
        switch evparams.walker.d3.borderRating
          case 'method1' %to be used with d3.rotationMode=='inchworm'
            if usehoh
              error('border rating method not compatible with hand over hand configuration!!!!');
            end
        
            %get vector displacement for each pair
            vectorAB = ps(pairs(:,1),:)-ps(pairs(:,2),:);
            %for each pair, get the point closest to the X axis.
            [distclosest idxclosest] = min(realsqrt([sum(realpow(ps(pairs(:,1),[2 3]), 2), 2), ...
                                                     sum(realpow(ps(pairs(:,2),[2 3]), 2), 2)]), [], 2);
            closest2 = idxclosest==2;
            pairs(closest2,:) = pairs(closest2, [2 1]);
            %get the vector from the closest point to the X axis
            xvec     = [zeros(size(pairs,1),1), ps(pairs(:,1),[2 3])];
            %get the cross product of vectorAB and xvec
            crossxv  = [xvec(:,2).*vectorAB(:,3)-xvec(:,3).*vectorAB(:,2), ...
                        xvec(:,3).*vectorAB(:,1)-xvec(:,1).*vectorAB(:,3), ...
                        xvec(:,1).*vectorAB(:,2)-xvec(:,2).*vectorAB(:,1)];
            %get the absolute cosine between cross1   and the X axis
            cosxv    = abs(crossxv(:,1))./realsqrt(sum(crossxv.*crossxv, 2));
            %get the absolute cosine between vectorAB and the X axis
            cosAB    = abs(vectorAB(:,1))./realsqrt(sum(vectorAB.*vectorAB, 2));
            %scale the distance between vectorAB and the X axis, then use a
            %sigmoid function to rate it
            ratiox = distclosest*evparams.walker.d3.distRatio;
            ratiox = 1-1./(1+exp(-(ratiox-0.05)*500));
            
            %the closest to 0 is the rating, the better
            rating   = cosxv+cosAB.*cosAB+ratiox;
          case 'handoverhand1'
            if ~usehoh
              error('border rating method not compatible with a configuration not being hand over hand!!!!');
            end
            
            psOriented = hoh.orientedPos(extended,:);
            
            %get the vectors from the center of mass of the leg to each
            %point in the pair
            vectorA = bsxfun(@minus, psOriented(pairs(:,1),:), cm);
            vectorB = bsxfun(@minus, psOriented(pairs(:,2),:), cm);
            %project them to the XZ plane
            vectorA(:,2) = 0;
            vectorB(:,2) = 0;
            %make them unitary
            vectorA = vectorA/realsqrt(sum(vectorA.*vectorA));
            vectorB = vectorB/realsqrt(sum(vectorB.*vectorB));
            
            %get their dot product with the X axis
            cosAX   = vectorA(:,1);
            cosBX   = vectorB(:,1);
            
            %PLEASE NOTE: we can ask whether we should test their direction
            %in +X or -X. We will test always in the +X, and the user must
            %tweak these variables:
            %    -evparams.walker.d3.hoh.swingAngle
            %    -evparams.walker.leg.stateInit
            %    -evparams.walker.leg.stateInitTime
            %in order to place the ATP bonds in the best possible
            %orientation, TAKING INTO ACCOUNT that THE FIRST LEG IS THE ONE
            %SWUNG IN A POSITIVE ANGLE.
            
            %the closest to 0 is the rating, the better
            rating = 1-cosAX.*cosBX;
            
          otherwise
            error('option evparams.walker.d3.borderRating==%s not understood!!!!\n', any2str(evparams.walker.d3.borderRating));
        end
        
            %TO IMPLEMENT
        %determine the susceptible pairs
        pairsI = idxps(pairs);
        if numel(pairsI)==2
          pairsI = reshape(pairsI, 1, 2);
        end
        switch evparams.walker.d3.borderMode
          case 'onlyconnected'
            if ~correctPairIndex
              pairsI = pairs;
            end
            validPairs = validPairs & ismember([min(pairsI,[],2), max(pairsI,[],2)], [min(spEnds,[],2), max(spEnds,[],2)], 'rows');
          case 'anyone'
            %validPairs = true(size(pairs,1),1);%error('This option is currently being evaluated!!!!');
          case 'notconnected'
            validPairs = validPairs & (~ismember([min(pairsI,[],2), max(pairsI,[],2)], [min(spEnds,[],2), max(spEnds,[],2)], 'rows'));
          case 'dist2' %points not connected but with at least one common neighbour
            cnns   = any(conns(:,pairsI(:,1)) & conns(:,pairsI(:,2)))';
            validPairs = validPairs & cnns & (~ismember([min(pairsI,[],2), max(pairsI,[],2)], [min(spEnds,[],2), max(spEnds,[],2)], 'rows'));
          case 'dist2Near' %points not connected but with at least one common neighbour, and nearer than a given distance
            cnns   = any(conns(:,pairsI(:,1)) & conns(:,pairsI(:,2)))';
            validPairs = validPairs & cnns & (~ismember([min(pairsI,[],2), max(pairsI,[],2)], [min(spEnds,[],2), max(spEnds,[],2)], 'rows'));
            validPairs(validPairs) = evparams.walker.d3.borderDist>=realsqrt(sum(realpow(pos(pairsI(validPairs,1),:)-pos(pairsI(validPairs,2),:), 2), 2));
          otherwise
            error('option evparams.walker.d3.borderMode==%s not understood!!!!\n', any2str(evparams.walker.d3.borderMode));
        end
        if isfield(evparams.walker.d3, 'filterATPTooClose')
          atppos = (pos(pairsI(validPairs,1),:)+pos(pairsI(validPairs,2),:))/2;
          atpdists = distanceMatrix(pos, atppos);
          validPairs(validPairs) = all(atpdists>=evparams.walker.d3.filterATPTooClose)';
        end
        pairs = pairs(validPairs,:);
        if isempty(pairs)
          legs = 'Error while calculating legs: no suitable pair was found';
          return;
        end
        rating = rating(validPairs);
        [pair pair] = min(rating);
        pair = pairs(pair,:);
        switch evparams.walker.d3.ATPConnMode
          case 'simplest'
            %ATP bonds in the middle of the members
            legs(k).ABCweight  = [1; 1];
            legs(k).ABCbinding = reshape(idxps(pair), [], 1);
          case 'triangle'
            %a triangle is implemented: ATP bonds to both members, and also
            %to their nearest neighbour
            %leg(k).dynLengthFactor is customized in order to allow for a
            %fine control of the modification: the members' bonds are
            %shortened, but the neighbour's isn't
            error('still not implemented!!!');
          case 'neighsPassive'
            %ATP bonds to both members, and also to all their neighbours
            %leg(k).dynLengthFactor is customized in order to allow for a
            %fine control of the modification: the members' bonds are
            %shortened, but the neighbours' aren't
            pair                    = idxps(pair);
            bothneighs              = conns(:,pair(1)) & conns(:,pair(2));
            bothneighs(pair)        = false; %just in case...
            bothneighs              = find(bothneighs);
            legs(k).ABCweight       = [ones(2,1); zeros(numel(bothneighs), 1)];
            legs(k).ABCbinding      = [reshape(pair, [], 1); reshape(bothneighs, [], 1)];
            legs(k).dynLengthFactor = [repmat(evparams.walker.leg.dynLengthFactor, 2, 1); ones(numel(bothneighs), 1)];
          otherwise
            error('option evparams.walker.d3.ATPConnMode==%s not understood!!!!\n', any2str(evparams.walker.d3.ATPConnMode));
        end
      else
        %and now, get the borderline edges in the hull (the ones having a point
        %from each region)
        borders = find(signos(1:end-1)~=signos(2:end));
        borders = [reshape(points(borders), [], 1), reshape(points(borders+1), [], 1)];
        switch size(borders, 1)
          case 1
            %this is freaking strange, surely signals a very badly folded
            %structure, not suitable for our purposes.
            legs = 'Too badly folded leg';
            return
          case 0
            %even freakier than before
            legs = 'Even worse folded leg';
            return
        end
        %if there are several such borders, get the longest. THE IDEA IS THAT
        %PERTURBING THIS BORDER WILL HAVE THE MAXIMAL EFFECT IN THE HINGE
        lengths = realsqrt(sum(realpow(ps(borders(:,1),:)-ps(borders(:,2),:), 2), 2));
        [b b]   = max(lengths);
        legs(k).ABCweight  = [1; 1];
        legs(k).ABCbinding = reshape(idxps(borders(b,:)), [], 1);
      end
%       else
%         signos = sign(autovec(idxps(points)));
%         POR AQUI!!!!!!!
%         -LO QUE PARECE MÁS SENSATO PARA CANALIZACIÓN:
%         -REPOSICIONAR LA ESTRUCTURA. DE HECHO, LA HEURISTICA DEL ABC DEBERIA SER DEPENDIENTE DEL METODO DE ROTACION Y NO AL REVÉS
%         -PARA CADA PIERNA
%         -HACER EL CONVEX HULL DEL INTERFACE
%         -TOMAR PARES DE VERTICES (A,B) NO CONTECTADOS Y DE REGIONES DIFERENTES
%         -CONSIDERAR EL VECTOR X, que representa el microtubulo
%         -CONSIDERAR EL PUNTO Q, que representa el centro de masas de los "toes" (alternativamente, un punto de intersección entre X y la estructura)
%         -POR CONVENCION, SEA A EL VERTICE MÁS CERCANO ENTRE A y B a Q.
%         -CONSIDERAR EL VECTOR V=ABxAQ
%         -LA IDEONEIDAD DEL PAR (A,B) debería estar condicionada a
%         -EL VECTOR V forma un angulo lo más recto posible con X (esto puede usarse como una heurística de que el "power stroke" hará moverse la pierna en la dirección de X)
%         -EL VECTOR AB tiene una longitud aproximada dada
%         -LOS PUNTOS AB ESTÁN AMBOS CONECTADOS A UN PUNTO C, QUE FORMARÁ PARTE DEL ABC (la longitud del muelle ATP-C no cambiará, o bien incluso aumentará)
%         
%         -MÁS COSAS QUE HACER:
%         -HACER QUE EN PROCESO DE MUTACION LOS MUELLES DE LA CADENA TENGAN UN K SUPERIOR
%         -SI, TRAS UNA MUTACIÓN, UN ESLABON DE LA CADENA ESTÁ SEPARADO EN UNA CANTIDAD QUE NO ESTÁ EN LONG_CANONICA +/- TOLERANCIA, RECHAZAR LA MUTACION
%         -LA PROBABILIDAD DE UNA MUTACIÓN ES PROPORCIONAL A LA LONGITUD DEL MUELLE
%         -...
%         -EVIDENTEMENTE, HAY QUE HACER EL CODIGO DE MUTACION, QUE NO ES TRIVIAL (NO SE DEBE MUTAR LA CADENA, POR EJEMPLO)
%         -COMPROBAR LA CANTIDAD DE MODOS DENTRO DE LA PRIMERA DECADA DEL GRAFICO DEL PAPER DEL PNAS. SE PUEDE HACER CON UN KDE
%         
%         -COSAS SOBRE LAS QUE HABLAR:
%         -DINÁMICA DE LA EVOLUCIÓN
%         -HECHO DE QUE SE HA EVOLUCIONADO
%         -ESTUDIO DEL SPECTRAL GAP
%       end
    end
    
    if usehoh
      if strcmp(evparams.walker.d3.hoh.fusionMode, 'justOnePoint')
        pivot = find(legs(idxlegs).ABCbinding==hoh.pivot);
        %make sure that the pivot is not present in the machinery: since it
        %is not replicated, it might cause problems
        if ~isempty(pivot)
          if isfield(legs(idxlegs), 'dynLengthFactor') && ( numel(legs(idxlegs).dynLengthFactor)==numel(legs(idxlegs).ABCbinding(pivot)) )
            legs(idxlegs).dynLengthFactor(pivot) = [];
          end
          legs(idxlegs).ABCweight(pivot)         = [];
          legs(idxlegs).ABCbinding(pivot)        = [];
        end
      end

      %the second leg will be the one swung in a negative angle
      otherleg = mod(idxlegs,2)+1;
      legs(otherleg) = legs(idxlegs);
      %the binding sites are displaced for the second leg
      legs(otherleg).ABCbinding = hoh.newIndexes(legs(otherleg).ABCbinding + hoh.idxoffset);
      legs(otherleg).toes       = hoh.newIndexes(legs(otherleg).toes       + hoh.idxoffset);
    end
    
  case 'meanpos'
    if d3
      error('This mode to calculate the ABC is not supported in 3D!!!! (In fact, it is entirely obsolete)');
    end
    for k=1:2
      %all points in the interface and their neighbours
      points = [interfaces(k).points1(:); interfaces(k).points2(:)];
      allCandidates = any(conns(points,:), 1)';
      allCandidates(points) = true;
      allCandidatesI = find(allCandidates);
      %absolute values of the components of the eigenvector
      vecabs  = 1./autovecabs;
      vecabsi = vecabs(points);
      %position of the ABC
      posABC  = sum(bsxfun(@times, pos(points,:), vecabsi))/sum(vecabsi);
      %distance from the ABC to candidate points
      dists   = realsqrt(sum(realpow(bsxfun(@minus, pos(allCandidates,:), posABC), 2), 2));
      %nearest point
      mindist = min(dists);
      %points relatively near
      near    = allCandidatesI(dists<=(mindist*evparams.walker.cutoff));
      if numel(near)<4
        %too few near points: extend to their neighbours inside the
        %"allCandidates" set
        binding = (any(conns(near,:))') & allCandidates;
        binding(near) = true;
        binding = find(binding);
      else
        binding = near;
      end
      legs(k).ABCweight  = vecabs(binding);
      legs(k).ABCbinding = binding;
    end
  otherwise
    error('ABCMode %s not understood!!!!!', evparams.walker.makeABCMode);
end
  
function toes = calculateToes(autovec, conns, interfaces, regions, evparams)
%now, let's calculate toes

autovecabs = abs(autovec);

%first, get each region's identifier (there must be 3 regions)
% r1 = [interfaces.region1];
% r2 = [interfaces.region2];
r1 = [interfaces(1).region1 interfaces(1).region2];
r2 = [interfaces(2).region1 interfaces(2).region2];

%detect the region having two interfaces and filter it: we want distal
%regions, having just one interface
if any(r1(1)==r2)
  r2 = r2(r2~=r1(1));
  r1 = r1(2);
else
  r2 = r2(r2~=r1(2));
  r1 = r1(1);
end
rs = [r1; r2];

toes = {[];[]};

%for each leg
for k=1:2
  %vertices of the region
  vs = find(regions==rs(k));
  
  if ~isempty(evparams.walker.minAmountOfToes)
    %get the minimal amount of toes
    minToes = evparams.walker.minAmountOfToes(conns, vs);
    %sort values in the region
    [sorted idx] = sort(autovecabs(vs), 'descend');
    %if there are other vertices with values nearly identical to the lowest
    %values among toes, also include them
    nToes = minToes;
    tol = sorted(minToes)*0.01; %eps(100);
    while abs(sorted(minToes)-sorted(nToes+1))<tol
      nToes = nToes+1;
    end
    toes{k} = vs(idx(1:nToes));
  else
    %this is to make toes all the leg's vertices. This way, we make sure it
    %can adhere to the microtubule
    toes{k} = vs;
  end
end


