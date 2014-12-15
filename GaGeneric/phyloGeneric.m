function pt = phyloGeneric(getBasicInfFun, printIndividualFun, poph, restrictgens, trimIdenticalToParent)
% PHYLOTREE A heavily modified form of TRIMTREEPLOT. 
%inputs:
%   -poph: population history, loaded with loadSimulation
%   -restrictgens: If absent or empty, the whole tree will be rendered.
%                  If numel(restrictgens)==1, the tree will be rendered for
%                  generations 0..restrictgens
%                  If numel(restrictgens)==2, the tree will be rendered for
%                  generations restrictgens(1)..restrictgens(2)
%   -trimIdenticalToParent: remove nodes that are identical to their
%                           parents. This may improve readability and 
%                           trazability if the tree is very big
%
%outputs:
%   -pt: struct. It contains the following fields:
%      +highlight, that is a function to highlight nodes. The function
%      takes these arguments: 
%       *mode: +If 'triplet', the second argument must be a Nx3 array of N
%               triplets [generation rangeid individual] to highlight
%              +If 'indexesTree', the second argument must be a list of
%               identifiers of nodes (as in poph.idx)
%              +If 'tracking', there are three more arguments:
%                   *individuals
%                   *startGen
%                   *endGen
%               where the individuals given by 'individuals' from
%               generation 'startGen' are tracked until generation 'endGen'
%               (please note that endGen can be either lower or greater
%               than startGen) 
%        After these first arguments, pairs (option, value) can be specified:
%           *option 'color': color specification for hihglighted nodes
%                            (default 'k')
%           *option 'removePrevious': flag to remove previous highlights
%                                     before drawing this one
%           *option 'plotTracking': if mode == 'tracking', this flag tells
%                                   whether to plot the amount of tracked
%                                   individuals by generation 
%      +highlightByAttr: 
%          Highlights individuals by attributes. The following arguments
%          are used to decide what individuals to highlight:
%            -attrName: a cell list of names of attributes of individuals
%                      (for example, {'fitness', 'nleafs'}).
%            -transformFun: if empty, this does nothing. If it is a cell
%                           array the same size as attrName, each position
%                           must be empty, or a function handle or a
%                           function name. If not empty, the i-th vector of 
%                           attribute values (corresponding to the i-th
%                           attrName) will be transformed by the 
%                           corresponding i-th function in transformFun,
%                           before being compared to the i-th attrRange.
%            -attrRange: a cell list (the same size as attrName) specifying
%                        the range for each attribute. If the range has 1
%                        element, only individuals with that precise
%                        attribute will be selected. If it has 2 elements,
%                        these define the lower and upper limits of the
%                        range of selected individuals
%            -Optional 4th argument: min and max generation to highlight
%            -Other optional arguments: pairs 'option', value:
%                                 -'removePrevious': flag to say if the
%                                  previous highlight should be removed
%                                 -'color': this is to say the ColorSpec to
%                                  use in the highlight
%      +highlightOnePerGeneration: 
%          Highlights an individual per generation, based on the value of
%          some attribute. The following arguments are used to decide what
%          individuals to highlight: 
%            -attrName: The name of an attribute (for example, 'fitness')
%            -transformFun: if empty, this does nothing. If a function
%                           handle or a function name, the vector of
%                           attribute values is transformed by this
%                           function before being passed to selecFun
%            -selecFun: A function to select an individual from a
%                       list of individuals (for example @min or @max; if
%                       the function is not any of these two, it
%                       must return as its first output the index of the
%                       selected individual. 
%                                  Optional argument: min and max 
%                                  generation to highlight
%                        Other arguments: pairs 'option', value:
%                                 -'removePrevious': flag to say if the
%                                  previous highlight should be removed
%                                 -'color': this is to say the ColorSpec to
%                                  use in the highlight
%            -Optional 4th argument: min and max generation to highlight
%            -Other optional arguments: pairs 'option', value:
%                                 -'removePrevious': flag to say if the
%                                  previous highlight should be removed
%                                 -'color': this is to say the ColorSpec to
%                                  use in the highlight
%      +switchShowSimmilarity: this is a function that accepts a flag. The
%           flag determines whether the trees will be shown when computing
%           simmilarity.
%      +setTreeRandomDevelopment: it accepts an argument that may have
%           three values: 
%               -true:         the trees are painted with randomized
%                              development 
%               -false:        the trees are painted with non-randomized
%                              development 
%               -empty matrix: the trees are painted using the simulation
%                              parameters for randomization
%             
%
% BEWARE: SOMETIMES, THE ALGORITHM TREELAYOUT ASSIGNS THE SAME POSITION IN
%         THE X AXIS TO TWO INDIVIDUALS IN THE SAME GENERATION, IF ONE
%         IS THE PARENT OF THE OTHER. TO PREVENT THEM FROM OVERLAPPING; WE
%         PLOT THE SON DISPLACED DOWNWARDS (ADDING 1/3 TO ITS Y POSITION,
%         ALTHOUGH THAT DOESN'T MEAN THAT ITS GENERATION HAS A FRACTIONAL
%         COMPONENT; OF COURSE!!!)

idxD = poph.tree.idxD;
idxA = poph.tree.idxA;
gD   = poph.tree.gD;

mingen = min(gD);
maxgen = max(gD);

% sparseIndexesPoph = sparse(1, poph.idx, 1:numel(poph.idx), 1, max(idxD));
% sparseIndexesTree = sparse(1, idxD,     1:numel(idxD),     1, max(idxD));

%remove nodes that are not within the specified span of generations
if nargin>3
  toUse = [];
  switch numel(restrictgens)
    case 1
      if restrictgens<maxgen
        toUse  = gD<=restrictgens;
        maxgen = restrictgens;
      end
    case 2
      if restrictgens(1)>restrictgens(2)
        error('Are you crazy? The second parameter is [startGen endGen], so keep startGen lower or equal than endGen, stupid!');
      end
      if (restrictgens(1)>mingen) || (restrictgens(2)<maxgen)
        toUse  = (gD>=restrictgens(1)) & (gD<=restrictgens(2));
        maxgen = min(restrictgens(2), maxgen);
        mingen = max(restrictgens(1), mingen);
      end
    case 0
      %do nothing
    otherwise
      error('Second argument must have either 1 or 2 elements, not %d!!!\n', numel(restrictgens));
  end
  if ~isempty(toUse)
    idxD   = idxD(toUse);
    idxA   = idxA(toUse);
    gD     = gD(toUse);
  end
end

%remove nodes that are identical to their parents
if (nargin>4) && trimIdenticalToParent
  chg = poph.tree.change;
  if ~isempty(toUse)
    chg = chg(toUse);
  end
  %find individuals not changing with respect to their parents
  nochange = cellfun(@(x)isempty(x) || ((numel(x)==1) && (x=='=')), chg);
  clear chg;
  %make sure that the individuals are not roots!!!
  nochange = find(nochange & (idxA~=0));
  for k=1:numel(nochange)
    toDelete    = idxD(nochange(k));
    itsAncestor = idxA(nochange(k));
    idxA(idxA==toDelete) = itsAncestor;
  end
  idxD(nochange) = [];
  idxA(nochange) = [];
    gD(nochange) = [];
  clear nochange;
end

clear toUse;

%make index matrix of parents
if (numel(poph.tree.idxD)==numel(idxD)) && isfield(poph.tree, 'parents')
  parents = poph.tree.parents;
else
%   try
%     %get parent indexes
%     sparseIndexes   = sparse(1,idxD+1,1:numel(idxD));
%     parents         = full(sparseIndexes(idxA+1)); 
%     clear sparseIndexes;
%   catch ME
%     clear sparseIndexes parents;
%     if strcmpi(ME.identifier, 'MATLAB:nomem')
%       %that sparse matrix was too much, poor old MATLAB! Let's do it the
%       %memory-safe way
      [nevermind, parents] = ismember(idxA, idxD); %#ok<ASGLU>
      parents = parents';
      clear nevermind;
%     else
%       rethrow(ME);
%     end
%   end
  if numel(poph.tree.idxD)==numel(idxD)
    poph.tree.parents = parents;
  end
end

useOneTreeRoot = true;

%find roots. If we set this to the empty matrix, individuals without
%offspring in the first generation will not be shown
roots             = find(parents==0);%[];%find(parents==0);

if useOneTreeRoot
  numTreeRoots  = 1;
else
  numTreeRoots  = numel(roots);
end
%setting super-roots' y coordinates to -INF, their associated lines are not
%drawn, but the roots are still shown
y               = [repmat(-inf, 1, numTreeRoots), gD'];%[gD(roots)'-1, gD'];
clear gD;
%make new roots for them
parents         = [zeros(1,numTreeRoots), parents+numTreeRoots];
%rewire old roots
roots           = roots+numTreeRoots;
parents(roots)  = 1:numTreeRoots;


x=trimtreelayout(parents);

%this bit of code is to disentangle nodes overlapped with their parents
pos  = [x((numTreeRoots+1):numel(parents))', y((numTreeRoots+1):numel(parents))'];
np   = parents((numTreeRoots+1):numel(parents))-numTreeRoots;
np(roots-numTreeRoots) = 0;
posp = repmat(inf, size(pos));
np0  = np>0;
posp(np0,:) = pos(np(np0),:);
overlapped  = all(posp==pos, 2)';
if any(overlapped)
  warning('phylotree:overlappedNodes', ...
          'PLEASE NOTE: Some nodes are overlapped with their parents (the cause is, likely, that one is the parent of the other, but they are nevertheless in the same generation). To fix this, the overlapping sons will be displaced 1/3 downwards from their parents, and also they will be given a very tiny shift to the right.'...
         );
  %do correction, displacing sons downwards in respect of their parents.
  y((numTreeRoots+1):numel(parents)) = y((numTreeRoots+1):numel(parents))+overlapped/3;%*min(1/size(pos,1)*5, 1/3);
  %we also displace it very gently in the x direction
  x((numTreeRoots+1):numel(parents)) = x((numTreeRoots+1):numel(parents))+overlapped*(1/size(pos,1)/3);
  clear pos posp overlapped;
  pos  = [x((numTreeRoots+1):numel(parents))', y((numTreeRoots+1):numel(parents))'];
  np   = parents((numTreeRoots+1):numel(parents))-numTreeRoots;
  np(roots-numTreeRoots) = 0;
  posp = repmat(inf, size(pos));
  np0  = np>0;
  posp(np0,:) = pos(np(np0),:);
  overlapped  = all(posp==pos, 2);
  if any(overlapped)
    warning('phylotree:overlappedNodes', ...
            'BEWARE: we have tried our best to prevent nodes to be placed in the same position as their parents, but, nevertheless, there remain some of them overlapped with their parents. So, please beware that the displayed tree is not entirely accurate, as some nodes will be hidden below others.'...
           );
  end
  clear pos posp overlapped;
end
clear np np0;

if numel(roots) == 0
  f = find(parents~=0);
  pp = parents(f);
  X = [x(f); x(pp); repmat(NaN,size(f))];
  Y = [y(f); y(pp); repmat(NaN,size(f))];
  X = X(:);
  Y = Y(:);
else
  %f = (numel(roots)+1):numel(parents);
  pp = parents((numTreeRoots+1):numel(parents));
  X = [x((numTreeRoots+1):numel(parents)); x(pp); repmat(NaN,1,numel(parents)-numTreeRoots)];
  Y = [y((numTreeRoots+1):numel(parents)); y(pp); repmat(NaN,1,numel(parents)-numTreeRoots)];
  X = X(:);
  Y = Y(:);
end

sepy    = 0.2;
minsepy = 0.01;

fig = figure;
ax = axes('XLim', [-0.05 1.05], 'YLim', [(mingen-1) (maxgen+1)], ...
          'YDir', 'reverse', 'YGrid', 'on', 'XTick', [], 'XAxisLocation', 'top', ...
          'Units', 'normalized', 'Position', [0.15 (sepy+minsepy) 0.8 (1-(sepy+2*minsepy))]);

%line(X, Y, 'Color', 'r', 'LineStyle', '-');
hold on;
h2 = line(X,...%x((numel(roots)+1):numel(parents)), ...
          Y,'Color', 'r',...%y((numel(roots)+1):numel(parents)), ...
          'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');%, 'MarkerSize', 10);

set(fig,'Toolbar','figure');
uinfo = uicontrol('Style', 'text', 'Units', 'normalized', ...
                  'FontName', 'FixedWidth', 'Position', [0 0 1 sepy], ...
                  'HorizontalAlignment', 'left', 'UserData', []);

edSelected = uicontrol('Style', 'edit', 'UserData', [], ...
                       'Units', 'normalized', 'Interruptible', 'off', 'BusyAction', 'cancel', ...
                       'Position', [0 0.93 0.1 0.07], 'String', sprintf('[gen rind ind]'), ...
                       'TooltipString', sprintf('Put here a [gen rangeid index] triple\nfor example, [12 1 87]\n(default rangeid is 1)'));
% edSimil    = uicontrol('Style', 'edit', 'UserData', [], ...
%                        'Units', 'normalized', 'Interruptible', 'off', 'BusyAction', 'cancel', ...
%                        'Position', [0 0.86 0.1 0.07], 'String', sprintf('[gen rind ind]'), ...
%                        'TooltipString', sprintf('Put here a [g1 r1 i1; g2 i2 r3] double triple\nfor example, [12 1 87; 10 1 23]\n(default rangeid is 1)'), ...
%                        'UserData', true);
closure           = struct('tree', [], 'widgets', [], 'getBasicInfFun', getBasicInfFun, 'printIndividualFun', printIndividualFun);
closure.data      = struct('x', [], 'y', [], 'genSpan', maxgen-mingen, ...
                           'minGen', mingen, 'maxGen', maxgen, ...
                           'poph', poph, 'idxD', idxD, 'idxA', idxA);
closure.widgets   = struct('ax', ax, 'tree', h2, 'uinfo', uinfo, ...
                           'edSelected', edSelected, ... %'edSimil', edSimil, ...
                           'fig', fig);
closure.data.x    = x((numTreeRoots+1):numel(parents)); clear x;
closure.data.y    = y((numTreeRoots+1):numel(parents)); clear y;
set(h2,          'ButtonDownFcn', {@showInf,          closure}, 'Interruptible', 'off');
set(edSelected,  'Callback',      {@makeHighLighting, closure}, 'Interruptible', 'off');
% set(edSimil,     'Callback',      {@showSimilarity,   closure}, 'Interruptible', 'off');
hold off;
ylabel(['generations = ' int2str(mingen) '...' int2str(maxgen)]);

pt = struct(...
  'highLight', highLightFUN(closure), ...
  'highLightOnePerGeneration', highLightOnePerGenerationFUN(closure), ...
  'highLightByAttr', highLightByAttrFUN(closure), ...
  'classifyByAtrr', classifyByAtrrFUN(closure), ...
  ...%'switchShowSimmilarity', showSimmilarityFUN(closure), ...
  'setRandomDevelopment', setTreeRandomDevelopmentFUN(closure) ...
  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minGen maxGen, otherArgs] = getGenerationRange(closure, otherArgs) %#ok<INUSL>
if (numel(otherArgs)==0) || ischar(otherArgs{1})
  minGen = closure.data.minGen;
  maxGen = closure.data.maxGen;
else
  genRange = otherArgs{1};
  switch numel(genRange)
    case 1
      genRange = [closure.data.minGen genRange];
    case 2
    otherwise
      error('genRange must have either 1 or 2 elements!!!');
  end
  minGen = max(closure.data.minGen, genRange(1));
  maxGen = min(closure.data.maxGen, genRange(2));
  otherArgs = otherArgs(2:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = highLightOnePerGenerationFUN(closure) %#ok<INUSL>
fun = @(varargin) doHighLightOnePerGeneration(closure, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectedINDEXES = doHighLightOnePerGeneration(closure, attrname, transformFun, selecFun, varargin) %#ok<INUSL>
otherArgs = varargin;
if (~isfield(closure.data.poph, attrname)) || (size(closure.data.poph.(attrname),1)~=size(closure.data.poph.fitness,1))
  fprintf('Sorry, but the record for this simulation does not hold information for <%s>.\nI cannot compute the highlighting', attrname);
  return;
end
[minGen maxGen otherArgs] = getGenerationRange(closure, otherArgs);
selected = zeros(maxGen-minGen+1,1);
g=minGen;
use2outputs = (ischar(selecFun)                 && any(strcmp(selecFun,           {'max', 'min'}))) || ...
              (isa(selecFun, 'function_handle') && any(strcmp(func2str(selecFun), {'max', 'min'})));
for k=1:numel(selected)
  thisGen = find(closure.data.poph.generation==g);
  toBeEvaluated = closure.data.poph.(attrname)(thisGen);
  if ~isempty(transformFun)
    toBeEvaluated = feval(transformFun, toBeEvaluated);
  end
  if use2outputs
    [nevermind toSelect] = feval(selecFun, toBeEvaluated); %#ok<ASGLU>
  else
               toSelect  = feval(selecFun, toBeEvaluated);
  end
  clear toBeEvaluated;
  if ~isempty(toSelect)
    selected(k) = thisGen(toSelect);
  end
  g = g+1;
end
selected = selected(selected>0);
selectedINDEXES = selected;
selected = closure.data.poph.idx(selected);
if ~isempty(selected)
  doHighLight(closure, 'indexestree', selected, otherArgs{:});
else
  fprintf('Sorry, no individuals were selected, so I will do nothing.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = highLightByAttrFUN(closure) %#ok<INUSL>
fun = @(varargin) doHighLightByAttr(closure, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toSelectINDEXES = doHighLightByAttr(closure, attrname, transformFun, attrRange, varargin) %#ok<INUSL>
otherArgs = varargin;
[minGen maxGen otherArgs] = getGenerationRange(closure, otherArgs);
toSelect = (closure.data.poph.generation>=minGen) & (closure.data.poph.generation<=maxGen);
for k=1:numel(attrname)
  if (~isfield(closure.data.poph, attrname{k})) || (size(closure.data.poph.(attrname{k}),1)~=size(closure.data.poph.fitness,1))
    fprintf('Sorry, but the record for this simulation does not hold information for <%s>.\nI cannot compute the highlighting', attrname{k});
    return;
  end
  toBeEvaluated = closure.data.poph.(attrname{k})(toSelect);
  if (~isempty(transformFun)) && (~isempty(transformFun{k}))
    toBeEvaluated = feval(transformFun{k}, toBeEvaluated);
  end
  switch numel(attrRange{k})
    case 1
      toSelect(toSelect) =  toBeEvaluated== attrRange{k};
    case 2
      toSelect(toSelect) = (toBeEvaluated >= attrRange{k}(1)) & ...
                           (toBeEvaluated <= attrRange{k}(2));
    otherwise
        error('attrRange must have either 1 or 2 elements in position %d!!!', k);
  end
  
end

toSelectINDEXES = find(toSelect);
toSelect = closure.data.poph.idx(toSelect);

if ~isempty(toSelect)
  doHighLight(closure, 'indexestree', toSelect, otherArgs{:});
else
  fprintf('Sorry, no individuals were selected, so I will do nothing.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = classifyByAtrrFUN(closure)
fun = @(varargin) doClassifyByAttr(closure, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doClassifyByAttr(closure, attrname, limits, colors, colornan)
limits = [-inf; limits(:); inf];
limits = [limits(1:end-1), limits(2:end)];
if size(limits,1)~=size(colors,1)
  error('The number of sets is %d, but the number of colors is %d!!!', size(limits,1), size(colors,1));
end
if exist('colornan', 'var')
  doHighLightByAttr(closure, {attrname}, {@isnan}, {true}, 'removeprevious', false, 'color', colornan);
end
for k=1:size(limits, 1)
  doHighLightByAttr(closure, {attrname}, [], {limits(k,:)}, 'removeprevious', false, 'color', colors(k,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = highLightFUN(closure) %#ok<INUSL>
%makes a closure to execute highlights
fun = @(varargin) doHighLight(closure, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indexesALL = doHighLight(closure, varargin)
% The function takes these arguments in varargin:
%     *mode: +If 'triplet', the second argument must be a Nx3 array of N
%             triplets [generation rangeid individual] to highlight
%            +If 'indexesTree', the second argument must be a list of
%             identifiers of nodes (as in poph.idx)
%            +If 'tracking', there are three more arguments:
%                 *individuals
%                 *startGen
%                 *endGen
%             where the individuals given by 'individuals' from generation
%             'startGen' are tracked until generation 'endGen' (please note
%             that endGen can be either lower or greater than startGen)
%      After these first arguments, pairs (option, value) can be specified:
%         *option 'color': color specification for hihglighted nodes
%                          (default 'k')
%         *option 'removePrevious': flag to remove previous highlights
%                                   before drawing this one
%         *option 'plotTracking': if mode == 'tracking', this flag tells
%                                 whether to plot the amount of tracked
%                                 individuals by generation 

idxv = 1;
mode = lower(varargin{idxv}); idxv=idxv+1;

doTracking      = false;
switch mode
  case 'triplets'
    %hand over the triplets to the callback
    triplets    = varargin{idxv}; idxv=idxv+1;
    toHighLight = {mode triplets [] [] []};
  case 'indexestree'
    indexes     = varargin{idxv}; idxv=idxv+1;
    toHighLight = {'idxs' indexes [] [] []};
  case 'indexespoph'
    indexes     = varargin{idxv}; idxv=idxv+1;
    toHighLight = {'poph' indexes [] [] []};
  case 'tracking'
    individuals = varargin{idxv}; idxv=idxv+1;
    startGen    = varargin{idxv}; idxv=idxv+1;
    endGen      = varargin{idxv}; idxv=idxv+1;
    doTracking  = true;
  otherwise
    error('Cannot understand mode %s!!!!!', any2str(mode));
end

color                 = 'k';
removePrevious        = true;
plotTracking          = false;
showFirst             = true;

while idxv<=numel(varargin)
  option              = varargin{idxv}; idxv=idxv+1;
  if idxv>numel(varargin)
    error('Option %s has no value!!!', any2str(option));
  end
  switch lower(option)
    case 'showfirst'
      showFirst       = varargin{idxv}; idxv=idxv+1;
    case 'color'
      color           = varargin{idxv}; idxv=idxv+1;
    case 'removeprevious'
      removePrevious  = varargin{idxv}; idxv=idxv+1;
    case 'plottracking'
      plotTracking    = varargin{idxv}; idxv=idxv+1;
    otherwise
      error('I cannot understand option %s!!!', any2str(option));
  end
end

if doTracking
  if plotTracking
    plotTracking = {'plot'};
  else
    plotTracking = {};
  end
  %track individuals in the history, and hand the resulting indexes over
  %to the callback
  indexesALL          = trackPopulation(closure.data.poph, individuals, startGen, endGen, 'all', plotTracking{:});
  toHighLight         = {'idxs' closure.data.poph.tree.idxD(indexesALL) [] [] []};
end

toHighLight{3}        = color;
toHighLight{4}        = removePrevious;
toHighLight{5}        = showFirst;

makeHighLighting(nan,toHighLight,closure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeHighLighting(src,event,closure)
%this function highlights a node or several nodes

%get current point
if isnumeric(src) && (numel(src)==1) && isnan(src)
  %this is a hack to invoke this function programmatically
  usingTriplet   = false;
  calculateIndex = true;
  switch event{1}
    case 'triplets'
      usingTriplet = true;
      triplet = event{2};
    case 'idxs'
      idx     = event{2};
    case 'poph'
      calculateIndex = false;
      index     = event{2};
    otherwise
      error('Third option %s not recognized!!!!', any2str(event{3}));
  end
  color          = event{3};
  removePrevious = event{4};
  showFirst      = event{5};
  programmatic   = true;
else
  %read triplet
  calculateIndex = true;
  triplet        = get(closure.widgets.edSelected, 'String');
  usingTriplet   = true;
  removePrevious = true;
  color          = 'k';
  showFirst      = true;
  programmatic   = false;
end

if usingTriplet
  idx = anytriplet2idx(triplet, closure);
  if any(isnan(idx));   return;                     end
  if programmatic;      triplet = any2str(triplet); end
  errorStr = sprintf('\n    Sorry, but no node in this dataset is identified by the triplet %s (idx=%%s)', triplet);
else
  errorStr =         '\n    Sorry, but no node in this dataset is identified by the index %s ';
end

if calculateIndex
  index = idx2index(idx, closure, closure.data.idxD, errorStr);
end
if any(isnan(index))
  return;
end

obj = get(closure.widgets.edSelected, 'UserData');
if removePrevious && (~isempty(obj))
  %delete previous highlight, if it exists
  try
    cellfun(@delete, obj);
  catch %#ok<CTCH>
    for k=1:numel(obj)
      try
        delete(obj);
      catch %#ok<CTCH>
      end
    end
  end
  set(closure.widgets.edSelected, 'UserData', []);
end

%paint it
figure(closure.widgets.fig);
if showFirst
  showInf(nan, [closure.data.x(index(1)) closure.data.y(index(1))], closure);
end
objnew = line(closure.data.x(index), closure.data.y(index), ...
              'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', ...
              'ButtonDownFcn', {@showInf, closure});
if removePrevious || isempty(obj)
  set(closure.widgets.edSelected, 'UserData', {objnew});
else
  set(closure.widgets.edSelected, 'UserData', [obj {objnew}]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make a closure to change the TreeRandomDevelopment flag for @showInf
function fun = setTreeRandomDevelopmentFUN(closure)
%makes a closure
fun = @(newflag) changeFlagTreeRandomDevelopment(closure, newflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%change the TreeRandomDevelopment flag for @changeFlagShowSimmilarity
function changeFlagTreeRandomDevelopment(closure, newflag)
if isempty(newflag) || (islogical(newflag) && (numel(newflag)==1))
  set(closure.widgets.uinfo, 'UserData', newflag);
else
  fprintf('Sorry, but I cannot understand this parameter: %s\n', any2str(newflag));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showInf(src,event,closure)
%This function shows the individual's information in the figure

%get current point
if isnumeric(src) && (numel(src)==1) && isnan(src)
  %this is a hack to invoke this function programmatically
  currentPoint= event;
  plotTree    = false;
else
  %this means that the function has been called by clicking on the tree
  currentPoint= get(closure.widgets.ax, 'CurrentPoint'); 
  plotTree    = true;
end
cx            = currentPoint(1,1);
cy            = currentPoint(1,2);
%find nearest individual
errs          = abs([closure.data.x'-cx, closure.data.y'-cy]);
%we divide error in Y axis by the generation span, so both axes are
%normalized. This prevents some very strange effects when clicking
[ignore idx]  = min(errs(:,1)+errs(:,2)/closure.data.genSpan);%#ok<ASGLU>
%display nearest individual's info
% if errs(idx,2)<=(1/3)
  %idxphD: index of selected individual in poph fields
  idxphD      = find(closure.data.poph.idx==closure.data.idxD(idx));
  %idxtrD: index of selected individual in poph.tree fields
  idxtrD      = find(closure.data.poph.tree.idxD==closure.data.idxD(idx));
  useGenome   = isfield(closure.data.poph, 'genome');
  if closure.data.idxA(idx)~=0 %identifier of the selected individual's ancestor
    %idxphA: index of selected individual's ancestor in poph fields
    idxphA    = find(closure.data.poph.idx==closure.data.idxA(idx));
    %ancestor's info
    if useGenome
      ancestor = sprintf('\n%s', closure.data.poph.genome{idxphA});
    else
      ancestor = '';
    end
    ancestor  = sprintf(['CHANGE WITH RESPECT TO ANCESTOR: %s\nANCESTOR:\n  %s' ancestor], ...
                 closure.data.poph.tree.change{idxtrD}, ...
                 closure.getBasicInfFun(closure.data.poph, idxphA) ...
                 );
  else
    ancestor  = sprintf('\n\n  THIS INDIVIDUAL HAS NO ANCESTOR!!!!');
  end
  %selected individual's info
  if useGenome
    str       = sprintf('%s\n', closure.data.poph.genome{idxphD});
  else
    str       = '';
  end
  str         = sprintf(['SELECTED:\n  %s\n' str '%s'], ...
                 closure.getBasicInfFun(closure.data.poph, idxphD), ...
                 ancestor);
   if plotTree
     if isempty(idxphD)
       str = sprintf('\n            I cannot paint this individual, it does not exist in the population!!!\nLikely, this is an individual from a the last generation in a simulation that was interrupted before properly finishing');
     else
       sel_typ = get(closure.widgets.fig,'SelectionType');
       switch sel_typ 
         %case 'normal'; case 'extend'
         case 'alt'
           %disp('User did a control-click')
           set(closure.widgets.uinfo, 'String', sprintf('\n     Please wait while painting the tree...'));
           modeRandom = get(closure.widgets.uinfo, 'UserData');
           if closure.data.poph.params.evaluateAlways
             idxFile = idxphD;
           else
             idxFile = getFirstUnchangedAncestor(closure.data.poph, idxphD);
           end
           closure.printIndividualFun(closure.data.poph, closure.getBasicInfFun, idxFile, idxphD);
% %            fname = [closure.data.poph.basedir filesep sprintf('i_%04d_%d', closure.data.poph.generation(idxFile), closure.data.poph.rangeid(idxFile)) filesep sprintf('ind%04d.png', closure.data.poph.individual(idxFile))];
% %            if exist(fname, 'file')
% %              figure; imshow(fname);
% %              title(closure.data.poph.genome{idxphD});
% %              xlabel(getBasicInf(closure.data.poph, idxphD));
% %            end
%            fname = getSnapshot(closure.data.poph, idxphD);
%            if ~isempty(fname)
%              figure; imshow(fname);
%              title(closure.data.poph.genome{idxphD});
%              xlabel(getBasicInf(closure.data.poph, idxphD));
%            end
        end
     end
   end
set(closure.widgets.uinfo, 'String', str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fname = getSnapshot(poph, ind)
fname = [poph.basedir filesep sprintf('i_%04d_%d', poph.generation(ind), poph.rangeid(ind)) filesep sprintf('ind%04d.png', poph.individual(ind))];
if ~exist(fname, 'file')
  if poph.tree.parents(ind)>0
    fname = getSnapshot(poph, poph.tree.parents(ind));
  else
    fname = '';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = anytriplet2idx(triplet, closure)
%translate from a possible text form of triplets to node identifiers
argchar   = ischar(triplet);
if argchar 
  str     = triplet;
  triplet = str2num(triplet); %#ok<ST2NM>
end
if any(any(isnan(triplet))) || ( (size(triplet,2)<2) || ((size(triplet,2)>3)) )
  if ~argchar
    str   = any2str(triplet);
  end
  set(closure.widgets.uinfo, 'String', sprintf('\n    Sorry, but I cannot understand ''%s'' as a triplet [generation, rangeid, individual] (if you omit rangeid, the default will be 1, wich might be ok for most applications)', str));
  idx = nan;
  return
end
if size(triplet,2)==2
  %put 1 as default range
  triplet = [triplet(:,1), ones(size(triplet,1),1), triplet(:,2)];
end
%find node corresponding to triplet
idx   = closure.data.poph.triplet2idx(triplet(:,1), triplet(:,2), triplet(:,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = idx2index(idx, closure, dataset, errorStr)
%translate from node identifiers to indexes

if numel(idx)>1
  index = find(ismember(dataset, idx));
else
  index = find(dataset==idx,1);
end
if isempty(index)
  set(closure.widgets.uinfo, 'String', sprintf(errorStr, any2str(idx)));
  index = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = trimtreelayout(parent,post)
%   TRIMTREELAYOUT A modified standard function TREELAYOUT. 
%     Produces no heights for the leaves

n = length(parent);

 pv = [];
 if (size(parent,1)>1), parent = parent(:)'; end
 if (nargin<2) && (~all(parent==0 | parent>(1:n)))
     % This does not appear to be in the form generated by ETREE.
     if (any(parent>n) || any(parent<0) || any(parent~=floor(parent)) ...
	       || any(parent==(1:n)))
        error('Bad vector of parent pointers.');
     end
     %[parent,pv] = fixparentCANONICAL(parent);
     %[parent,pv] = fixparent(parent);
     [parent,pv] = convert2postorder(parent);
 end

if nargin < 2,

    % Create the adjacency matrix A of the given tree,
    % and get the postorder with another call to etree.

    j = find(parent);
    A = sparse (parent(j), j, 1, n, n);
    A = A + A' + speye(n,n);
    [ignore, post] = etree(A); %#ok<ASGLU>
%    post
end;

% Add a dummy root node #n+1, and identify the leaves.

parent = rem(parent+n, n+1) + 1;  % change all 0s to n+1s
isaleaf = ones(1,n+1);
isaleaf(parent) = 0;%zeros(n,1);

% In postorder, compute heights and descendant leaf intervals.
% Space leaves evenly in x (in postorder).

xmin = n(1,ones(1,n+1)); % n+1 copies of n
xmax = zeros(1,n+1);
nleaves = 0;

for i = 1:n,
    node = post(i);
    if isaleaf(node),
        nleaves = nleaves+1;
        xmin(node) = nleaves;
        xmax(node) = nleaves;
    end;
    dad = parent(node);
    %RG
%    height(dad) = max (height(dad), height(node)+1);
    xmin(dad)   = min (xmin(dad),   xmin(node));
    xmax(dad)   = max (xmax(dad),   xmax(node));
end;

% Compute coordinates, leaving a little space on all sides.

deltax = 1/(nleaves+1);
x = deltax * (xmin+xmax)/2;

% Omit the dummy node.

x = x(1:n);

if ~isempty(pv)
   x(pv) = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is a recursive version of fixparent that, quite unexpectedly, is far
%faster than its non-recursive cousins. The reason, I think, is a clever
%pre-allocation and pre-calculation of nodes' sons.
function [parents2, pv] = convert2postorder(parents)
parents  = int32(parents);
parents2 = repmat(int32(-1), size(parents));
pv       = repmat(int32(-1), size(parents));
roots    = int32(find(parents==0));
idx2     = int32(1);
%make index matrix of sons
sons     = cellfun(@(x)int32(x), cell(size(parents)), 'uniformoutput', false);
for z=int32(1):numel(parents)
  pz = parents(z);
  if pz>0
    sons{pz}(end+1,1) = z;
  end
end
for k=int32(1):numel(roots)
  [parents2,pv,idx2,newroot] = postorder(sons, parents2, pv, idx2,roots(k));
  parents2(newroot) = 0;
end
parents2 = double(parents2);
pv       = double(pv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The recursive algorithm to order the tree in postorder
function [parents2,pv,idx2,root] = postorder(sons, parents2, pv, idx2, root)

childs    = sons{root};
childsnew = zeros(size(childs), 'int32');
for k=int32(1):numel(childs)
  [parents2,pv,idx2,childsnew(k)] = postorder(sons, parents2, pv, idx2, childs(k));
end
parents2(childsnew) = idx2;
pv(idx2) = root;
root     = idx2;
idx2     = idx2 + 1;