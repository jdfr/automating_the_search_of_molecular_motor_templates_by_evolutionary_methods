function genome = DAcreateRandomGenome(mode, varargin)
%create a genome whose genes' parameters are randomly chosen. Parameters:
%   -mode: if 'geneTypesList', the remaining parameters are:
%           -geneTypesList: list of gene types
%           -geneTypes: optional, structure containing geneType codes
%          if 'allowedTypesTable', the remaining parameters are:
%           -numGenes: length of the genome
%           -allowedTypesTable: optional, table of allowed gene types  
%                               and the relative probability of occurrence  
%                               for each one (a default one is used if no
%                               table is provided)
%           -geneTypes: optional, structure containing geneType codes


switch mode(1)
  case 'g'%'geneTypesList'
    geneTypesList = varargin{1};
    geneTypes = extractGeneTypes(varargin, 2);
    doLoop = false;
  case 'a'%'allowedTypesTable'
    numGenes = varargin{1};
    mutationArgs = varargin{2};
    allowedTypesTable = mutationArgs.allowedTypesTable;
    geneTypes = mutationArgs.geneTypes;
    geneDistribution = mutationArgs.geneDistribution;
%       allowedTypesTable = [...
%         geneTypes.switchForcement                           1; ...
%         geneTypes.mitosisSpringShift                        2; ...
%         geneTypes.mitosisStatic                             1.5; ...
%         geneTypes.mirrorSprings                             1; ...
%         geneTypes.permuteShiftMatrix                        1; ...
%         geneTypes.changeGeneDistribution                    1; ...
%         geneTypes.changeAttachedSpringsDistribution         0.5; ...
%         geneTypes.changeSpringR                             0.1; ...
%         ...
%         ...%geneTypes.changeMitosisTime, ...
%         ...%geneTypes.waitSomeTime, ...
%         ...%geneTypes.changeSpringK, ...
%         ...%geneTypes.changeEdgeTime, ...
%         ];
    geneTypesList = makeGeneTypesList(allowedTypesTable, numGenes);
    doLoop = numel(varargin)>2;
    if doLoop
      maxCells = varargin{3};
      %if the max number of cells is declared, first of all, do not
      %generate more mitosis genes than allowed cells.
      geneTypesList = testGenomeTypesList(geneTypesList, geneDistribution, geneTypes, maxCells, allowedTypesTable, numGenes);
    end
  otherwise
    error('unrecognized mode (%s)!!!', mode);
end

goon = true;

while goon
  baregen = struct('type', {[]}, 'args', {{}}); 

  genome = cell(1, numGenes);

  for k=1:numGenes
    gen = baregen;
    gen.type = geneTypesList(k);
    %for each gene type with parameters, make random choices for the
    %parameters
    switch gen.type
        case geneTypes.switchForcement
          spring      = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2'}, ones(1,6), 1);
          both        = pickOption([true false], [1 1], 1);
          generations = pickOption(2:7, ones(1,6),1);
          switchtype  = pickOption([true false], [1 1], 1);
          gen.args = {spring{1} both generations switchtype};
        case geneTypes.permuteShiftMatrix
          gen.args = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2', 'V', 'H'}, ones(1,8), 1);
        case geneTypes.changeGeneDistribution
          gen.args = pickOption({[] -1 0 1}, ones(1,4), 1);
        case geneTypes.changeAttachedSpringsDistribution
          gen.args = {pickOption([-1 0 1], ones(1,3), 1)};
        case geneTypes.changeMitosisTime
          gen.args = {pickOption([-1 -2 1 2], ones(1,4), 1)};
        case geneTypes.waitSomeTime
          gen.args = {pickOption([0 1 2 3], ones(1,4), 1)};
        case geneTypes.changeSpringR
          opt = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2'}, ones(1,6), 1);
          if pickOption(0:1, [1 1], 1)
            ratio = affinTransform(rand, [0.2 0.9], [0 1]); %pick a ratio from 0.2...0.9
          else
            ratio = affinTransform(rand, [1.1 3], [0 1]); %pick a ratio from 1.1...3
          end
          gen.args = {opt{1} ratio};
        case geneTypes.changeSpringK
          opt = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2'}, ones(1,6), 1);
          if pickOption(0:1, [1 1], 1)
            ratio = affinTransform(rand, [0.2 0.9], [0 1]); %pick a ratio from 0.2...0.9
          else
            ratio = affinTransform(rand, [1.1 3], [0 1]); %pick a ratio from 1.1...3
          end
          gen.args = {opt{1} ratio};
        case geneTypes.changeEdgeTime
          gen.args = {affinTransform(rand, [0.5 2], [0 1]);};
    end

    genome{k} = gen;
  end
  
  goon = doLoop;
  if goon
    %if max number of cells is provided, make sure we do not surpass it
    goon = DAnumCellsInGenome(genome, geneDistribution, geneTypes, maxCells)>maxCells;
    if goon
      %if so, let's try again;
      geneTypesList = makeGeneTypesList(allowedTypesTable, numGenes);
      geneTypesList = testGenomeTypesList(geneTypesList, geneDistribution, geneTypes, maxCells, allowedTypesTable, numGenes);
    end
  end
  
end

function geneTypes = extractGeneTypes(VARARGIN, n)
if numel(VARARGIN)<n
  ss = AtomSystem;
  geneTypes = ss.geneTypes;
  clear ss;
else
  geneTypes = VARARGIN{n};
end

function geneTypesList = makeGeneTypesList(allowedTypesTable, numGenes)
geneTypesList = pickOption(allowedTypesTable(:,1), allowedTypesTable(:,2), numGenes);

function geneTypesList = testGenomeTypesList(geneTypesList, geneDistribution, geneTypes, maxCells, allowedTypesTable, numGenes)
while DAnumCellsInGenome(geneTypesList, geneDistribution, geneTypes)>maxCells
  geneTypesList = makeGeneTypesList(allowedTypesTable, numGenes);
end
