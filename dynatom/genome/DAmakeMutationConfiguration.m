function [genomeParams mutateFun makeupParamsFun createRndFun] = DAmakeMutationConfiguration(geneTypes)

genomeParams.initNumGenes      = 30;
genomeParams.geneTypes         = geneTypes;
mutationArgs.poissons = {0.5, 0:10, poisspdf(0:10, 0.5); ...
                         1,   0:12, poisspdf(0:12, 1);};  
genomeParams.lambdaMutationIdx = 1;
genomeParams.maxCells          = 60;30;
genomeParams.maxGenes          = 80;
genomeParams.geneDistribution  = nan;
mutationArgs.geneTypes         = geneTypes;
mutationArgs.mutationTable     = {...
  'dpG',   020, @mutationDupGen;          ...
  'dpS',   000, @mutationDupSegment;      ...
  'ins',   040, @mutationInsGen;          ...
  'rpl',   020, @mutationReplaceGen;      ...
  'del',   020, @mutationDelSegment;      ...
  'shf',   020, @mutationShiftSegment;    ...
  'chg',   100, @mutationChgGeneParam     ...
  };
mutationArgs.allowedTypesTable = [...
  geneTypes.mitosisSpringShift                        2; ...
  ...%geneTypes.mitosisStatic                             1.5; ...
  geneTypes.mirrorSprings                             1; ...
  geneTypes.permuteShiftMatrix                        1; ...
  geneTypes.changeGeneDistribution                    1; ...
  geneTypes.changeAttachedSpringsDistribution         1;%0.5; ...
  geneTypes.changeSpringR                             1;%0.1; ...
  ...
  ...%geneTypes.switchForcement                           1; ...
  ...%geneTypes.changeMitosisTime, ...
  ...%geneTypes.waitSomeTime, ...
  ...%geneTypes.changeSpringK, ...
  ...%geneTypes.changeEdgeTime, ...
  ];

genomeParams.mutationArgs      = mutationArgs;
mutateFun                      = @DAmutageGenome;
makeupParamsFun                = @DAmakeupParams;
createRndFun                   = @DAcreateRandomGenomes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = DAmakeupParams(params)
params.genome.mutationArgs.geneDistribution = params.evparams.ss.atomDefaultVals{4}{1}.geneDistribution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function genomes = DAcreateRandomGenomes(tam_poblacion, params)
genomes = cell(tam_poblacion, 1);
initNumGenes = params.genome.initNumGenes;
mutationArgs = params.genome.mutationArgs;
maxCells = params.genome.maxCells;

for k=1:numel(genomes)
  % 'a' stands for 'allowedTypesTable'
  genomes{k} = DAcreateRandomGenome('a', initNumGenes, mutationArgs, maxCells);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome cad maskNewThis] = DAmutageGenome(idx, G, params, P, idxInP, maskNew) 
%[G.genome{j} cad maskNew(j)] = params.mutateGenome(j, G, params, P, idxInP, maskNew);

genome = G.genome{idx};

[genome cad maskNewThis] = mutateGenome(genome, params);


function [genome cad masknew] = mutateGenome(genome, params) %#ok<*INUSL>

geneTypes         = params.genome.geneTypes;
lambdaMutationIdx = params.genome.lambdaMutationIdx;
mutationArgs      = params.genome.mutationArgs;
maxCells          = params.genome.maxCells;
maxGenes          = params.genome.maxGenes;
geneDistribution  = mutationArgs.geneDistribution;
genome2str        = params.genome2str;

%piece the mutation table
mutationNames     =              mutationArgs.mutationTable(:,1);
mutationWeights   = cell2arrayMX(mutationArgs.mutationTable(:,2));
mutationFunctions =              mutationArgs.mutationTable(:,3);

goon = true;

oldgenome = genome;

  %numMutations
  %n = 3;
n = pickOption(mutationArgs.poissons{lambdaMutationIdx, [2 3]}, 1);
%n =  poissrnd(lambdaMutation, 1, 1);
if n==0
  cad = [];
  masknew = false;
  return;
end
nts = 0;

while goon
  genome = oldgenome;
  
  %make the mutation choices
  mutations = pickOption(1:numel(mutationNames), mutationWeights, n);

  indexes = cell(n,1);
  segments = cell(n,1);
  changeds = false(n,1);
  
  %perfom the mutations
  for k=1:n
    [genome indexes{k} segments{k} changeds(k)] = mutationFunctions{mutations(k)}(genome, [], mutationArgs); %mutationNames{k}
  end
  nts = nts+1;
  
  goon = (numel(genome)>maxGenes) || (DAnumCellsInGenome(genome, geneDistribution, geneTypes, maxCells)>maxCells);
  if goon
    if nts>200
      error('Too many trials!!!');
    end
    n = pickOption(mutationArgs.poissons{lambdaMutationIdx, [2 3]}, 1);
    %n = poissrnd(lambdaMutation, 1, 1);
    if n==0
      n=1; %at least a mutation
    end
  end
end

masknew = any(changeds);
if masknew
  for k=1:numel(indexes)
    indexes{k}  = mat2str(indexes{k});
    segments{k} = genome2str(segments{k});
  end
  cad = [mutationNames(mutations), indexes, segments];
  % cad = [mutationNames(mutations), ...
  %        cellfun(@mat2str, indexes, 'uniformoutput', false), ...
  %        cellfun(genome2str, segments, 'uniformoutput', false)]';
  cad = horzcat(cad{:});
else
  cad = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome geneIndex segment chgd] = mutationDupGen(genome, mutationName, mutationArgs) %#ok<*INUSD>
ng = numel(genome);
if ng>0
  geneIndex = pickOption(1:ng, ones(1, ng),      1);
  numDups   = pickOption(1:7,  mutationArgs.poissons{1, 3}(1:7), 1);
  segment   = repmat(genome(geneIndex), 1, 1+numDups);
  if size(genome,1)<=size(genome,2)
    genome = [genome(1:geneIndex-1), segment, genome(geneIndex+1:end)];
  else
    genome = [genome(1:geneIndex-1); segment'; genome(geneIndex+1:end)];
  end
  chgd = true;
else
  geneIndex = 0;
  segment = {};
  chgd = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome geneIndex segment chgd] = mutationDupSegment(genome, mutationName, mutationArgs)
ng = numel(genome);
if ng>0
  indexStart = pickOption(1:ng, ones(1, ng),       1);
  numIns     = pickOption(2:13, mutationArgs.poissons{2, 3}(1:12), 1);
  %genome(indexStart:min(indexStart+numIns-1, numel(genome))) = [];
  indexBounds = [indexStart, min(indexStart+numIns-1, ng)];
  segment   = genome(indexBounds(1):indexBounds(2));
  %indexBounds = sort(pickOption(1:ng, ones(1, ng), 2));
  geneIndex = [mat2str(indexStart) ',' mat2str(numIns)];
  if size(genome,1)<=size(genome,2)
    genome = [genome(1:indexBounds(2)), segment, genome((indexBounds(2)+1):end)];
  else
    genome = [genome(1:indexBounds(2)); segment'; genome((indexBounds(2)+1):end)];
  end
  chgd = true;
else
  geneIndex = 0;
  segment = {};
  chgd = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome position gen chgd] = mutationInsGen(genome, mutationName, mutationArgs, position)
ng = numel(genome);
if ng==0
  position = 1;
else
  if nargin<4
      position = pickOption(1:(ng+1), ones(1, ng+1), 1);
  end
end

gen = DAcreateRandomGenome('a', 1, mutationArgs);

if size(genome,1)<=size(genome,2)
  genome = [genome(1:(position-1)), gen{1}, genome(position:end)];
else
  genome = [genome(1:(position-1)); gen{1}; genome(position:end)];
end

chgd = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome geneIndex segment chgd] = mutationReplaceGen(genome, mutationName, mutationArgs)
ng = numel(genome);
position = pickOption(1:ng, ones(1, ng), 1);
geneType = genome{position}.type;
genome(position) = [];
%do not let the new gene to be the same type as the replaced one
mutationArgs.allowedTypesTable(mutationArgs.allowedTypesTable(:,1)==geneType, :) = [];
[genome geneIndex segment chgd] = mutationInsGen(genome, mutationName, mutationArgs, position);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome indexStart segment chgd] = mutationDelSegment(genome, mutationName, mutationArgs)
ng = numel(genome);
if ng>1
  indexStart = pickOption(1:ng, ones(1, ng),      1);
  numDels    = pickOption(1:7,  mutationArgs.poissons{2, 3}(1:7), 1);
  segment    = genome(indexStart:min(indexStart+numDels-1, ng));
  if numel(segment)==numel(genome)
    indexStart = 0;
    segment = {};
    chgd = false;
  else    
    genome(indexStart:min(indexStart+numDels-1, ng)) = [];
    chgd = true;
  end
else
  indexStart = 0;
  segment = {};
  chgd = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome indexPos segment chgd] = mutationShiftSegment(genome, mutationName, mutationArgs)
ng = numel(genome);
if ng>1
  indexStart  = pickOption(1:ng, ones(1, ng),      1);
  numShifts   = pickOption(1:7,  mutationArgs.poissons{2, 3}(1:7), 1);
  indexBounds = [indexStart, min(indexStart+numShifts-1, ng)];
  segment     = genome(indexBounds(1):indexBounds(2));
  if numel(segment)==ng
    indexPos = 0;
    segment = {};
    chgd = false;
    return;
  end
  restGenome  = genome;
  restGenome(indexBounds(1):indexBounds(2)) = [];
  nrg         = numel(restGenome);
  range       = indexStart+[-5 5];%*numShifts;
  if range(2)>nrg
    range = [range(1)-(range(2)-nrg), nrg];
  end
  if range(1)<1
    range = [1, range(2)+(1-range(1))];
  end
  if range(2)>nrg %if so, it means that the range is broader than the rest of the genome itself
    range = [1 nrg];
  end
  indexPos    = round(affinTransform(rand, range, [0 1]));
  if size(genome,1)<=size(genome,2)
    genome = [restGenome(1:(indexPos-1)), segment, restGenome(indexPos:end)];
  else
    genome = [restGenome(1:(indexPos-1)); segment; restGenome(indexPos:end)];
  end
  chgd = indexBounds(1)~=indexPos;
else
  indexPos = 0;
  segment = {};
  chgd = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genome geneIndex segment chgd] = mutationChgGeneParam(genome, mutationName, mutationArgs)
ng = numel(genome);
geneIndex = 0;
segment = {};
chgd = false;
if (ng>0) 
  geneTypes = mutationArgs.geneTypes;
  elegibles = find(genesWithParams(genome, geneTypes));
  ne = numel(elegibles);
  if ne>0
    geneIndex = elegibles(pickOption(1:ne, ones(1, ne), 1));
    oldargs = genome{geneIndex}.args;
    switch genome{geneIndex}.type
      case geneTypes.switchForcement
        argToMutate = pickOption(1:4, ones(1, 4), 1);
        switch argToMutate
          case 1
            %spring index, relative to shiftMatrix
            opt = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2'}, ones(1,6), 1);
            genome{geneIndex}.args{1} = opt{1};
          case 2
            %force indexed spring and its dual
            genome{geneIndex}.args{2} = ~genome{geneIndex}.args{2};
          case 3
            %generations while the effect will last
            num = genome{geneIndex}.args{3};
            while (num==genome{geneIndex}.args{3})
              num = num + pickOption([-1 1 2 3], ones(1,4), 1);
              num = max(min(num, 7), 1); %TODO: completely arbitrary!!!!!
            end
            genome{geneIndex}.args{3} = num;
          case 4
            %also switch type
            genome{geneIndex}.args{4} = ~genome{geneIndex}.args{4};
        end
      case geneTypes.permuteShiftMatrix
        opt = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2', 'V', 'H'}, ones(1,8), 1);
        genome{geneIndex}.args{1} = opt{1};
      case geneTypes.changeGeneDistribution
        opt = pickOption({[] -1 0 1}, ones(1,4), 1);
        genome{geneIndex}.args{1} = opt{1};
      case geneTypes.changeAttachedSpringsDistribution
        genome{geneIndex}.args{1} = pickOption([-1 0 1], ones(1,3), 1);
      case geneTypes.changeMitosisTime
        genome{geneIndex}.args{1} = genome{geneIndex}.args{1}+pickOption([-1 -2 1 2], ones(1,4), 1);
      case geneTypes.waitSomeTime
        genome{geneIndex}.args{1} = genome{geneIndex}.args{1}+pickOption([-1 -2 1 2], ones(1,4), 1);
      case geneTypes.changeSpringR
        opt = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2'}, ones(1,6), 1);
        genome{geneIndex}.args{1} = opt{1};
        if pickOption(0:1, [1 1], 1)
          ratio = affinTransform(rand, [0.2 0.9], [0 1]); %pick a ratio from 0.2...0.9
        else
          ratio = affinTransform(rand, [1.1 3], [0 1]); %pick a ratio from 1.1...3
        end
        genome{geneIndex}.args{2} = ratio;
      case geneTypes.changeSpringK
        opt = pickOption({'V1' 'V2' 'H1', 'H2', 'D1', 'D2'}, ones(1,6), 1);
        genome{geneIndex}.args{1} = opt{1};
        if pickOption(0:1, [1 1], 1)
          ratio = affinTransform(rand, [0.2 0.9], [0 1]); %pick a ratio from 0.2...0.9
        else
          ratio = affinTransform(rand, [1.1 3], [0 1]); %pick a ratio from 1.1...3
        end
        genome{geneIndex}.args{2} = ratio;
      case geneTypes.changeEdgeTime
        ratio = affinTransform(rand, [0.5 2], [0 1]);
        genome{geneIndex}.args{1} = ratio;
    end
    chgd = ~isequal(oldargs, genome{geneIndex}.args);
    segment = genome(geneIndex);
%     if chgd
%       segment = genome(geneIndex);
%     else
%       segment = {};
%     end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conds = genesWithParams(genome, geneTypes)
conds = false(size(genome));
enesWithParams = [...
  geneTypes.switchForcement, ...
  geneTypes.permuteShiftMatrix, ...
  geneTypes.changeGeneDistribution, ...
  geneTypes.changeAttachedSpringsDistribution, ...
  geneTypes.changeMitosisTime, ...
  geneTypes.waitSomeTime, ...
  geneTypes.changeSpringR, ...
  geneTypes.changeSpringK, ...
  geneTypes.changeEdgeTime ...
  ];
for k=1:numel(genome)
  conds(k) = any(genome{k}.type==enesWithParams);
end
  


