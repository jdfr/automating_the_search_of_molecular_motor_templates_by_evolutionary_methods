function nc = DAnumCellsInGenome(genome, geneDistribution, geneTypes, maxCells)

if isempty(genome)
  nc = 1;
elseif isnumeric(genome)
  %quick and dirty approximation
  ms = [geneTypes.mitosisSpringShift, geneTypes.mitosisStatic];
  nc = 1 + sum((genome==ms(1)) | (genome==ms(2)));
elseif isstruct(genome{1})
  %genesMitosis = [geneTypes.mitosisSpringShift, geneTypes.mitosisStatic];
  if nargin<4
    maxCells = inf;
  end
  %nc2 = numCellsAux(1, genome, geneTypes, genesMitosis, geneDistribution, maxCells);
  nc = numCellsAuxIter(genome, geneTypes, geneDistribution, maxCells);
%   if nc~=nc2
%     error('heyyyy!');
%   end
  
%   ms = [geneTypes.mitosisSpringShift, geneTypes.mitosisStatic];
%   nc = 1 + sum(cellfun(@(x)any(x.type==ms), genome));
else
  error('parameter not understood!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc = numCellsAuxIter(genome, geneTypes, distr, maxCells)
types = zeros(size(genome));
for k=1:numel(genome)
  types(k) = genome{k}.type;
end
%nc = numCellsAuxIterMX(types, geneTypes.mitosisSpringShift, geneTypes.mitosisStatic
isMitosis  = (types==geneTypes.mitosisSpringShift) | (types==geneTypes.mitosisStatic);
isChgDistr = types==geneTypes.changeGeneDistribution;
nc         = 1;

stack      = [1, numel(genome) distr; zeros(sum(isMitosis),3)];
idxS       = 1;

while idxS>0
  k   = stack(idxS, 1);
  ind = stack(idxS, 2);
  if k>ind
    idxS = idxS - 1;
    continue;
  end
  if isMitosis(k)
    nc = nc+1;
    if nc>maxCells
      return
    end
    [iA iB] = DAsplitGenome(genome(k+1:ind), stack(idxS, 3), true);
    iA      = iA+k;
    iB      = iB+k;
    if isempty(iB)
      stack(idxS,[1 2]) = iA;
    elseif isempty(iA)
      stack(idxS,[1 2]) = iB;
    else
      if idxS==size(stack, 1)
        stack  = [stack,  zeros(size(stack))]; %#ok<AGROW>
      end
      stack(idxS,[1 2]) = iA;
      idxS              = idxS + 1;
      stack(idxS,:)     = [iB stack(idxS-1,3)];
    end
  elseif isChgDistr(k)
    distr = genome{k}.args{1};
    if isempty(distr)
      distr = nan;
    end
    stack(idxS, [1 3]) = [k+1 distr];
  else
    stack(idxS, 1) = k+1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc = numCellsAux(nc, genome, geneTypes, genesMitosis, distr, maxCells)
if nc>maxCells;
  return;
end
while numel(genome)>0
  gen = genome{1};
  genome = genome(2:end);
  if gen.type==geneTypes.changeGeneDistribution
    distr = gen.args{1};
  elseif ismember(gen.type, genesMitosis)
    [gA gB] = DAsplitGenome(genome, distr);
    nc = numCellsAux(nc+1, gA, geneTypes, genesMitosis, distr, maxCells);
    if nc>maxCells;
      return;
    end
    nc = numCellsAux(nc, gB, geneTypes, genesMitosis, distr, maxCells);
    if nc>maxCells;
      return;
    end
    genome = {};
  end
end



