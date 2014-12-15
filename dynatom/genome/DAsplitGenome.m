%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gA gB] = DAsplitGenome(gA, geneDistribution, mode)
%split a genome

%ratioIndicator is either [], -1, 0, or 1
if isempty(geneDistribution) || isnan(geneDistribution)
  if (nargin<3) || (~mode)
    gB = gA;
  else
    gA = [1 numel(gA)];
    gB = gA;
  end
elseif isnan(geneDistribution)
  error('geneDistribution cannot be NAN!!! Surely, it has not been initialized...');
else
  ratios = [0 0.5 1];
  ratio = ratios(sign(geneDistribution)+2);
  % switch sign(geneDistribution)
  %   case +1
  %     ratio = 0;
  %   case 0
  %     ratio = 0.5;
  %   case -1
  %     ratio = 1;
  % end

  numGenes = numel(gA);
  numGenesRetenidos = round(ratio*numGenes);
  if (nargin<3) || (~mode)
    gB = gA((numGenesRetenidos+1):end);
    gA = gA(1:numGenesRetenidos);
  else
    gB = [(numGenesRetenidos+1) numel(gA)];
    gA = [1 numGenesRetenidos];
  end
end
end
