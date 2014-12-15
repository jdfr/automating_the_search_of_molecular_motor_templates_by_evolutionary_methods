%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = getFirstUnchangedAncestor(poph, ind, isparent)
if nargin<3
  isparent = false;
end
parent = poph.tree.parents(ind);
if strcmp(poph.tree.change{ind}, '=') && (parent>0) 
  if ( isfield(poph, 'rndSeed')     && strcmp( poph.rndSeed{ind},    poph.rndSeed{parent}    ) ) || ...
     ( isfield(poph, 'rndSeedDev')  && strcmp( poph.rndSeedDev{ind}, poph.rndSeedDev{parent} ) )
    idx = getFirstUnchangedAncestor(poph, parent, true);
  elseif (poph.individual(ind)==1) && poph.params.elitism
    idx = ind;
  end
elseif strcmp(poph.tree.change{ind}, '[REVERTED]') && (parent>0)
  grandfather = poph.tree.parents(parent);
  idx = getFirstUnchangedAncestor(poph, grandfather, true);
elseif isparent && isnan(poph.fitness(ind)) && (parent>0) && isfield(poph.params, 'rewindNans') && poph.params.rewindNans
  idx = getFirstUnchangedAncestor(poph, parent, true);
else
  idx = ind;
end
