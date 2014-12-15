function [summary indexes] = reapIndividualGenealogy(poph, generation, rangeid, individual, complete)
%reap the genealogy of a individual specified by its generation, range and
%index in the population. basedir is the directory of the simulation

indexes    = zeros(0,2);
continuar  = true;
genome     = '';
identifier = poph.triplet2idx(generation, rangeid, individual);
while continuar
  idxTree = find(poph.tree.idxD==identifier, 1);
  idxPop  = find(      poph.idx==identifier, 1);
  continuar = ~isempty(idxPop);
  if continuar
    identifier    = poph.tree.idxA(idxTree);
    if complete || (~strcmp(genome, poph.genome{idxPop}))
      indexes     = [indexes; idxPop, idxTree]; %#ok<AGROW>
      genome      = poph.genome{idxPop};
    end
  end
end

summary                    = struct;
fn                         = fieldnames(poph);
for k=1:poph.numFieldsPop
  summary.(fn{k})          = poph.(fn{k})(indexes(:,1));
end
summary.changes            = cell(size(summary.genome));
summary.changes{end}       = '';
summary.changes(1:(end-1)) = poph.tree.change(indexes(1:(end-1),2));
