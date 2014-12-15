function poph = addNSons(poph)
validTree = (poph.tree.gD>=min(poph.generation)) & (poph.tree.gD<=max(poph.generation));
nsons = cellfun('prodofsize', poph.tree.sons(validTree));
[nevermind order] = sort(poph.tree.idxD(validTree)); clear nevermind; %#ok<ASGLU>
poph.nsons = nsons(order)';
