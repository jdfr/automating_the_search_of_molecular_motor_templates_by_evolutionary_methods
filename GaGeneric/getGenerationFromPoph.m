function P = getGenerationFromPoph(poph, gen)

fns = fieldnames(poph);

gen = poph.generation == gen;

P = makeEmptyStructure(fns(1:poph.numFieldsPop));

for k=1:poph.numFieldsPop
  P.(fns{k}) = poph.(fns{k})(gen);
end

if isfield(poph, 'nsons')
  P.nsons = poph.nsons(gen);
end