%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%equality test for genomes, data types (logical/numerical) aside
function areEqual = DAgenomesAreEqual(genomeA, genomeB)

if (~iscell(genomeA)) || (~iscell(genomeB))
  error('genomes are cell arrays!');
end

areEqual = numel(genomeA)==numel(genomeB);

if areEqual
  areEqual = all(cellfun(@genesAreEqual, reshape(genomeA, 1, []), reshape(genomeB, 1, [])));
end
end
%%%%%%%%%%%%%
function areEqual = genesAreEqual(genA, genB)

if (~isstruct(genA)) || (~isstruct(genB))
  error('genes are structs!');
end

areEqual = genA.type == genB.type;

if areEqual
  hasargsA = isfield(genA, 'args') && (~isempty(genA.args));
  hasargsB = isfield(genB, 'args') && (~isempty(genB.args));
  if (hasargsA && (~hasargsB)) || (hasargsB && (~hasargsA))
    error('genes of same type must have same number of args! (1)');
  end
  
  if hasargsA && hasargsB
    if numel(genA.args) ~= numel(genB.args)
      error('genes of same type must have same number of args! (2)');
    end

    areEqual = all(cellfun(@isequal, reshape(genA.args, 1, []), reshape(genB.args, 1, [])));
  end
end
end
