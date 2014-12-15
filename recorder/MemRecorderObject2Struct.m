function recS = MemRecorderObject2Struct(recO, ss)

if nargin<2
  ss = BasicSystem(1);
end

recS = memRecordingNoNewClasses.MemRecorder(ss, 1, 1, 1, 1);

names = fieldnames(recS);

for k=1:numel(names)
  nam = names{k};
  switch nam
    case {'emptyState', 'fieldSets', 'TYPE'}
    otherwise
      recS.(names{k}) = recO.(names{k});
  end
end
