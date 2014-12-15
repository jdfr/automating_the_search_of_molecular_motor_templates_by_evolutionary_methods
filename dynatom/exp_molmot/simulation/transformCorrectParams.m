%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evparams = transformCorrectParams(evparams)
%this function is here for backwards compatibility with older parameter
%configurations, which didn't group all corrections in a single
%substructure
wnames      = fieldnames(evparams.walker);

corrections = strmatch('correctMistake', wnames);

if isempty(corrections)
  if ~isfield(evparams.walker, 'correct')
    evparams.walker.correct = struct;
  end
else

  wnames                  = reshape(wnames(corrections), 1, []);

  newnames                = cellfun(@(x)[lower(x(15)) x(16:end)], wnames, 'uniformoutput', false);

  values                  = cellfun(@(x)evparams.walker.(x), wnames, 'uniformoutput', false);
  
  evparams.walker         = rmfield(evparams.walker, wnames);
  
  aux                     = [newnames; values];
  
  evparams.walker.correct = struct(aux{:});
end
