%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = showError(error)
if isempty(error)
  str='<EMPTY ERROR>';
elseif (isa(error, 'MException') || isstruct(error)) && all(ismember({'message', 'identifier', 'stack'}, lower(fieldnames(error))))
  lines = cell(size(error.stack));
  for m=1:numel(error.stack)
      lines{m} = sprintf('    File: %s\n  Name: %s\n  Line: %g\n', error.stack(m).file, error.stack(m).name, error.stack(m).line);
  end
  str = horzcat(...
      sprintf('ERROR:\n  Message:    <<%s>>\n  Identifier: <<%s>>\n  Stack:\n', error.message, error.identifier), ...
      lines{:});
else
  str=['<NO STRUCT ERROR: ' any2str(error) '>'];
end
str = strrep(str, '\', '\\');
