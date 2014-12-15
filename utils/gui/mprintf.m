%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make a simultaneous fprint to screen and one or several files
function mprintf(file, varargin)
outstring = sprintf(varargin{:});
fprintf(outstring);
if iscell(file)
  for k=1:numel(file)
    fprintf(file{k}, outstring);
  end
else
  fprintf(file, outstring);
end
