%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vprintf(verbose, varargin)
if verbose
  fprintf(varargin{:});
end
end
