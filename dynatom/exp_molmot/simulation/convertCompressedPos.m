function [output snapped err] = convertCompressedPos(input, evparams, failWithError)

% if isstruct(input)
%   qtz   = input.quant;
%   qtz.q = quantizer('ufixed', [qtz.nbits qtz.nfrac]);
%   input = input.pos;
% else
%   qtz = evparams.genome.quant;
% end

if isnumeric(input)
  qtz = evparams.genome.quant;
  if isfloat(input)
    isfixed      = ~isempty(qtz.q);
    tipo         = ['uint' num2str(qtz.nbits)];
    if isfixed
      [qmin qmax]= range(qtz.q);
    else
      qmin       = 0;
      qmax       = double(cast(inf, tipo))/realpow(2, qtz.nfrac);
    end
    mins         = min(input);
    maxs         = max(input);
    desps        = (qmax+qmin-mins-maxs)/2;
    err          = any((mins+desps)<qmin);
    if ((nargin<3) || failWithError) && err
      error('This genome could not be compressed in the range [%d, %d]!!!\n%s', qmin, qmax, any2str(input));
    end
    for k=1:3
      input(:,k) = input(:,k)+desps(k);
    end
    if isfixed
      snapped    = snapToGrid(input, qtz.q);
    else
      snapped    = input;
    end
    output       = cast(snapped*realpow(2, qtz.nfrac), tipo);
  else
    output       = double(input)/realpow(2, qtz.nfrac);
  end
end
