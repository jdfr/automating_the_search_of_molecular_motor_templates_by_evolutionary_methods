function signal = eraseSaltAndPepper(signal, unit)
%erases each point i that is just one unit above i+1 and i-1

if nargin<2
  unit=1;
end

ds = diff(signal);

pepper = find((ds(1:end-1)==unit)  & (ds(2:end)==-unit)) + 1;
signal(pepper) = signal(pepper) - unit;
ds(pepper) = 0;

salt   = find((ds(1:end-1)==-unit) & (ds(2:end)==unit))  + 1;
signal(salt)   = signal(salt)   + unit;
