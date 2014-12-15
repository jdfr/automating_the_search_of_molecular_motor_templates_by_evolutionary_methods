function sowSeed(rndSeed)

if (~isempty(rndSeed))
  if iscell(rndSeed) && (numel(rndSeed)==1)
    rndSeed = rndSeed{1};
  end
  if ischar(rndSeed)
    rndSeed = hex2num(rndSeed);
  end
  rand('twister', rndSeed);
end
