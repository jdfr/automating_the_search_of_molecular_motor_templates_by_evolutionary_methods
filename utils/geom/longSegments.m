%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function longs = longSegments(segments)
  longs = realsqrt(sum(realpow(segments, 2), 2));
end
