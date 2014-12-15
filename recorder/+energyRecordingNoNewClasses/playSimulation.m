%specially tailored "playSimulation" since we only have energy
%information, only energy is passed down to the function
function playSimulation(rec, Ts, fun)
  es = rec.e(1:(rec.index-1));
  TS = rec.ts(1:(rec.index-1));
  [nevermind, loc] = ismember(Ts, TS);
  clear nevermind;
  loc = loc(loc~=0);
  arrayfun(fun, TS(loc), es(loc));
end
