function rec = recordAllState(rec, ss)
  if (rec.index>1) && (ss.t==rec.ts(rec.index-1))
    rec.e = calculateEnergy(ss);
  end
  %make sure that room is available
  if rec.index>numel(rec.e)
    rec.e(end+rec.inc,:)  = 0;
    rec.ts(end+rec.inc,:) = 0;
  end
  %extract state and record it
  rec.ts(rec.index) = ss.t;
  rec.e(rec.index)  = calculateEnergy(ss);
  rec.index = rec.index+1;
  %to harvest energy information is straightforward if all of these
  %hold:
  %   -nº of dimensions is 2
  %   -all points are dynamical_p
  %   -no point is dynamical_m, no spring is neither dynamical_r nor dynamical_k  
  rec.straightforward = (size(ss.pos,2)==2) && all(ss.dynamical_p) && (~(any(ss.dynamical_r) || any(ss.dynamical_k) || any(ss.dynamical_m) ));
  recordingCallback(rec, ss); %do not forget to call recordingCallback!!!
end
