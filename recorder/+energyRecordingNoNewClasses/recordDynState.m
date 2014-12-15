function rec = recordDynState(rec, ss, Ts, states)
  %make sure that room is available
  while (rec.index+numel(Ts))>numel(rec.e)
    rec.e(end+rec.inc,:)  = 0;
    rec.ts(end+rec.inc,:) = 0;
  end
  indexes1 = rec.index;
  indexes2 = rec.index+numel(Ts)-1;
  if (rec.index>1) && (Ts(1)==rec.ts(rec.index-1))
    indexes1 = indexes1+1;
    states = states(:,2:end);
  end
  rec.ts(indexes1:indexes2) = Ts;
  if rec.straightforward
    %fast way
    npoints = size(ss.pos,1);
    svel1   = ss.selector.vel(1,1);
    svel2   = ss.selector.vel(end,2);
    %now, extract posX and posY
    selPX = ss.selector.pos(1,1):ss.selector.pos(1,2);
    selPX = reshape(selPX, numel(selPX)/2, 2);
    selPY = selPX(:,2);
    selPX = selPX(:,1);
    %make them relative to spring ends
    selPX1 = selPX(ss.springEnds(:,1));
    selPX2 = selPX(ss.springEnds(:,2));
    selPY1 = selPY(ss.springEnds(:,1));
    selPY2 = selPY(ss.springEnds(:,2));
    %now, calculate energy in a vectorized way
    vp = realpow(states(svel1:svel2, :), 2);
    e_kinetic = sum(bsxfun(@times, vp(1:npoints,:)+vp((npoints+1):end,:), ss.m));
    springLengths = realsqrt(realpow(states(selPX2,:)-states(selPX1,:), 2) + realpow(states(selPY2,:)-states(selPY1,:), 2));
    springDisplacements = bsxfun(@minus, ss.r, springLengths);
    e_potential = sum(bsxfun(@times, ss.k, realpow(springDisplacements,2)));
    rec.e(indexes1:indexes2) = (e_kinetic+e_potential)/2;
    rec.index = indexes2+1;
  else
    %do it the veeeeeery slow way
    for k=1:numel(Ts)
      ssOtro = dumpEvolvedState(ss, states(:,k));
      rec.e(rec.index) = calculateEnergy(ssOtro);
      rec.index = rec.index+1;
    end
  end
  recordingCallback(rec, ss, Ts, states); %do not forget to call recordingCallback!!!
end
