function [etot es springLengths] = calculatePotentialEnergy(ss)
  springVectors = ss.pos(ss.springEnds(:,2),:)-ss.pos(ss.springEnds(:,1),:);
  springLengths = realsqrt(sum(realpow(springVectors,2), 2));
  springDisplacements = ss.r-springLengths;
  es   = ss.k.*realpow(springDisplacements,2)/2;
  etot = sum(es);

