function e = calculateEnergy(ss)
  springVectors = ss.pos(ss.springEnds(:,2),:)-ss.pos(ss.springEnds(:,1),:);
  springLengths = realsqrt(sum(realpow(springVectors,2), 2));
  springDisplacements = ss.r-springLengths;
  e_pot = sum( ss.k.*realpow(springDisplacements,2) )/2;
  e_vel = sum(sum(realpow(ss.vel,2),2).*ss.m)/2;
  e = e_pot+e_vel;
end
