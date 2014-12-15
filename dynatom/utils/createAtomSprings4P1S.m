function [springR, springK] = createAtomSprings4P1S(pos, springToFix, fixedForceDensity, springsSpec)

eqsystem = makeAtomEquilibriumEquations;
localSprings = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

[pos, forceDensities] = solveSystem(eqsystem, (1:4)', pos, springToFix, fixedForceDensity);

%calculate rest length and rigidity parameter K for new springs
diffs = (pos(localSprings(:,1),:)-pos(localSprings(:,2),:));
springLengths = realsqrt(sum(diffs.*diffs, 2));
forces = forceDensities.*springLengths;
[springR, springK] = calculateSpringSpecs(springLengths, forces, springsSpec);
