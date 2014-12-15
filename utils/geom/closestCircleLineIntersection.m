function distancesIntersec = closestCircleLineIntersection(circleCenter, circleRadius, points, vector)
%given a circle, a matrix NX2 of points, and a vector, returns the closest point of
%intersection between the sphere and the lines containing the points, all
%of them in the direction given by vector

%equation of circle: (x-cx)^2+(y-cy)^2=r^2
%equations of line: x=x1+u*(x2-x1), y=y1+u*(y2-y1)

%fill points1 and points2
points1         = points;
points2         = bsxfun(@plus, points, vector);
segmentsL       = points2-points1;
segmentsCL      = bsxfun(@minus, points1, circleCenter);
%substituting the line equations in the circle equation, we get a equation
%of type a*u^2+b*u+c=0. The parameters are
a               = sum((segmentsL).^2,2);
b               = 2*sum(segmentsL.*segmentsCL,2);
c               = sum(points1.^2,2) + sum(circleCenter.^2,2) - 2*sum(bsxfun(@times, points1, circleCenter),2) - circleRadius.^2;
roots           = sqrt(b.^2-4*a.*c);
u1s             = (-b+roots)/2./a;
u2s             = (-b-roots)/2./a;
reallyIntersec  = imag(u1s)==0;
%find closest points (the closest point among u1 and u2 will have the minimal U parameter)  
closestUs       = zeros(size(u1s));
u1sIsClosest    = abs(u1s)<abs(u2s);
closestUs(u1sIsClosest)   = u1s(u1sIsClosest);
closestUs(~u1sIsClosest)  = u2s(~u1sIsClosest);
%calculate intersection vectors (a vector that, if his origin is placed on
%point1, will reach the intersection at its tip). To do so, we use the
%parametric line equation
segmentsIntersec  = bsxfun(@times, closestUs, segmentsL);
distancesIntersec = repmat(nan, size(points,1),1);
distancesIntersec(reallyIntersec) = realsqrt(sum(segmentsIntersec(reallyIntersec,:).^2, 2));


%http://local.wasp.uwa.edu.au/~pbourke/geometry/sphereline/