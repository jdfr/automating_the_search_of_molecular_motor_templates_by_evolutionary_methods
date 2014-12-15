function d = pointPolyDist(point, poly)
%point-to-polygon distance. The distance is negative if the point is inside
%the polygon. I picked up a code from the file exchange, but it was faulty:
%the programmer was careless enough to use a line representation tha
%didn't mixed well with segments parallel to X/Y axis between points and
%polygons. So I have had to rewrite it from scratch. Still, I owe him an
%idea: to use 'inpolygon'

%if polygon is not closed, close it
if any(poly(1,:)~=poly(end,:))
  poly = [poly; poly(1,:)];
end

%calculate intersections between polygonal lines and lines perpendicular to
%them passing by the point

%vectors parallel to polygonal lines
vecpar = diff(poly);
% %vectors perpendicular to polygonal lines. This is not strictly necessary
% vecper = [-vecpar(:,2), vecpar(:,1)];

%polygonal lines are of the form ANYPOINT = POINT1 + U * (POINT2-POINT1)
%polygonal lines' perpendiculars are of the form ANYPOINT = POINT0 + V *PERPENDICULAR(POINT2-POINT1) 
%where point2, point1 define a polygonal line, and point0 is the point
%whose distance to the polygon we are seeking. So, we solve a system of
%equations for each line, equally ANYPOINT in both equations and for both
%variables X and Y
U = -((vecpar(:,1).*(poly(1:end-1,1) - point(1))) + (vecpar(:,2).*(poly(1:end-1,2) - point(2))) ) ./ sum(realpow(vecpar, 2), 2);

% find all cases where projected point is inside the segment: the parameter
% will be between 0 and 1 (ranging from one point to the other in the
% segment)
idx = (U>=0) & (U<=1);

%find projection point for those that intersect inside the segment
proj = poly(idx,:);
if ~isempty(proj)
  proj(:,1) = proj(:,1) + U(idx).*vecpar(idx,1);
  proj(:,2) = proj(:,2) + U(idx).*vecpar(idx,2);
end
%add the polygonal vertices to the projection points, since the closest
%point might be also a vertex
proj = [proj; poly(1:end-1,:)];


%find distance from polygon points (projections and vertices) to point 0
dist = realsqrt( realpow(proj(:,1)-point(1), 2) + realpow(proj(:,2)-point(2), 2) );

%find minimal distance
d = min(dist);
%[mind id] = min(dist);
%idx = find(idx);

%invert distance if the point is inside the polygon
if(inpolygon(point(1), point(2), poly(:,1), poly(:,2))) 
   d = -d;
end
