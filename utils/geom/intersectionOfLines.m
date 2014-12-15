%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uA,uB] = intersectionOfLines(lineA, lineB)
%lineA=[x1,y1,x2,y2], lineB = [x31,y31,x41,y41;x31,y31,x41,y41;...].
%
% Lines are interpreted as in these vectorized equations:
%  P(lineA) = P1+u*(P2-P1)
%  P(lineB) = P3+u*(P4-P3)
%
%This function calculates intersection parameters between line specified in
%lineA and lines specified in lineB, so intersection points will be
%
% [x1+uA*(x2-x1), y1+uA*(y2-y1)]
% and, equivalently:
% [x3+uB*(x4-x3), y3+uB*(y4-y3)]

%extract and name variables, to avoid a crazy long equation
x1 = lineA(:,1);
y1 = lineA(:,2);
x2 = lineA(:,3);
y2 = lineA(:,4);

x3 = lineB(:,1);
y3 = lineB(:,2);
x4 = lineB(:,3);
y4 = lineB(:,4);

%name differences, also to avoid crazy long equations
x2x1 = x2-x1;
y2y1 = y2-y1;
x4x3 = x4-x3;
y4y3 = y4-y3;

x1x3 = x1-x3;
y1y3 = y1-y3;

%name component common to both equations
denominator = (x2x1.*y4y3 - y2y1.*x4x3);

%calculate line parameters
uA = - (y4y3.*x1x3 - x4x3.*y1y3) ./ denominator;
uB = - (y2y1.*x1x3 - x2x1.*y1y3) ./ denominator;

%                       x1 y4y3 - x3 y4y3 - x4x3 y1 + x4x3 y3
%                uA=  - -------------------------------------
%                               x2x1 y4y3 - y2y1 x4x3
%  
%                      -x2x1 y1 + x2x1 y3 + y2y1 x1 - y2y1 x3
%                uB=  - --------------------------------------
%                              x2x1 y4y3 - y2y1 x4x3

end
