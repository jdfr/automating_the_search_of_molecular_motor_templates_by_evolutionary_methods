function x = antisigmoid01(y, slope)

excess = 2/(exp(slope/2)-1); %to be sure that antisigmoid(0)=0 and antisigmoid(1)=1

x = 0.5 - log( (1+excess)./(y+excess/2) -1 ) / slope;


