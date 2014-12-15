function y = sigmoid01(x, slope)

excess = 2/(exp(slope/2)-1); %to be sure that sigmoid(0)=0 and sigmoid(1)=1

y = (1+excess)./(1+exp(-slope*(x-0.5))) - excess/2;
