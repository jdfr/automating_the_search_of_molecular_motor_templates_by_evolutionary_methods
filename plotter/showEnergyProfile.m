function [h e ee ts] = showEnergyProfile(rec, limitA, limitB)

h = figure;
[e ts] = energyProfile(rec, 0, limitA, limitB);
ee = diff(e)./diff(ts);
%[mini maxi] = localMM(ee, 500, true);
L = extr(ee,true);
maxi = L{1};
mini = L{2};

plot(ts,e, ...
     ts(1:end-1),ee, ...
     ts(maxi), ee(maxi), 'or', ...
     ts(mini), ee(mini), 'oc', ...
     ts(maxi), e(maxi), 'or', ...
     ts(mini), e(mini), 'oc' ...
   );
grid on;
axis square;