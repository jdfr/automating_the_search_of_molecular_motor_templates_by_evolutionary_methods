%wl and wk control the profile of the intensity depending on the distance.
%K is the coupling constant for this force. tol is a minimal limit for
%defining a cutoff distance. It is guarranteed that, at cutoff distance,
%the normalized weibull function will be about the same order of magnitud
%as tol (it may differ, but not too much)
function chargeParams = ChargeWParams(wl, wk, K, tol)

chargeParams = struct('wl', wl, 'wk', wk, 'K', K, 'tol', tol, 'moda', modaWeibull(wl,wk), 'cutoff', K*weibullMaxDist(wl,wk,K,tol));

cs = 0:0.01:1; d=0:0.01:10; hold on; for k=1:numel(cs); y=1./(((d).^2)+c(k)); plot(d,y); end; grid on; axis square; hold off;