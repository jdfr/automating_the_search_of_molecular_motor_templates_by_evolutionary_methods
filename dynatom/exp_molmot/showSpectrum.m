function showSpectrum(ss)

[autovecs autovecs autovals] = ANMdecomposition(ss);

autovals       = diag(autovals);

bigEnough      = (autovals>1e-12) ;%| ((autovals>1e-4) & (collectivityB>0.1));
bigEnough(1:3) = false;

autovals       = autovals(bigEnough);

logautovals    = log10(autovals/autovals(1));

ys             = ones(size(autovals));

figure;

subplot(2,1,1);
plot(logautovals,ys,'+');
title('non trivial eigenvalues in normalized log space');
grid on;

subplot(2,1,2);
plot(autovals,ys,'+');
title('actual non trivial eigenvalues');
grid on;

