function k=curvature(ts, e)

norm = @(x) (x-min(x))/(max(x)-min(x));

%e is in the range 0,1
e = norm(e);
%ts is in the range 0,1
ts = (ts-ts(1))/(ts(end)-ts(1));

step = ts(2)-ts(1);

%first derivative
de        = diff(e)/step;
ts_de     = ts(1:(end-1))+step/2;
%second derivative
dde       = diff(de);%/step;

% %smoothing
% n=1;
% ddef = dde;
% ddef = [repmat(ddef(1), n, 1); ddef; repmat(ddef(end), n, 1)];
% ddef = medfilt2(ddef, [(1+n*2), 1]);
% ddef = ddef((n+1):(end-n));

ts_dde    = ts(2:(end-1));
%first derivative, interpolated to be defined at the same timesteps as the
%second derivative
center_de = interp1(ts_de,de,ts_dde,'spline');
%curvature
k = dde./realsqrt(realpow(1+realpow(center_de,2), 3));

% %normalize
% center_de=norm(-center_de);
% dde=dde/(max(dde)-min(dde));
% k=k/(max(k)-min(k));
% 
% % plot(ts, e, 'b', ts_dde, center_de, 'r', ts_dde, dde, 'g', ts_dde, k, 'k'); grid on;
% % legend({'función', 'derivada', 'segunda derivada', 'curvatura'});
% plot(ts, e, 'b', ts_dde, center_de, 'r', ts_dde, k, 'k'); grid on;
% legend({'función', 'derivada', 'curvatura'});
% % plot(ts, e, 'b', ts_dde, k, 'k'); grid on;
% % legend({'función', 'curvatura'});
