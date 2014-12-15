function axiswindow = getAxisWindow(pos, margen)
mins = min(pos);
maxs = max(pos);
intervs = maxs-mins;
meds = (mins+maxs)/2;
if nargin<2
  margen = 1.1;
end
if margen<0
  maxi = (max(intervs)-margen)/2;
else
  maxi = max(intervs)*(1+margen)/2;
end
if size(pos, 2)<3
  axiswindow = [meds(1)+[-1 1]*maxi, meds(2)+[-1 1]*maxi];
else
  axiswindow = [meds(1)+[-1 1]*maxi, meds(2)+[-1 1]*maxi, meds(3)+[-1 1]*maxi];
end
end
