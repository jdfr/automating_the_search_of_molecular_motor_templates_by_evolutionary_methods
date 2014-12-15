function d = distanceMatrixSquared(puntos, puntos2)

if (nargin == 1)
    pp       = sum(puntos.*puntos, 2);
    resto    = puntos*puntos';
    auxiliar = repmat(pp, [1 size(pp,1)]);
    d        = abs( auxiliar + auxiliar' - 2*resto );
else
    puntos  = puntos';
    puntos2 = puntos2';
    aa=sum(puntos.*puntos,1); bb=sum(puntos2.*puntos2,1); ab=puntos'*puntos2; 
    d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
end
