function d = distanceMatrix(puntos, puntos2)
%Distancias entre conjuntos de vértices
%
% Devuelve las distancias entre los puntos de un solo conjunto (matriz
% NX2), o entre los puntos de dos conjuntos distintos (matrices NX2 y MX2).
% Hecho a partir de una utilidad del fileexchange en matlabcentral, a la
% cual hemos añadido opción de hallar las distancias entre todos los
% puntos de un conjunto, y para pares de conjuntos, hemos
% ajustado el código para que trabaje con vectores columna

if (nargin == 1)
    pp       = sum(puntos.*puntos, 2);
    resto    = puntos*puntos';
    auxiliar = repmat(pp, [1 size(pp,1)]);
    d        = realsqrt( abs( auxiliar + auxiliar' - 2*resto ));
else
    puntos  = puntos';
    puntos2 = puntos2';
    aa=sum(puntos.*puntos,1); bb=sum(puntos2.*puntos2,1); ab=puntos'*puntos2; 
    d = realsqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
end
