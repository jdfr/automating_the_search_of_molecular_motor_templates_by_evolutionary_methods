function map = owncolormap(n)

c = [0     0     0     ; ...
     0.75  0.75  0     ; ...
     0     0.75  0.75  ; ...
     0     0.5   0     ; ...
     0.75  0     0.75  ; ...
    ];
     

map = c(rem(0:n-1,size(c,1))+1,:);
end
