%return weibull's mode location
function moda = modaWeibul(l,k)
moda = l*(((k-1)/k).^(1/k)); 
