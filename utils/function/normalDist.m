function densities = normalDist(points, mean, std)
densities = 1/std/sqrt(2*pi)*exp(-realpow((points-mean)/std, 2)/2);