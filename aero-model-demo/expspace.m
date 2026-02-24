function space = expspace(min, max, n_points)

gf = (max / min)^(1 / (n_points - 1));
space = min * gf.^((1:n_points) - 1);

end