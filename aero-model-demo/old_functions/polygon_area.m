function area = polygon_area (x_points, y_points)

    % Make the coordinates circle back to the first point
    x_coords = [x_points, x_points(1)];
    y_coords = [y_points, y_points(1)];
    
    % Find the Surface areas using Shoelace's method
    area = 0;
    for j = 1:(length(x_coords) - 1)
        area = area + (x_coords(j) * y_coords(j + 1) - y_coords(j) * x_coords(j + 1));
    end
    area = area / 2;

end