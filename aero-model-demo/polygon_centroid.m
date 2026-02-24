function [x_bar, y_bar] = polygon_centroid (x_points, y_points, area)

    % Make the coordinates circle back to the first point
    x_coords = [x_points, x_points(1)];
    y_coords = [y_points, y_points(1)];

    % Find the X and Y centroid positions
    x_bar = 0;
    y_bar = 0;
    for j = 1:(length(x_coords) - 1)
        x_bar = x_bar + (x_coords(j) + x_coords(j + 1)) * (x_coords(j) * y_coords(j + 1) - y_coords(j) * x_coords(j + 1));
        y_bar = y_bar + (y_coords(j) + y_coords(j + 1)) * (x_coords(j) * y_coords(j + 1) - y_coords(j) * x_coords(j + 1));
    end
    x_bar = x_bar / (6 * area);
    y_bar = y_bar / (6 * area);

end