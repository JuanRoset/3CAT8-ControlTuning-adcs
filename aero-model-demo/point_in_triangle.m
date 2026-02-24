function inside = point_in_triangle(triangle, point)

% Assume point is outside
inside = false;

% Only do area check if point coordinates li within bounding box of triangle
if point(1) < max(triangle(:, 1)) && point(1) > min(triangle(:, 1)) && point(2) < max(triangle(:, 2)) && point(2) > min(triangle(:, 2))

    % Set default tolerance value
    tolerance = 1e-6;
    
    % Extraxt triangle vertex coordinates
    vertex_A = triangle(1, :);
    vertex_B = triangle(2, :);
    vertex_C = triangle(3, :);
    
    % Compute the three different areas composed of the sides and the point
    Area_A = triangle_area([vertex_A; vertex_B; point]);
    Area_B = triangle_area([vertex_B; vertex_C; point]);
    Area_C = triangle_area([vertex_C; vertex_A; point]);
    
    % Compute the triangle's area and the point's total area
    Area_triangle = triangle_area(triangle);
    Area_point = Area_A + Area_B + Area_C;
    
    % If point area is larger than triangle area, point is outside
    inside = ~(Area_point - Area_triangle >= tolerance);
end

end

function area = triangle_area(triangle)

% Extraxt triangle vertex coordinates
point_1 = triangle(1, :);
point_2 = triangle(2, :);
point_3 = triangle(3, :);

% Apply the determinant area formula
double_area_A = point_1(1) * (point_2(2) - point_3(2));
double_area_B = point_2(1) * (point_3(2) - point_1(2));
double_area_C = point_3(1) * (point_1(2) - point_2(2));
area = 0.5 * abs(double_area_A + double_area_B + double_area_C);

end