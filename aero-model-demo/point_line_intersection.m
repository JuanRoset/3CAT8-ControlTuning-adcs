function point = point_line_intersection(line_point, line_vector, plane_point, plane_normal)

% Normalize vectors
line_vector = line_vector / norm(line_vector);
plane_normal = plane_normal / norm(plane_normal);

% Calculate the plane's offset 'd'
d = - dot(plane_point, plane_normal);

% Check if the line is parallel to the plane
if dot(plane_normal, line_vector) == 0
    error('Line is parallel to the plane. No intersection.');
end

% Calculate the line parameter and the resulting point
t = - (dot(plane_normal, line_point) + d) / dot(plane_normal, line_vector);
point = line_point + t .* line_vector;

end