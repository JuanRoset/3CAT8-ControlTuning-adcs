clear
close all
clc

% Set dimensions
width = 0.1;
height = 0.2;
long = 0.3;
inner_radius = 0.225;
outer_radius = 0.6;

% Create vortexes
vertexs = [-long/2, -width/2,  height/2;
           -long/2,  width/2,  height/2;
           -long/2,  width/2, -height/2;
           -long/2, -width/2, -height/2;
            long/2, -width/2,  height/2;
            long/2,  width/2,  height/2;
            long/2,  width/2, -height/2;
            long/2, -width/2, -height/2];

[n_points, ~] = size(vertexs);
            
% Define surfaces
surfaces = [1, 2, 3;
            3, 4, 1;
            5, 6, 7;
            7, 8, 5;
            2, 6, 7;
            7, 3, 2;
            7, 3, 4;
            4, 8, 7;
            8, 4, 1;
            1, 5, 8;
            1, 5, 6;
            6, 2, 1];

% Define ray
ray = [1; 1; 1];
ray = ray ./ norm(ray);

% Define new coordinates
if abs(ray(3)) == 1
    ex = sign(ray) .* [1; 0; 0];
    ey = sign(ray) .* [0; 1; 0];
else 
    ek = [0; 0; 1];
    ex = cross(ray, ek);
    ey = cross(ex, ray);
end
E = [ex, ey, ray];
E_warn = inv(E); % Double assign no to get matlab's orange warning >:/
E_inv = E_warn;

[n_surfaces, n_sides] = size(surfaces);

%project points
projected_points = zeros(n_points, 3);
points_2D = zeros(n_points, 2);
centroids_2D = zeros(n_surfaces, 2);

% Project points from 3D space to 2D coordinates
for i = 1:n_points
    projected_points(i, :) = (vertexs(i, :)' - ray * (ray' * vertexs(i, :)'))';
    twoD_point = E_inv * projected_points(i, :)';
    points_2D(i, :) = [twoD_point(1), twoD_point(2)];
end

% Rasterization parameters
delta_x = 0.01;
delta_y = 0.01;
delta_S = delta_x * delta_y;

% Find min and max discretized coordinates
x_minimum = min(points_2D(:, 1));
x_maximum = max(points_2D(:, 1));
y_minimum = min(points_2D(:, 2));
y_maximum = max(points_2D(:, 2));
x_maximum = delta_x * (floor((x_maximum - x_minimum) / delta_x) + 1) + x_minimum;
y_maximum = delta_y * (floor((y_maximum - y_minimum) / delta_y) + 1) + y_minimum;

% Determine number of cells per axis
n_cells_x = (x_maximum - x_minimum) / delta_x;
n_cells_y = (y_maximum - y_minimum) / delta_y;

% Determine positions of cell centers
p_x = linspace(x_minimum + delta_x / 2, x_maximum - delta_x / 2, n_cells_x);
p_y = linspace(y_minimum + delta_y / 2, y_maximum - delta_y / 2, n_cells_y);

segments = [];

% Find centroids of all surfaces in 2D
for i = 1:n_surfaces
    indices = surfaces(i, :);
    
    % Circular shift the indices array to include the wrap-around edge
    x_points = points_2D(indices, 1)';
    y_points= points_2D(indices, 2)';
    
    Area = polygon_area(x_points, y_points);
    [x_bar, y_bar] = polygon_centroid(x_points, y_points, Area);

    centroids_2D(i, :) = [x_bar, y_bar];
end

% Find the 2D segments that make up the projected mesh
for i = 1:n_surfaces
    indices = surfaces(i, :);
    
    % Circular shift the indices array to include the wrap-around edge
    indices_shifted = [indices, indices(1)];
    
    % Extract the x and y coordinates of points using indices
    x_coords = points_2D(indices_shifted, 1);
    y_coords = points_2D(indices_shifted, 2);
    
    % Add all of the segments
    for j = 1:(length(x_coords) - 1)
        current_segment = [x_coords(j), y_coords(j), x_coords(j + 1), y_coords(j + 1)];
        segments = cat(1, segments, current_segment);
    end    
end

% Remove duplicated element from array
segments = unique(segments, 'rows');

% Print segments
figure()
hold on
for i = 1:n_points
    scatter(points_2D(i, 1), points_2D(i, 2))
end

figure()
hold on
colors = lines(n_surfaces); % Generate a set of unique colors
for i = 1:n_surfaces
    indices = surfaces(i, :);
    
    % Circular shift the indices array to include the wrap-around edge
    indices_shifted = [indices, indices(1)];
    
    % Extract the x and y coordinates of points using indices
    x_coords = points_2D(indices_shifted, 1);
    y_coords = points_2D(indices_shifted, 2);
    
    % Plot the surface edges with unique color
    %plot(x_coords', y_coords', '-', 'Color', colors(i, :));

    plot(x_coords', y_coords', '-', 'Color', 'k');

    % Scatter plot with unique color for each surface
    scatter(centroids_2D(i, 1), centroids_2D(i, 2), 10, 'b', 'filled');
end

%hold on
%for i = 1:n_cells_x
%    for j = 1:n_cells_y
%        scatter(p_x(i), p_y(j))
%    end
%end

% Plot vertices
figure()
scatter3(vertexs(:,1), vertexs(:,2), vertexs(:,3), 'filled');
hold on;

% Plot edges
for i = 1:size(surfaces, 1)
    edge_points = vertexs(surfaces(i,:), :);
    plot3(edge_points(:,1), edge_points(:,2), edge_points(:,3), 'k', 'LineWidth', 1.5);
end

% Initialize points that interact with the medium
direct_points = [];
direct_surfaces = [];
direct_normals = [];
direct_areas = [];

% Loop through rasterized points
for i = 1:n_cells_x
    for j = 1:n_cells_y
        % Fetch point
        point = [p_x(i), p_y(j)];
        possible_surfaces = [];
        possible_points = [];
        possible_distances = [];
        possible_normals = [];
        possible_areas = [];

        % Check intersections with surfaces
        for k= 1:n_surfaces
            indices = surfaces(k, :);

            % Extract the x and y coordinates of points using indices
            triangle = points_2D(indices, :);
            x_coords = points_2D(indices, 1);
            y_coords = points_2D(indices, 2);
            
            % Check if the point is inside a projected surface
            if point_in_triangle(triangle, point)
                possible_surfaces = cat(1, possible_surfaces, k);

                % Convert 2D coordinates to 3D points in projection plane
                ray_point = E * [point, 0]';

                % Compute the plane normal
                triangle_points = vertexs(indices, :);
                triangle_vectors = [triangle_points(3, :) - triangle_points(1, :);
                                    triangle_points(2, :) - triangle_points(1, :)];
                plane_normal = cross(triangle_vectors(1, :), triangle_vectors(2, :));
                plane_normal = plane_normal / norm(plane_normal);
                %plane_normal = surface_normals(k, :);
                possible_normals = cat(1, possible_normals, plane_normal);

                % Compute the intersection point
                possible_point = point_line_intersection(ray_point, ray, triangle_points(1, :), plane_normal);
                possible_points = cat(1, possible_points, possible_point');

                % Compute the distance to the plane
                changed_point = E_inv * possible_point;
                distance = changed_point(3);
                possible_distances = cat(1, possible_distances, distance);

                % Compute the unprojected concentrated areas
                ray_unit = - ray / norm(ray);
                cos_theta = dot(ray_unit, plane_normal) / norm(plane_normal);
                unprojected_area = delta_S / abs(cos_theta);
                possible_areas = cat(1, possible_areas, unprojected_area);
            end
        end

        % 

        if ~isempty(possible_surfaces)
            [~, first_index] = min(possible_distances);
            final_point = possible_points(first_index, :);
            final_surface = surfaces(possible_surfaces(first_index, :));
            final_normal = possible_normals(first_index, :);
            final_area = possible_areas(first_index);

            direct_points = cat(1, direct_points, final_point);
            direct_normals = cat(1, direct_normals, final_normal);
            direct_areas = cat(1, direct_areas, final_area);
            direct_surfaces = cat(1, direct_surfaces, final_surface);

            scatter3(final_point(1), final_point(2), final_point(3), 'filled', 'MarkerFaceColor',[0 .75 .75]);
        end

        %disp(possible_surfaces);

    end
end

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Ray intersections');

% Set aspect ratio to equal
axis equal;

%%

% Close all existing figures
%close all

% Define velocity
Velocity = 7800;

% Define CG
position_CG = [long/2, 0, 0];

% Define velocity vector
V = Velocity * ray';

% Define constants
rho = 6.03e-14;
sigma_n = 0.8;
sigma_t = 0.8;

% Initialize arrays to store aerodynamic coefficients
%CL = zeros(length(beta_array), 1);
%CD = zeros(length(beta_array), 1);

% initialize total force coefficients
Force = [0, 0, 0];
Torque = [0, 0, 0];

% Iterate through each beta value
for i = 1:size(direct_points, 1)
    % Extract current beta
    normal = direct_normals(i, :);
    normal = normal / norm(normal);
    
    % Compute area vector based on current beta
    A = normal .* direct_areas(i)';
    
    % Normalize velocity and area vectors
    V_unit = (V ./ Velocity)';
    A_unit = normal';
    
    % Compute aerodynamic force components
    F = - rho * (V' * V) * ((2 - sigma_n - sigma_t) * (V_unit' * A_unit)^2 * A_unit + sigma_t * (V_unit' * A_unit) * V_unit) * direct_areas(i);
    %Lift = F(2);
    %Drag = -F(1);

    lever_arm = direct_points(i) - position_CG;

    Force = Force + F';
    Torque = Torque + cross(lever_arm, F);
    
    % Compute aerodynamic coefficients and store them
    %CL(i) = 2 * Lift / (rho * Surface * Velocity^2);
    %CD(i) = -2 * Drag / (rho * Surface * Velocity^2);
end

disp(Force)
disp(Torque)