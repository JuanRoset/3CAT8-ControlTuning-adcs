clear
close all
clc

% Load satellite geometry
load('3cat8_geometry.mat')

% Define CG
position_CG = [0.0, 0.0, 0.0];

% Specify aerodynamic properties
Velocity = 7672;
rho = 3.8e-12;

% Define sunlight parameters
c = 2.99792458e8; % Speed of light in vacuum
Ls = 3.828*10^26; % Solar luminosity constant
AU = 149597870.7*10^3; % Astronomical unit
P = Ls / (4 * pi * AU^2) / c; % Pressure constant

% Define rasterization parameters
delta_x = 0.02;
delta_y = 0.02;

%% Generate rays systematically
yaw = linspace(0, 180, 12) .* pi / 180;
yaw(1)   = [];
yaw(end) = [];
roll = linspace(0, 360, 11) .* pi / 180;
roll(end) = [];

ray = [0; 0; 1]';
ray_array = zeros(length(yaw) * length(roll) + 2, 3);
ray_array(1,   :) =  ray;
ray_array(end, :) = -ray;

for i = 1:length(yaw)

    yawed_ray = euler2rot(yaw(i), 0, 0) * ray';

    for j = 1:length(roll)

        counter = (i - 1) * length(roll) + j;
        final_ray = euler2rot(0, 0, roll(j)) * yawed_ray;
        ray_array(counter + 1, :) = final_ray';
        scatter3(ray_array(counter + 1, 1), ray_array(counter + 1, 2), ray_array(counter + 1, 3), 'b.', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

    end
end

figure()
scatter3(ray_array(:, 1), ray_array(:, 2), ray_array(:, 3), 'b.', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
axis equal;

%% Generate rays randomly
n_rays = 100;
long_lat = rand([n_rays, 2]);
long_lat(:, 1) = long_lat(:, 1) .* 2 * pi;
long_lat(:, 2) = asin((long_lat(:, 2) - 0.5) .* 2);
%ray_array = [cos(long_lat(:, 2)) .* [cos(long_lat(:, 1)), sin(long_lat(:, 1))], sin(long_lat(:, 2))];

%figure()
%scatter(long_lat(:, 1), long_lat(:, 2), 'b.', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

%figure()
%scatter3(ray_array(:, 1), ray_array(:, 2), ray_array(:, 3), 'b.', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%axis equal;

%% Perform numerical approximation
Forces_srp = zeros(size(ray_array, 1), 3);
Forces_air = zeros(size(ray_array, 1), 3);
Torques_srp = zeros(size(ray_array, 1), 3);
Torques_air = zeros(size(ray_array, 1), 3);

for a = 1:size(ray_array, 1)

    % Print
    disp(['Ray number: ', num2str(a), '    -    Progress: ', num2str(100 * a / size(ray_array, 1)), ' %'])

    % Store current ray
    ray = ray_array(a, :)';

    % Define new coordinates
    if abs(ray(3)) == 1
        ex = sign(ray(3)) .* [1; 0; 0];
        ey = sign(ray(3)) .* [0; 1; 0];
    else 
        ek = [0; 0; 1];
        ex = cross(ray, ek);
        ey = cross(ex, ray);
    end
    E = [ex, ey, ray];
    E_warn = inv(E); % Double assign no to get matlab's orange warning >:/
    E_inv = E_warn;
    
    [n_points, ~] = size(vertexs);
    [n_surfaces, n_sides] = size(surfaces);
    
    %project points
    projected_points = zeros(n_points, 3);
    points_2D = zeros(n_points, 2);
    
    % Project points from 3D space to 2D coordinates
    for i = 1:n_points
        projected_points(i, :) = (vertexs(i, :)' - ray * (ray' * vertexs(i, :)'))';
        twoD_point = E_inv * projected_points(i, :)';
        points_2D(i, :) = [twoD_point(1), twoD_point(2)];
    end
    
    % Find maximum and minimum x-y coordinates of projected points
    x_minimum = min(points_2D(:, 1));
    x_maximum = max(points_2D(:, 1));
    y_minimum = min(points_2D(:, 2));
    y_maximum = max(points_2D(:, 2));
    x_maximum = delta_x * (floor((x_maximum - x_minimum) / delta_x) + 1) + x_minimum;
    y_maximum = delta_y * (floor((y_maximum - y_minimum) / delta_y) + 1) + y_minimum;
    
    % Define the number of cells per axis
    n_cells_x = (x_maximum - x_minimum) / delta_x;
    n_cells_y = (y_maximum - y_minimum) / delta_y;
    
    % Set the positions of the centers of each cell
    p_x = linspace(x_minimum + delta_x / 2, x_maximum - delta_x / 2, n_cells_x);
    p_y = linspace(y_minimum + delta_y / 2, y_maximum - delta_y / 2, n_cells_y);
    
    % initialize total force and torques as zero
    Force_srp = [0, 0, 0];
    Force_air = [0, 0, 0];
    Torque_srp = [0, 0, 0];
    Torque_air = [0, 0, 0];
    Total_area = 0;
    
    % Loop through rasterized points
    for i = 1:n_cells_x
        for j = 1:n_cells_y
            
            % Fetch point
            point = [p_x(i), p_y(j)];
            
            % Initialize the found intersection points
            possible_surfaces = [];
            possible_points = [];
            possible_distances = [];
            possible_normals = [];
    
            % Check intersections with surfaces
            for k= 1:n_surfaces
                indices = surfaces(k, :);
    
                % Extract the x and y coordinates of points using indices
                triangle = points_2D(indices, :);
                x_coords = points_2D(indices, 1);
                y_coords = points_2D(indices, 2);
    
                % Check if ray intersects surface and filter via cosines
                plane_normal = surface_normals(k, :);
                coincident = false;
                if dot(plane_normal, ray) < 0
                    coincident = point_in_triangle(triangle, point);
                end
                
                % Check if the point is inside a projected surface
                if coincident
    
                    % Store normal and surface points
                    triangle_points = vertexs(indices, :);
                    possible_normals = cat(1, possible_normals, plane_normal);
    
                    % Store index of possible surfaces
                    possible_surfaces = cat(1, possible_surfaces, k);
    
                    % Convert 2D coordinates to 3D points in projection plane
                    ray_point = E * [point, 0]';
    
                    % Compute the intersection point
                    possible_point = point_line_intersection(ray_point, ray, triangle_points(1, :), plane_normal);
                    possible_points = cat(1, possible_points, possible_point');
    
                    % Compute the distance to the plane
                    changed_point = E_inv * possible_point;
                    distance = changed_point(3);
                    possible_distances = cat(1, possible_distances, distance);
    
                    
                end
            end
    
            % If there was one or more intersecting surfaces, store first one
            if ~isempty(possible_surfaces)
    
                % Find the smallest distance as the first one to collide
                [~, first_index] = min(possible_distances);
                final_point = possible_points(first_index, :);
                surface_idx = possible_surfaces(first_index, :);
                final_normal = possible_normals(first_index, :);
    
                % Unproject areas
                plane_normal = surface_normals(surface_idx, :);
                
                % Find edge points of cell in 2D plane
                point_a = [point - [delta_x, delta_y] ./ 2, 0];
                point_b = point_a + [delta_x, 0, 0];
                point_c = point_a + [0, delta_y, 0];
                
                % Find edge points of cell to 3D plane
                point_a = E * point_a';
                point_b = E * point_b';
                point_c = E * point_c';
                
                % Re-project edje points onto surface
                point_a = point_line_intersection(point_a, ray, final_point', plane_normal');
                point_b = point_line_intersection(point_b, ray, final_point', plane_normal');
                point_c = point_line_intersection(point_c, ray, final_point', plane_normal');
    
                % Compute real area of cell
                vec_a = point_b - point_a;
                vec_b = point_c - point_a;
                final_area = norm(cross(vec_a, vec_b));
    
                % Extract current surface normal
                normal = plane_normal;
                
                % Normalize velocity and area vectors
                V_unit = ray;
                A_unit = normal';
                Area = final_area * surface_correction(surface_idx);
    
                rho_a = surface_rho_a(surface_idx);
                rho_s = surface_rho_s(surface_idx);
                rho_d = surface_rho_d(surface_idx);
    
                sigma_n = surface_sigma_n(surface_idx);
                sigma_t = surface_sigma_t(surface_idx);
                
                % Compute aerodynamic force
                F_srp = - P * Area * dot(A_unit, ray) .* ((1 - rho_s) * ray + 2 * (rho_s * dot(A_unit, ray) - rho_d / 3) * A_unit);
                F_air = - rho * Velocity^2 * Area .* ((2 - sigma_n - sigma_t) * (V_unit' * A_unit)^2 * A_unit + sigma_t * (V_unit' * A_unit) * V_unit);
            
                % Compute surface differential position w.r.t. the CG
                lever_arm = final_point - position_CG;
                
                % Add the force and torque contributions
                Force_srp = Force_srp + F_srp';
                Force_air = Force_air + F_air';
                Torque_srp = Torque_srp + cross(lever_arm, F_srp);
                Torque_air = Torque_air + cross(lever_arm, F_air);

                % Store the forces and torques
                Forces_srp(a, :) = Force_srp;
                Forces_air(a, :) = Force_air;
                Torques_srp(a, :) = Torque_srp;
                Torques_air(a, :) = Torque_air;
    
                %Total_area = Total_area + Area;
    
                
            end
    
        end
    end

end

save('train_data_100.mat', 'ray_array', 'Forces_srp', 'Forces_air', 'Torques_srp', 'Torques_air')