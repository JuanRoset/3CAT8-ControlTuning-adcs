clear
close all
clc

% Define rasterization parameters
delta_array = [0.5, 0.05, 0.005, 0.0005];
delta_array = expspace(0.25, 0.025, 30)';

% Define velocity ray
ray = [1; 1; 1];
ray = ray ./ norm(ray);

% Define aerodynamic parameters
Velocity = 7672;
rho = 3.8e-12;
sigma_n = 0.8;
sigma_t = 0.8;

% Define CG
position_CG = [0.05, 0.0025, 0.05];

%% Perform analytical computation of the forces and torques in the simple case
area_a = 0.03;
area_b = 0.06;
area_c = 0.02;

normal_a = [-1, 0, 0]';
normal_b = [0, -1, 0]';
normal_c = [0, 0, -1]';

p_a = [-0.1,  0, 0];
p_b = [0, -0.05, 0];
p_c = [0, 0, -0.15];

% Compute aerodynamic force components
F_a = - rho * Velocity^2 * area_a .* ((2 - sigma_n - sigma_t) * (ray' * normal_a)^2 * normal_a + sigma_t * (ray' * normal_a) * ray);
F_b = - rho * Velocity^2 * area_b .* ((2 - sigma_n - sigma_t) * (ray' * normal_b)^2 * normal_b + sigma_t * (ray' * normal_b) * ray);
F_c = - rho * Velocity^2 * area_c .* ((2 - sigma_n - sigma_t) * (ray' * normal_c)^2 * normal_c + sigma_t * (ray' * normal_c) * ray);

% Compute aerodynamic torques of all plates
r_a = p_a - position_CG;
r_b = p_b - position_CG;
r_c = p_c - position_CG;
tau_a = cross(r_a, F_a);
tau_b = cross(r_b, F_b);
tau_c = cross(r_c, F_c);

% Compute total forces and torques
F_analytical = (F_a + F_c + F_b)';
T_analytical = tau_a + tau_b + tau_c;

%% Perform numerical approximation
load('3cat8_geometry.mat')

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

approx_forces = zeros(size(delta_array, 1), 3);
approx_torques = zeros(size(delta_array, 1), 3);
num_cells = zeros(size(delta_array, 1), 1);

for n = 1:size(delta_array, 1)

    % Define rasterization parameters
    delta_x = delta_array(n);
    delta_y = delta_x;

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
                
                % Check if the point is inside a projected surface
                if point_in_triangle(triangle, point)
                    % Store index of possible surfaces
                    possible_surfaces = cat(1, possible_surfaces, k);
    
                    % Convert 2D coordinates to 3D points in projection plane
                    ray_point = E * [point, 0]';
    
                    % Compute the plane normal
                    triangle_points = vertexs(indices, :);
                    %triangle_vectors = [triangle_points(3, :) - triangle_points(1, :);
                    %                    triangle_points(2, :) - triangle_points(1, :)];
                    %plane_normal = cross(triangle_vectors(1, :), triangle_vectors(2, :));
                    plane_normal = surface_normals(k, :);
                    possible_normals = cat(1, possible_normals, plane_normal);
    
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
                %point_d = point_a + [delta_x, delta_y, 0];
                
                % Find edge points of cell to 3D plane
                point_a = E * point_a';
                point_b = E * point_b';
                point_c = E * point_c';
                %point_d = E * point_d';
                
                % Re-project edje points onto surface
                point_a = point_line_intersection(point_a, ray, final_point', plane_normal');
                point_b = point_line_intersection(point_b, ray, final_point', plane_normal');
                point_c = point_line_intersection(point_c, ray, final_point', plane_normal');
                %point_d = point_line_intersection(point_d, ray, final_point', plane_normal');
    
                % Compute real area of cell
                vec_a = point_b - point_a;
                vec_b = point_c - point_a;
                final_area = norm(cross(vec_a, vec_b));
    
                % Store relevant found parameters
                direct_points = cat(1, direct_points, final_point);
                direct_normals = cat(1, direct_normals, final_normal);
                direct_areas = cat(1, direct_areas, final_area);
                direct_surfaces = cat(1, direct_surfaces, possible_surfaces(first_index, :));
    
                
            end
    
        end
    end

    % Approximate Aerodynamic calculation
    
    % Define velocity vector
    V = Velocity * ray';
    
    % initialize total force and torques as zero
    Force = [0, 0, 0];
    Torque = [0, 0, 0];
    
    % Iterate through each beta value
    for i = 1:size(direct_points, 1)
    
        % Extract current surface normal
        normal = direct_normals(i, :);
        
        % Normalize velocity and area vectors
        V_unit = ray;
        A_unit = normal';
        
        % Compute aerodynamic force
        F = - rho * Velocity^2 * direct_areas(i) .* ((2 - sigma_n - sigma_t) * (V_unit' * A_unit)^2 * A_unit + sigma_t * (V_unit' * A_unit) * V_unit);
    
        % Compute surface differential position w.r.t. the CG
        lever_arm = direct_points(i, :) - position_CG;
        
        % Add the force and torque contributions
        Force = Force + F';
        Torque = Torque + cross(lever_arm, F);
    end

    approx_forces(n, :) = Force;
    approx_torques(n, :) = Torque;
    num_cells(n) = n_cells_x * n_cells_y;

end

%% Process the results

close all
norm_forces_approx = sqrt(sum(approx_forces.^2, 2));
norm_torques_approx = sqrt(sum(approx_torques.^2, 2));
force_errors = abs(norm_forces_approx - norm(F_analytical)) ./ norm(F_analytical);
torque_errors = abs(norm_torques_approx - norm(T_analytical)) ./ norm(T_analytical);
log_Ncells = log10(num_cells);

% Find error model
log_Ferr = log10(100 .* force_errors);
log_Terr = log10(100 .* torque_errors);
Ferr_model = polyfit(log_Ncells, log_Ferr, 1);
Terr_model = polyfit(log_Ncells, log_Terr, 1);
log_Xcoords = linspace(min(log_Ncells), max(log_Ncells));

% Plot the obtained magnitudes

% Define colors
cool_red = [0.7020    0.2941    0.4039];
cool_green = [0.5216    0.7412    0.4392];
cool_yellow = [1.0000    0.8510    0.5020];
cool_blue = [0.4510    0.7412    0.8314];
cool_purple = [0.9020    0.5020    0.9490];

line_width = 1.0;

% Plot results and errors
hfig = figure('Position',[100, 100, 1200, 800]);
tiledlayout(2, 2)

nexttile
loglog(num_cells, norm_forces_approx, 'k', 'LineWidth', line_width, 'DisplayName', 'Approximate');
hold on
loglog([min(num_cells), max(num_cells)], norm(F_analytical) .* [1, 1], ':k', 'LineWidth', line_width, 'DisplayName', 'Exact')
%grid on
xlim([min(num_cells), max(num_cells)])
ylabel('Force modulus (N)')
xlabel('Number of cells')
legend('Location','northeast')
set(findall(hfig,'-property','FontSize'),'FontSize',12)
set(findall(hfig,'-property','Box'),'Box','off')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

nexttile
loglog(num_cells, 100 .* force_errors, 'k', 'LineWidth', line_width, 'DisplayName', 'Real error');
hold on
plot(10.^log_Xcoords, 10.^(Ferr_model(2) + Ferr_model(1) .* log_Xcoords), '-.k', 'LineWidth', line_width, 'DisplayName', 'Model')
%grid on
xlim([min(num_cells), max(num_cells)])
ylim(100 .* [min(force_errors), max(force_errors)])
ylabel('Force relative error ($\%$)')
xlabel('Number of cells')
legend('Location','northeast')
set(findall(hfig,'-property','FontSize'),'FontSize',12)
set(findall(hfig,'-property','Box'),'Box','off')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

nexttile
loglog(num_cells, norm_torques_approx, 'k', 'LineWidth', line_width, 'DisplayName', 'Approximate');
hold on
loglog([min(num_cells), max(num_cells)], norm(T_analytical) .* [1, 1], ':k', 'LineWidth', line_width, 'DisplayName', 'Exact')
%grid on
xlim([min(num_cells), max(num_cells)])
ylabel('Torque modulus (Nm)')
xlabel('Number of cells')
legend('Location','northeast')
set(findall(hfig,'-property','FontSize'),'FontSize',12)
set(findall(hfig,'-property','Box'),'Box','off')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

nexttile
loglog(num_cells, 100 .* torque_errors, 'k', 'LineWidth', line_width, 'DisplayName', 'Real error');
hold on 
plot(10.^log_Xcoords, 10.^(Terr_model(2) + Terr_model(1) .* log_Xcoords), '-.k', 'LineWidth', line_width, 'DisplayName', 'Model')
%grid on
xlim([min(num_cells), max(num_cells)])
ylim(100 .* [min(torque_errors), max(torque_errors)])
ylabel('Torque relative error ($\%$)')
xlabel('Number of cells')
legend('Location','northeast')
set(findall(hfig,'-property','FontSize'),'FontSize',12)
set(findall(hfig,'-property','Box'),'Box','off')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%% Find desired number of cells
desired_error = 2.5; % In '%'
n_cells = floor(10^((log10(desired_error )- Terr_model(2)) / Terr_model(1))) + 1;

clc
disp('- Convergence results -')
disp(['For a torque relative error of: ', num2str(desired_error), ' %'])
disp(['A number of: ', num2str(n_cells), ' cells is necessary'])