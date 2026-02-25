% Set dimensions
body_width = 0.1;                   % Width of the body
body_height = 0.2;                  % Height of the body
body_length = 0.3;                  % Length of the body

panel_length = 0.25;                % Length of the solar panel
panel_height = 0.175;               % Height of the solar panel

inner_radius = 0.323;               % Inner radius of the antenna
outer_radius = 1.152;                % Outer radius of the antenna
antenna_separation = body_length / 2 + 0.5; % Distance from the origin to antenna
antenna_phase = (60 + 7.55) * pi / 180;     % Phase angle of the antenna

antenna_lowPoly_correction = 0.5103; % Correct for the more detailed antenna geometry

% Specify surface radiation properties
rho_a = 0.910;
rho_s = 0.072;
rho_d = 0.018;

% Specify surface aerosynamic properties
sigma_n = 0.8;
sigma_t = 0.8;

% Create body vertexes
vertexs = [ body_height/2, -body_width/2, -body_length/2;
            body_height/2,  body_width/2, -body_length/2;
           -body_height/2,  body_width/2, -body_length/2;
           -body_height/2, -body_width/2, -body_length/2;
            body_height/2, -body_width/2,  body_length/2;
            body_height/2,  body_width/2,  body_length/2;
           -body_height/2,  body_width/2,  body_length/2;
           -body_height/2, -body_width/2   body_length/2];  % Vertices of the body

% Define body surfaces
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
            6, 2, 1];    % Faces of the body

% Define body surface normals
surface_normals = [ 0,  0, -1;
                    0,  0, -1;
                    0,  0,  1;
                    0,  0,  1;
                    0,  1,  0;
                    0,  1,  0;
                   -1,  0,  0;
                   -1,  0,  0;
                    0, -1,  0;
                    0, -1,  0;
                    1,  0,  0;
                    1,  0,  0];    % Normals of the body surfaces

% Add surface properties
surface_correction = ones(size(surfaces, 1), 1);
surface_rho_a = rho_a .* ones(size(surfaces, 1), 1);
surface_rho_s = rho_s .* ones(size(surfaces, 1), 1);
surface_rho_d = rho_d .* ones(size(surfaces, 1), 1);
surface_sigma_n = sigma_n .* ones(size(surfaces, 1), 1);
surface_sigma_t = sigma_t .* ones(size(surfaces, 1), 1);

% Add deployable solar panel
panel_surfaces = [1, 2, 3;
                  2, 3, 4];    % Faces of the solar panel
panel_surfaces = panel_surfaces + size(vertexs, 1); % Update surface indices
surfaces = cat(1, surfaces, panel_surfaces);    % Add panel faces to surfaces
surfaces = cat(1, surfaces, panel_surfaces);    % Add panel faces again to surfaces

% Update surface properties
surface_correction = cat(1, surface_correction, ones(size(panel_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_rho_a = cat(1, surface_rho_a, rho_a .* ones(size(panel_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_rho_s = cat(1, surface_rho_s, rho_s .* ones(size(panel_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_rho_d = cat(1, surface_rho_d, rho_d .* ones(size(panel_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_sigma_n = cat(1, surface_sigma_n, sigma_n .* ones(size(panel_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_sigma_t = cat(1, surface_sigma_t, sigma_t .* ones(size(panel_surfaces, 1) * 2, 1));    % Add panel faces to surfaces

panel_vertexs = [-body_height/2,  body_width/2,  panel_length/2;
                 -body_height/2,  body_width/2, -panel_length/2;
                 -body_height/2 - panel_height,  body_width/2,  panel_length/2;
                 -body_height/2 - panel_height,  body_width/2, -panel_length/2];    % Vertices of the solar panel
vertexs = cat(1, vertexs, panel_vertexs);    % Add panel vertices to vertexs
surface_normals = cat(1, surface_normals, repmat([0, -1, 0], size(panel_surfaces, 1), 1));    % Update surface normals
surface_normals = cat(1, surface_normals, repmat([0,  1, 0], size(panel_surfaces, 1), 1));    % Update surface normals

% Add deployable antenna
antenna_surfaces = [1, 9, 4;
                    1, 4, 5;
                    1, 5, 2;
                    2, 5, 6;
                    2, 6, 7;
                    2, 7, 3;
                    3, 7, 8;
                    3, 8, 9;
                    3, 9, 1];    % Faces of the antenna
antenna_surfaces = antenna_surfaces + size(vertexs, 1);    % Update surface indices

surfaces = cat(1, surfaces, antenna_surfaces);    % Add antenna faces to surfaces
surfaces = cat(1, surfaces, antenna_surfaces);    % Add antenna faces again to surfaces

% Update surface properties
surface_correction = cat(1, surface_correction, antenna_lowPoly_correction .* ones(size(antenna_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_rho_a = cat(1, surface_rho_a, rho_a .* ones(size(antenna_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_rho_s = cat(1, surface_rho_s, rho_s .* ones(size(antenna_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_rho_d = cat(1, surface_rho_d, rho_d .* ones(size(antenna_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_sigma_n = cat(1, surface_sigma_n, sigma_n .* ones(size(antenna_surfaces, 1) * 2, 1));    % Add panel faces to surfaces
surface_sigma_t = cat(1, surface_sigma_t, sigma_t .* ones(size(antenna_surfaces, 1) * 2, 1));    % Add panel faces to surfaces

surface_normals = cat(1, surface_normals, repmat([0, 0, 1], size(antenna_surfaces, 1), 1));    % Update surface normals
surface_normals = cat(1, surface_normals, repmat([0, 0, -1], size(antenna_surfaces, 1), 1));    % Update surface normals

% Calculate vertices for outer antenna
outer_antenna_vertexs = [];
for i = 1:3
    new_point = [outer_radius * cos(antenna_phase + (i - 1) * 2 / 3 * pi), outer_radius * sin(antenna_phase + (i - 1) * 2 / 3 * pi), antenna_separation];
    outer_antenna_vertexs = cat(1, outer_antenna_vertexs, new_point);
end
vertexs = cat(1, vertexs, outer_antenna_vertexs);    % Add outer antenna vertices to vertexs

% Calculate vertices for inner antenna
inner_antenna_vertexs = [];
for i = 1:6
    new_point = [inner_radius * cos(antenna_phase + (i - 1) * pi / 3 + pi / 6), inner_radius * sin(antenna_phase + (i - 1) * pi / 3 + pi / 6), antenna_separation];
    inner_antenna_vertexs = cat(1, inner_antenna_vertexs, new_point);
end
vertexs = cat(1, vertexs, inner_antenna_vertexs);    % Add inner antenna vertices to vertexs

% Calculate surface centers and values
surface_values = zeros(size(surfaces,1),1);   % preallocate
surface_centers = zeros(size(surfaces,1),3);   % preallocate
for i = 1:size(surfaces,1)
    vertex_idxs = surfaces(i, :);
    a = vertexs(vertex_idxs(1), :);
    b = vertexs(vertex_idxs(2), :);
    c = vertexs(vertex_idxs(3), :);
    surface_centers(i,:) = (a + b + c)/3;

    % --- centroid ---
    surface_centers(i,:) = (a + b + c)/3;
    
    % --- triangle area ---
    ab = b - a;
    ac = c - a;
    surface_values(i) = 0.5 * norm(cross(ab, ac));
end

% Plotting
close all
hfig = figure();
for i = 1:size(surfaces, 1)
    indices = [surfaces(i, :), surfaces(i, 1)];
    if surface_correction(i) == 1
        color = [0.2588    0.3961    0.4902];
    else
        color = [0.9294    0.6941    0.1255];
    end
    fill3(vertexs(indices, 1), vertexs(indices, 2), vertexs(indices, 3), color)    % Fill surface with red color
    hold on

    for j = 1:length(indices) - 1
        point_a = vertexs(indices(j), :);
        point_b = vertexs(indices(j + 1), :);

        plot3([point_a(1), point_b(1)], [point_a(2), point_b(2)], [point_a(3), point_b(3)], '-k', 'LineWidth', 1.0)    % Plot edges
    end
end
antenna_middle = [-body_height/4, 0, body_length/2];
plot3([vertexs(13, 1), antenna_middle(1)], [vertexs(13, 2), antenna_middle(2)], [vertexs(13, 3), antenna_middle(3)], '-k', 'LineWidth', 1.0)    % Plot edges
plot3([vertexs(14, 1), antenna_middle(1)], [vertexs(14, 2), antenna_middle(2)], [vertexs(14, 3), antenna_middle(3)], '-k', 'LineWidth', 1.0)    % Plot edges
plot3([vertexs(15, 1), antenna_middle(1)], [vertexs(15, 2), antenna_middle(2)], [vertexs(15, 3), antenna_middle(3)], '-k', 'LineWidth', 1.0)    % Plot edges
view([180 + 45, 17])

%xlabel('X-axis (m)');    % Label x-axis
%ylabel('Y-axis (m)');    % Label y-axis
%zlabel('Z-axis (m)');    % Label z-axis
%title('Simplified body geometry');    % Set title

% ---- Plot surface centers ----
scatter3(surface_centers(:,1), ...
         surface_centers(:,2), ...
         surface_centers(:,3), ...
         18, 'r', 'filled', ...
         'MarkerEdgeColor','k');

set(gca,'XColor', 'none','YColor','none','ZColor','none')
set(gca, 'color', 'none');
axis equal;    % Set equal scale for all axes

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',12) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

save('3cat8_geometry.mat', 'vertexs', 'surfaces', 'surface_values', 'surface_centers', 'surface_normals', 'surface_correction', 'surface_rho_a', 'surface_rho_s', 'surface_rho_d', 'surface_sigma_n', 'surface_sigma_t')    % Save geometry to file