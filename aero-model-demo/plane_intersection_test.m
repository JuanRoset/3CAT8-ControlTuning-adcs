% Define the line (point and direction vector)
line_point = [1, 1, 1];
line_vector = [1, 2, 3];

% Define the plane (point and normal vector)
plane_point = [4, 4, 0];
plane_normal = [1, 1, 1];

% Calculate the plane's offset 'd'
d = -dot(plane_normal, plane_point);

% Calculate the intersection point
intersection_point = point_line_intersection(line_point, line_vector, plane_point, plane_normal);

% Visualization
figure;
hold on;

% Plot the line
line_end = line_point + 10 * line_vector; % Extend the line for visualization
plot3([line_point(1), line_end(1)], ...
      [line_point(2), line_end(2)], ...
      [line_point(3), line_end(3)], 'b');

% Plot the plane
[X,Y] = meshgrid(-5:5,-5:5);
Z = -(plane_normal(1)*X + plane_normal(2)*Y + d) / plane_normal(3);
surf(X,Y,Z,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','r');

% Plot the intersection point
plot3(intersection_point(1), intersection_point(2), intersection_point(3), 'ro', 'MarkerSize', 10);

% Set labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Intersection of Line and Plane');

axis equal; % Ensure equal scaling on all axes
view(3);    % Set view to 3D

hold off;
