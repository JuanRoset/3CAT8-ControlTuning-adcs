% Define the vertices of the triangle
triangle = [1, 1; 4, 5; 7, 2];

% Generate grid of points
[X, Y] = meshgrid(1:0.1:8, 1:0.1:6);
points = [X(:), Y(:)];

% Initialize array to store results
inside = zeros(size(points, 1), 1);

% Check each point if it's inside the triangle
for i = 1:size(points, 1)
    inside(i) = point_in_triangle(triangle, points(i, :));
end

% Plotting
figure;
hold on;

% Plot the triangle
plot([triangle(:, 1); triangle(1, 1)], [triangle(:, 2); triangle(1, 2)], 'b-', 'LineWidth', 2);

% Plot points inside the triangle in green
plot(points(inside == 1, 1), points(inside == 1, 2), 'g.', 'MarkerSize', 1);

% Plot points outside the triangle in red
plot(points(inside == 0, 1), points(inside == 0, 2), 'r.', 'MarkerSize', 1);

title('Points inside and outside the triangle');
xlabel('X');
ylabel('Y');
legend('Triangle', 'Points inside', 'Points outside');
axis equal;
grid on;
hold off;
