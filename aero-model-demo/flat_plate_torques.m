% Close all existing figures
close all
clear
clc

% Generate an array of beta values from pi to 0
beta_array = linspace(pi, 0);

% Define surface area and velocity
Surface = 1;
Velocity = 7800;

% Define velocity vector
V = [-Velocity; 0.0];

% Define constants
rho = 6.03e-14;
sigma_n = 0.8;
sigma_t = 0.8;

% Initialize arrays to store aerodynamic coefficients
CL = zeros(length(beta_array), 1);
CD = zeros(length(beta_array), 1);

% Iterate through each beta value
for i = 1:length(beta_array)
    % Extract current beta
    beta = beta_array(i);
    
    % Compute area vector based on current beta
    A = Surface .* [sin(beta); cos(beta)];
    
    % Normalize velocity and area vectors
    V_unit = V ./ Velocity;
    A_unit = A ./ Surface;
    
    % Compute aerodynamic force components
    F = rho * (V' * V) * ((2 - sigma_n - sigma_t) * (V_unit' * A_unit)^2 * A_unit + sigma_t * (V_unit' * A_unit) * V_unit) * Surface;
    Lift = F(2);
    Drag = -F(1);
    
    % Compute aerodynamic coefficients and store them
    CL(i) = 2 * Lift / (rho * Surface * Velocity^2);
    CD(i) = -2 * Drag / (rho * Surface * Velocity^2);
end

% Plot aerodynamic coefficients against angle of incidence
figure('Position', [0, 0, 700, 225])
hold on
plot(beta_array .* (180 / pi), CD, '-.k', 'LineWidth', 1.5)
plot(beta_array .* (180 / pi), CL, '--k', 'LineWidth', 1.25)

% Approximate values
%plot(beta_array.*180/pi, sin(beta_array).*max(CD))
%plot(beta_array.*180/pi, sin(beta_array.*2).*max(CL))

grid off
box off
title(sprintf('Flat plate with $\\sigma_n = %.2f$, $\\sigma_t = %.2f$', sigma_n, sigma_t), 'Interpreter', 'latex');
xlabel('Angle of incidence ($^\circ$)', 'Interpreter', 'latex')
ylabel('Aerodynamic coefficients', 'Interpreter', 'latex')
legend('$C_D$', '$C_L$', 'Location', 'best', 'Interpreter', 'latex')
legend('boxoff');

% Set X axis to be at the origin
ax = gca;
ax.XAxisLocation = 'origin';

% Set tick labels to LaTeX font
ax.TickLabelInterpreter = 'latex';

% Set custom x-axis tick marks
xticks(0:30:max(beta_array .* (180 / pi)))

% Position the x-label manually to avoid overlap
xlabelPos = get(gca, 'Position');
xlabelHeight = 0.05; % Adjust this value according to your plot
xlabel('Angle of incidence ($^\circ$)', 'Position', [200, -xlabelHeight*8, 0], 'Interpreter', 'latex');
