function [force_body, torque_body] = srp_torque(sat_shape, position_CG, light_ray, pressure)

load(sat_shape)

% Initialize torque
torque_body = [0, 0, 0];

% Specify air ray
light_ray = light_ray ./ norm(light_ray);

for i = 1:max(size(surfaces))

    S_normal = surface_normals(i, :);
    ray_dot = dot(S_normal, light_ray);

    if sign(ray_dot) <0 
    
        S_value = surface_values(i);
        rho_s = surface_rho_s(i);
        rho_d = surface_rho_d(i);

        lever_arm = surface_centers(i, :) - position_CG;
        
        % Compute aerodynamic force and torque
        F = - pressure * S_value * ray_dot .* ((1 - rho_s) * light_ray + 2 * (rho_s * ray_dot - rho_d / 3) * S_normal);
        T = cross(lever_arm, F);

        % Add contribution to the total torque
        force_body = force_body + F;
        torque_body = torque_body + T;

    end

end

end

% void srp_torque(int config, double ray[3], double pressure, double * torque_body){
%     // Function to compute the aerodynamic torques on the satellite
% 
%     // Find the ray unit vector of air
%     normalize_vector(ray, 3);
% 
%     // Compute pitch and yaw angles
%     double pitch, yaw;
%     ray_to_angles(ray, pitch, yaw);
% 
%     // Select the right coefficients for the regression model
%     double * coefs_x, * coefs_y, * coefs_z;
%     int length_x, length_y, length_z;
%     switch (config) {
%     case 0: // Case of Solar Panel + PocketCubes
% 
%         coefs_x = (double *) &coefficients_theta_x_srp_pq[0];
%         coefs_y = (double *) &coefficients_theta_y_srp_pq[0];
%         coefs_z = (double *) &coefficients_theta_z_srp_pq[0];
% 
%         length_x = coefficients_theta_x_srp_pq_length;
%         length_y = coefficients_theta_y_srp_pq_length;
%         length_z = coefficients_theta_z_srp_pq_length;
% 
%         break;
% 
%     case 1: // Case of Solar Panel without PocketCubes
% 
%         coefs_x = (double *) &coefficients_theta_x_srp_nopq[0];
%         coefs_y = (double *) &coefficients_theta_y_srp_nopq[0];
%         coefs_z = (double *) &coefficients_theta_z_srp_nopq[0];
% 
%         length_x = coefficients_theta_x_srp_nopq_length;
%         length_y = coefficients_theta_y_srp_nopq_length;
%         length_z = coefficients_theta_z_srp_nopq_length;
% 
%         break;
% 
%     case 2: // Case of Solar Panel + Antenna
% 
%         coefs_x = (double *) &coefficients_theta_x_srp_antenna[0];
%         coefs_y = (double *) &coefficients_theta_y_srp_antenna[0];
%         coefs_z = (double *) &coefficients_theta_z_srp_antenna[0];
% 
%         length_x = coefficients_theta_x_srp_antenna_length;
%         length_y = coefficients_theta_y_srp_antenna_length;
%         length_z = coefficients_theta_z_srp_antenna_length;
% 
%         break;
% 
%     default:
%         break;
%     }
% 
%     // Apply the regression model
%     torque_body[0] = pressure * Xaxis_torqueC(pitch, yaw, config, coefs_x, length_x);
%     torque_body[1] = pressure * Yaxis_torqueC(pitch, yaw, config, coefs_y, length_y);
%     torque_body[2] = pressure * Zaxis_torqueC(pitch, yaw, config, coefs_z, length_z);
% 
% }