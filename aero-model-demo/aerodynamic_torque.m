function [force_body, torque_body] = aerodynamic_torque(config, air_velocity, density)

load('3cat8_geometry.mat')

% Initialize torque
force_body = [0, 0, 0];
torque_body = [0, 0, 0];

% Define CG
position_CG = [0.0, 0.0, 0.0];

% Specify air ray
Airspeed_val = norm(air_velocity);
Airspeed_dir = air_velocity ./ Airspeed_val;

for i = 1:max(size(surfaces))

    S_normal = surface_normals(i, :);
    ray_dot = dot(S_normal, Airspeed_dir);

    if sign(ray_dot) <0 
    
        S_value = surface_values(i);
        sigma_t = surface_sigma_t(i);
        sigma_n = surface_sigma_n(i);

        lever_arm = surface_centers(i, :) - position_CG;
        
        % Compute aerodynamic force and torque
        F = - density * Airspeed_val^2 * S_value * ((2 - sigma_n - sigma_t) * ray_dot^2 * S_normal + sigma_t * ray_dot * Airspeed_dir);
        T = cross(lever_arm, F);

        % Add contribution to the total torque
        force_body = force_body + F;
        torque_body = torque_body + T;

    end

end

end

% void aerodynamic_torque(int config, double velocity[3], double density, double * torque_body){
%     // Function to compute the aerodynamic torques on the satellite
% 
%     // Find the ray unit vector of air
%     double ray[3] = {velocity[0], velocity[1], velocity[2]};
%     normalize_vector(ray, 3);
% 
%     // Find norm of velocity
%     double velocity_norm = vector_norm(velocity, 3);
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
%         coefs_x = (double *) &coefficients_theta_x_air_pq[0];
%         coefs_y = (double *) &coefficients_theta_y_air_pq[0];
%         coefs_z = (double *) &coefficients_theta_z_air_pq[0];
% 
%         length_x = coefficients_theta_x_air_pq_length;
%         length_y = coefficients_theta_y_air_pq_length;
%         length_z = coefficients_theta_z_air_pq_length;
% 
%         break;
% 
%     case 1: // Case of Solar Panel without PocketCubes
% 
%         coefs_x = (double *) &coefficients_theta_x_air_nopq[0];
%         coefs_y = (double *) &coefficients_theta_y_air_nopq[0];
%         coefs_z = (double *) &coefficients_theta_z_air_nopq[0];
% 
%         length_x = coefficients_theta_x_air_nopq_length;
%         length_y = coefficients_theta_y_air_nopq_length;
%         length_z = coefficients_theta_z_air_nopq_length;
% 
%         break;
% 
%     case 2: // Case of Solar Panel + Antenna
% 
%         coefs_x = (double *) &coefficients_theta_x_air_antenna[0];
%         coefs_y = (double *) &coefficients_theta_y_air_antenna[0];
%         coefs_z = (double *) &coefficients_theta_z_air_antenna[0];
% 
%         length_x = coefficients_theta_x_air_antenna_length;
%         length_y = coefficients_theta_y_air_antenna_length;
%         length_z = coefficients_theta_z_air_antenna_length;
% 
%         break;
% 
%     default:
%         break;
%     }
% 
%     // Apply the regression model
%     double test_torque[3];
%     test_torque[0] = density * pow(velocity_norm, 2) * Xaxis_torqueC(pitch, yaw, config, coefs_x, length_x);
%     test_torque[1] = density * pow(velocity_norm, 2) * Yaxis_torqueC(pitch, yaw, config, coefs_y, length_y);
%     test_torque[2] = density * pow(velocity_norm, 2) * Zaxis_torqueC(pitch, yaw, config, coefs_z, length_z);
% 
%     // Assign results if the obtained torque is not Nan
%     for(int i = 0; i < 3; i++){
%         if(test_torque[i] == test_torque[i]) torque_body[i] = test_torque[i];
%         else{
%             //std::cout<<"\nAero torque returns a Nan\n";
%             //printf("Velocity: %.12f, %.12f, %.12f\n", velocity[0], velocity[1], velocity[2]);
%             //printf("rho %.12f, vel %.12f, Coefs: %.12f, %.12f, %.12f", density, velocity_norm, Xaxis_torqueC(pitch, yaw, config, coefs_x, length_x), Yaxis_torqueC(pitch, yaw, config, coefs_y, length_y), Zaxis_torqueC(pitch, yaw, config, coefs_z, length_z));
%         }
%     }
% }