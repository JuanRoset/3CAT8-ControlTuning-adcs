#include "../include/disturbance_torques.hpp"

void gravitational_gradient_torque(double inertia[3][3], double attitude_quaternion[4], double position[3], double position_norm, double *torque_body){
    // Function for computing the gravitational disturbance torque on body axes
    // The input and output torque should be a 1D array of length 3

    // Normalize the position vector
    double unit_position[3];
    for(int i = 0; i < 3; i++)
        unit_position[i] = position[i] / position_norm;

    // Rotate the unit position vector to be expressed with respect to the body frame
    quaternion_rotate(attitude_quaternion, unit_position, "passive");

    // Compute the product between the inertia matrix and the rotated unit position vector
    double inertia_unit_product[3];
    for(int i = 0; i < 3; i++){
        double result = 0.0;
        for(int m = 0; m < 3; m++)
            result += inertia[i][m] * unit_position[m];
        inertia_unit_product[i] = result;
    }

    // Compute the cross product between the position unit vector and I
    double pos_inertia_pos[3];
    cross_product_3(unit_position, inertia_unit_product, pos_inertia_pos);

    // Disturbance-free torque
    torque_body[0] = 3.0 * pos_inertia_pos[0] * mu_Earth / pow(position_norm, 3);
    torque_body[1] = 3.0 * pos_inertia_pos[1] * mu_Earth / pow(position_norm, 3);
    torque_body[2] = 3.0 * pos_inertia_pos[2] * mu_Earth / pow(position_norm, 3);
}


void aerodynamic_torque(int config, double velocity[3], double density, double * torque_body){
    // Function to compute the aerodynamic torques on the satellite

    // Find the ray unit vector of air
    double ray[3] = {velocity[0], velocity[1], velocity[2]};
    normalize_vector(ray, 3);

    // Find norm of velocity
    double velocity_norm = vector_norm(velocity, 3);

    // Compute pitch and yaw angles
    double pitch, yaw;
    ray_to_angles(ray, pitch, yaw);

    // Select the right coefficients for the regression model
    double * coefs_x, * coefs_y, * coefs_z;
    int length_x, length_y, length_z;
    switch (config) {
    case 0: // Case of Solar Panel + PocketCubes

        coefs_x = (double *) &coefficients_theta_x_air_pq[0];
        coefs_y = (double *) &coefficients_theta_y_air_pq[0];
        coefs_z = (double *) &coefficients_theta_z_air_pq[0];
        
        length_x = coefficients_theta_x_air_pq_length;
        length_y = coefficients_theta_y_air_pq_length;
        length_z = coefficients_theta_z_air_pq_length;

        break;

    case 1: // Case of Solar Panel without PocketCubes
    
        coefs_x = (double *) &coefficients_theta_x_air_nopq[0];
        coefs_y = (double *) &coefficients_theta_y_air_nopq[0];
        coefs_z = (double *) &coefficients_theta_z_air_nopq[0];
        
        length_x = coefficients_theta_x_air_nopq_length;
        length_y = coefficients_theta_y_air_nopq_length;
        length_z = coefficients_theta_z_air_nopq_length;

        break;

    case 2: // Case of Solar Panel + Antenna

        coefs_x = (double *) &coefficients_theta_x_air_antenna[0];
        coefs_y = (double *) &coefficients_theta_y_air_antenna[0];
        coefs_z = (double *) &coefficients_theta_z_air_antenna[0];
        
        length_x = coefficients_theta_x_air_antenna_length;
        length_y = coefficients_theta_y_air_antenna_length;
        length_z = coefficients_theta_z_air_antenna_length;

        break;
    
    default:
        break;
    }

    // Apply the regression model
    double test_torque[3];
    test_torque[0] = density * pow(velocity_norm, 2) * Xaxis_torqueC(pitch, yaw, config, coefs_x, length_x);
    test_torque[1] = density * pow(velocity_norm, 2) * Yaxis_torqueC(pitch, yaw, config, coefs_y, length_y);
    test_torque[2] = density * pow(velocity_norm, 2) * Zaxis_torqueC(pitch, yaw, config, coefs_z, length_z);

    // Assign results if the obtained torque is not Nan
    for(int i = 0; i < 3; i++){
        if(test_torque[i] == test_torque[i]) torque_body[i] = test_torque[i];
        else{
            //std::cout<<"\nAero torque returns a Nan\n";
            //printf("Velocity: %.12f, %.12f, %.12f\n", velocity[0], velocity[1], velocity[2]);
            //printf("rho %.12f, vel %.12f, Coefs: %.12f, %.12f, %.12f", density, velocity_norm, Xaxis_torqueC(pitch, yaw, config, coefs_x, length_x), Yaxis_torqueC(pitch, yaw, config, coefs_y, length_y), Zaxis_torqueC(pitch, yaw, config, coefs_z, length_z));
        }
    }
}

void srp_torque(int config, double ray[3], double pressure, double * torque_body){
    // Function to compute the aerodynamic torques on the satellite

    // Find the ray unit vector of air
    normalize_vector(ray, 3);

    // Compute pitch and yaw angles
    double pitch, yaw;
    ray_to_angles(ray, pitch, yaw);

    // Select the right coefficients for the regression model
    double * coefs_x, * coefs_y, * coefs_z;
    int length_x, length_y, length_z;
    switch (config) {
    case 0: // Case of Solar Panel + PocketCubes

        coefs_x = (double *) &coefficients_theta_x_srp_pq[0];
        coefs_y = (double *) &coefficients_theta_y_srp_pq[0];
        coefs_z = (double *) &coefficients_theta_z_srp_pq[0];
        
        length_x = coefficients_theta_x_srp_pq_length;
        length_y = coefficients_theta_y_srp_pq_length;
        length_z = coefficients_theta_z_srp_pq_length;

        break;

    case 1: // Case of Solar Panel without PocketCubes
    
        coefs_x = (double *) &coefficients_theta_x_srp_nopq[0];
        coefs_y = (double *) &coefficients_theta_y_srp_nopq[0];
        coefs_z = (double *) &coefficients_theta_z_srp_nopq[0];
        
        length_x = coefficients_theta_x_srp_nopq_length;
        length_y = coefficients_theta_y_srp_nopq_length;
        length_z = coefficients_theta_z_srp_nopq_length;

        break;

    case 2: // Case of Solar Panel + Antenna

        coefs_x = (double *) &coefficients_theta_x_srp_antenna[0];
        coefs_y = (double *) &coefficients_theta_y_srp_antenna[0];
        coefs_z = (double *) &coefficients_theta_z_srp_antenna[0];
        
        length_x = coefficients_theta_x_srp_antenna_length;
        length_y = coefficients_theta_y_srp_antenna_length;
        length_z = coefficients_theta_z_srp_antenna_length;

        break;
    
    default:
        break;
    }

    // Apply the regression model
    torque_body[0] = pressure * Xaxis_torqueC(pitch, yaw, config, coefs_x, length_x);
    torque_body[1] = pressure * Yaxis_torqueC(pitch, yaw, config, coefs_y, length_y);
    torque_body[2] = pressure * Zaxis_torqueC(pitch, yaw, config, coefs_z, length_z);

}