// Include the header files with the necessary functions

#include "../include/control_functions.h"

// Define functions

void pd_control_torque(double error_quaternion[4], double Kp[3], double angular_rate[3], double Kd[3], double torque_body[3]){
    // Function for computing the control torques on body axes
    // The input and output torque should be a 1D array of length 3

    // Disturbance-free torque
    for(int i = 0; i < 3; i++)
        torque_body[i] = - (Kp[i] * error_quaternion[i + 1] + Kd[i] * angular_rate[i]);
        
}

void pid_control_torque(double error_quaternion[4], double Kp[3], double int_quaternion[4], double Ki[3], double angular_rate[3], double Kd[3], double torque_body[3]){
    // Disturbance-free torque
    for(int i = 0; i < 3; i++){
        //torque_body[i] = - (sign(error_quaternion[0]) * (Kp[i] + Ki[i] * int_quaternion[0]) * error_quaternion[i + 1] + Kd[i] * angular_rate[i]);
        torque_body[i] = - (sign(error_quaternion[0]) * Kp[i] * error_quaternion[i + 1] + sign(int_quaternion[0]) * Ki[i] * int_quaternion[i + 1] + Kd[i] * angular_rate[i]);
    }
}

void magnetorquer_torque_to_moment(double magnetic_field_body[3], double desired_torque[3], double resulting_moment[3]){

    // Find the resulting magnetic dipole moment
    cross_product_3(magnetic_field_body, desired_torque, resulting_moment);

    // Normalize torque figures
    double norm_B_squared = pow(vector_norm(magnetic_field_body, 3), 2);
    for(int i = 0; i < 3; i++)
        resulting_moment[i] /= norm_B_squared;
}

void magnetorquer_detumbling_moment(double magnetic_field_body[3], double angular_rate[3], double kw, double resulting_moment[3]){
    // Function for computing the magnetic moment to be applied for satellite detumbling
    // The input and output torque should be a 1D array of length 3

    // Find the resulting magnetic dipole moment
    cross_product_3(magnetic_field_body, angular_rate, resulting_moment);

    // Normalize torque figures
    double norm_B_squared = pow(vector_norm(magnetic_field_body, 3), 2);
    for(int i = 0; i < 3; i++)
        resulting_moment[i] *= - kw / norm_B_squared;
    

}

void estimate_disturbance_torque(const double inertia[3][3], double angular_acceleration[3], double angular_rate[3], double control_torque[3], double *output){

    double term_Ia[3];
    for(int r = 0; r < 3; r++){
        double sum = 0.0;
        for(int c = 0; c < 3; c++){
            sum += inertia[r][c] * angular_acceleration[c];
        }
        term_Ia[r] = sum;
    }

    double term_Iw[3];
    for(int r = 0; r < 3; r++){
        double sum = 0.0;
        for(int c = 0; c < 3; c++){
            sum += inertia[r][c] * angular_rate[c];
        }
        term_Iw[r] = sum;
    }

    double term_cross[3];
    cross_product_3(angular_rate, term_Iw, term_cross);
    
    for(int r = 0; r < 3; r++){
        output[r] = term_Ia[r] + term_cross[r] - control_torque[r];
    }
}