/*        _   _   _ _             _         _____            _             _    _____ _                 _       _                      
     /\  | | | | (_) |           | |       / ____|          | |           | |  / ____(_)               | |     | |                     
    /  \ | |_| |_ _| |_ _   _  __| | ___  | |     ___  _ __ | |_ _ __ ___ | | | (___  _ _ __ ___  _   _| | __ _| |_ ___  _ __          
   / /\ \| __| __| | __| | | |/ _` |/ _ \ | |    / _ \| '_ \| __| '__/ _ \| |  \___ \| | '_ ` _ \| | | | |/ _` | __/ _ \| '__|         
  / ____ \ |_| |_| | |_| |_| | (_| |  __/ | |___| (_) | | | | |_| | | (_) | |  ____) | | | | | | | |_| | | (_| | || (_) | |            
 /_/    \_\__|\__|_|\__|\__,_|\__,_|\___|  \_____\___/|_| |_|\__|_|  \___/|_| |_____/|_|_| |_| |_|\__,_|_|\__,_|\__\___/|_|      

  ____            _   ___             _    _ _____   _____   _   _                   _____       _     _           _                   
 |___ \          | | / _ \           | |  | |  __ \ / ____| | \ | |                 / ____|     | |   | |         | |                  
   __) | ___ __ _| || (_) |  ______  | |  | | |__) | |      |  \| | __ _ _ __   ___| (___   __ _| |_  | |     __ _| |__                
  |__ < / __/ _` | __> _ <  |______| | |  | |  ___/| |      | . ` |/ _` | '_ \ / _ \\___ \ / _` | __| | |    / _` | '_ \               
  ___) | (_| (_| | || (_) |          | |__| | |    | |____  | |\  | (_| | | | | (_) |___) | (_| | |_  | |___| (_| | |_) |              
 |____/ \___\__,_|\__\___/            \____/|_|     \_____| |_| \_|\__,_|_| |_|\___/_____/ \__,_|\__| |______\____|___*/

// ===================================================================== //
// This is the main simulation code for the ADCS of the 3cat-8 satellite //
// Author: Pau Climent Salazar - pau.climent@estudiantat.upc.edu         //
// ===================================================================== //


// Include libraries

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <vector>

// Include the source files with the necessary functions

#include "../include/actuator_functions.hpp"
#include "../include/simulation_utils.hpp"
#include "../include/disturbance_torques.hpp"
#include "../include/control_functions.h"
#include "../include/algebraic_functions.h"

// Specify config file paths

std::string filePathOutput            = "../data/simulation_output.csv";
std::string orbitPath                 = "../data/orbit_file.csv";
std::string satelliteParamFilePath    = "../parameters/satellite_parameters.json";
std::string astronomicalParamFilePath = "../parameters/astronomical_parameters.json";

// Declare satellite parameters

int configuration; // 0 = pq, 1 = nopq, 2 = antenna
double inertia[3][3];
double inertia_principal_min;
double inertia_wheels[3][3] = {0.0};

// ꩜ Declare Detumbling control parameters ꩜

bool detumbled = false;
double detumbling_threshold = 2 * pi / 180;

// ⊗ Declare Pointing control parameters ⊗

double Kp[3];
double Kd[3];
bool active_force_control;
bool RW_unloading;
double RW_relaxation; // RW relaxation constant - set to 1
double Ku; // unloading constant

// ⁂ Set the desired disturbances ⁂

bool gravitational_disturbance;
bool aerodynamic_disturbance;
bool solarPressure_disturbance;

// Declare initial offset parameters

double euler_axis[3];
double theta;


// Define the nadir rotation matrix
double rotation_matrix_nadir[3][3] = {
    {- 1.0,   0.0,   0.0},
    {  0.0,   0.0, - 1.0},
    {  0.0, - 1.0,   0.0}
};
double rotation_matrix_nadir_inverse[3][3];

// Define simulation parameters
double simulation_orbits;
double delta_time;
bool propagate_orbit;

std::string integration_method = "Euler";

// Declare variables
double angular_rate_body[3]   = {0.0};
double attitude_quaternion[4] = {0.0};
double quaternion_o2t[4]      = {0.0};

// Declare functions
void euler_equation_body_wheels(double torque_body[3], double angular_rate[3], double wheels_angular_rate[3], double inertia_sat[3][3], double inertia_wheels[3][3], double inverse_inertia[3][3], double angular_acceleration[3]);
void numerical_propagation(double state[], double derivative[], int length, const double delta_time, double* next_state, std::string method);


// Main int

int main(){

    // ✓ Declare control mode ✓
    std::string control_mode; // Mission = 0 -> no contr, 1 -> ideal, 2 -> detumbling, 3 -> mtq, 4 -> hybrid

    // Print status
    std::cout << "Retrieving orbital infromation" << std::endl;
    simulation_fetch(simulation_orbits, delta_time, propagate_orbit, gravitational_disturbance, aerodynamic_disturbance, solarPressure_disturbance, euler_axis, theta, satelliteParamFilePath);

    // Read orbital parameters from TLE file
    orbitalElements elements;
    orbit_fetch(elements, satelliteParamFilePath);

    // Read sun activity parameters from astronomical parameters file
    sunActivity sun_indices;
    sun_fetch(sun_indices, astronomicalParamFilePath);

    // Initialize the actuators' parameters
    magnetorquersParameters magnetorquers;
    reactionWheelParameters reactionWheel;

    //printf("Orbital period is %.12f seconds\n", elements.period);
    double simulation_time = simulation_orbits * elements.period;
    int iteration_steps = int(simulation_time / delta_time);

     // Preallocate and store propagation of orbit positions for simulation time
    std::vector<std::vector<double>> velocity_array(3, std::vector<double>(iteration_steps));
    std::vector<std::vector<double>> position_array(3, std::vector<double>(iteration_steps));
    std::vector<std::vector<double>> quaternion_rsw_array(4, std::vector<double>(iteration_steps));
    std::vector<double> time_array(iteration_steps);   
    std::vector<std::vector<double>> longitude_latitude(2, std::vector<double>(iteration_steps));
    std::vector<std::vector<double>> magnetic_field_array(3, std::vector<double>(iteration_steps));
    std::vector<std::vector<double>> sun_position(3, std::vector<double>(iteration_steps));
    std::vector<double> eclipse(iteration_steps);
    std::vector<double> geomagnetic_inclination(iteration_steps);
    std::vector<double> air_density(iteration_steps);

    // If stated, propagate and store orbit state
    if(propagate_orbit == true){
        write_orbit_state(delta_time, elements, sun_indices, orbitPath, iteration_steps);
    }

    // Print status
    std::cout << "Retrieving satellite parameters" << std::endl;

    // Populate the inertia matrix of the reaction wheel
    inertia_wheels[3][3] = reactionWheel.mainInertia;

    // Fetch all the satellita and sim parameters
    satellite_fetch(inertia, inertia_principal_min, control_mode, Kp, Kd, RW_relaxation, Ku, active_force_control, RW_unloading, configuration, satelliteParamFilePath);

    // Print status
    std::cout << "Reading orbital state" << std::endl;

    // Store necessary data from the orbit variables
    read_orbit_state(delta_time, simulation_time, orbitPath, velocity_array, position_array, quaternion_rsw_array, time_array, longitude_latitude, magnetic_field_array, sun_position, eclipse, geomagnetic_inclination, air_density);
    
    // Compute mean angular motion for the detumbling
    double omega_o = elements.mean_motion * 2 * pi / (24.0 * 3600.0);


    // Find inverse of inertia matrix
    double total_inertia[3][3] = {0.0};
    for(int i = 0; i < 3; i++){
        double result = 0;
        for(int j = 0; j < 3; j++)
            total_inertia[i][j] = inertia[i][j] + inertia_wheels[i][j];
    }

    double total_inverse_inertia[3][3] = {0};
    inverse_3x3_Matrix(total_inertia, total_inverse_inertia);

    // Declare the initial angular rate in the body reference frame
    if (control_mode == "detumble"){
        detumbled = false;
        angular_rate_body[0] = 0.506145 / sqrt(3.0) / 1.0;
        angular_rate_body[1] = 0.506145 / sqrt(3.0) / 1.0;
        angular_rate_body[2] = 0.506145 / sqrt(3.0) / 1.0;
    }
    else{
        angular_rate_body[0] = 0.0000;
        angular_rate_body[1] = - omega_o;
        angular_rate_body[2] = 0.0000;
    }

    // Open the file for writing
    std::ofstream outFile(filePathOutput);

    // Check if the file is opened successfully
    if (!outFile.is_open()) {
        std::cerr << "Error opening the file: " << filePathOutput << std::endl;
        return 1; // Exit with an error code
    }

    // Write output data header
    outFile << "Time,AngularRateX,AngularRateY,AngularRateZ,Quaternion_e2b_W,Quaternion_e2b_X,Quaternion_e2b_Y,Quaternion_e2b_Z,Quaternion_e2o_W,Quaternion_e2o_X,Quaternion_e2o_Y,Quaternion_e2o_Z,Quaternion_error_W,Quaternion_error_X,Quaternion_error_Y,Quaternion_error_Z,Current_mtq_X,Current_mtq_Y,Current_mtq_Z,Power_mtq_X,Power_mtq_Y,Power_mtq_Z,Torque_mtq_X,Torque_mtq_Y,Torque_mtq_Z,Torque_RW,AngRate_RW,Current_RW,Power_RW,Gg_Torque_X,Gg_Torque_Y,Gg_Torque_Z,Air_Torque_X,Air_Torque_Y,Air_Torque_Z,SRP_Torque_X,SRP_Torque_Y,SRP_Torque_Z\n";

    // Preallocate memory for state variables
    double torque_body[3]                    = {0.0};
    double angular_acceleration_body[3]      = {0.0};
    double angular_acceleration_fixed[3]     = {0.0};
    double angular_rate_f2b[3]               = {0.0};
    double angular_rate_f2o[3]               = {0.0};
    double angular_rate_o2b[3]               = {0.0};
    double angular_rate_o2b_body[3]          = {0.0};
    double attitude_quaternion_derivative[4] = {0.0};
    double disturbance_torque_body[3]        = {0.0};
    double control_torque_body[3]            = {0.0};
    double magnetic_field_body[3]            = {0.0};
    double quaternion_b2f[4]                 = {0.0};
    double target_torque_body[3]             = {0.0};
    double current_draw[3]                   = {0.0};
    double power_draw[3]                     = {0.0};
    // Variables for the pending updates
    double quaternion_f2o[4]                 = {0.0};
    double quaternion_f2t[4]                 = {0.0};
    double quaternion_t2b[4]                 = {0.0};
    double quaternion_t2b_int[4]             = {0.0};
    double position_f[3]                     = {0.0};
    double velocity_f[3]                     = {0.0};
    // Variables for actuators
    double target_magnetic_moment[3]         = {0.0};
    double magnetorquer_torque_body[3]       = {0.0};
    double magnetorquer_current[3]           = {0.0};
    double magnetorquer_power[3]             = {0.0};
    double magnetorquer_duty_cycle [3]       = {0.0};
    double reactionWheel_rate                = 0.0;
    double reactionWheel_target              = 0.0;
    double reactionWheel_acceleration        = 0.0;
    double reactionWheel_current             = 0.0;
    double reactionWheel_power               = 0.0;
    double reactionWheel_torque              = 0.0;
    double average_kw                        = 0.0;
    double gravity_disturbance_torque[3]     = {0.0};
    double aerodynamic_disturbance_torque[3] = {0.0};
    double srp_disturbance_torque[3]         = {0.0};
    double disturbance_torques_estimate[3]   = {0.0};

    // Set the past Reaction Wheel angular rates and reference rates as constants with the given initial values
    for(int i = 0; i < 3; i++){
        reactionWheel.past_angular_rates[i]  = reactionWheel_rate;
        reactionWheel.past_references[i + 1] = reactionWheel_rate;
    }
    reactionWheel.past_references[0] = reactionWheel_rate;

    inverse_3x3_Matrix(rotation_matrix_nadir, rotation_matrix_nadir_inverse);
    rotation_matrix_to_quaternion(rotation_matrix_nadir_inverse, quaternion_o2t);

    // Update Orbit frame of reference quaternion
    quaternion_f2o[0] = quaternion_rsw_array[0][0];
    quaternion_f2o[1] = quaternion_rsw_array[1][0];
    quaternion_f2o[2] = quaternion_rsw_array[2][0];
    quaternion_f2o[3] = quaternion_rsw_array[3][0];

    // Compute the initial quaternion -- Now in NADIR
    quaternion_product(quaternion_f2o, quaternion_o2t, attitude_quaternion);

    // Define parameters for rotation
    double perturbance_quaternion[4];
    rotation_quaternion(theta, euler_axis, perturbance_quaternion);

    // Compute the initial attitude quaternion with an initial perturbance
    quaternion_product(attitude_quaternion, perturbance_quaternion, attitude_quaternion);

    // Print status
    std::cout << "Simulating attitude" << std::endl;

    // Loop through each time step
    for(int t = 0; t < iteration_steps; t++){

        // Update Orbit frame of reference quaternion
        quaternion_f2o[0] = quaternion_rsw_array[0][t];
        quaternion_f2o[1] = quaternion_rsw_array[1][t];
        quaternion_f2o[2] = quaternion_rsw_array[2][t];
        quaternion_f2o[3] = quaternion_rsw_array[3][t];

        // Compute the target quaternion -- Now in NADIR
        quaternion_product(quaternion_f2o, quaternion_o2t, quaternion_f2t);

        // Compute the quaternion error (inputs must be the attitude quaternion, and the target quaternion wrt. to the same frame - ECI)
        quaternion_error(quaternion_f2t, attitude_quaternion, quaternion_t2b);

        quaternion_t2b_int[0] += delta_time * quaternion_t2b[0];
        quaternion_t2b_int[1] += delta_time * quaternion_t2b[1];
        quaternion_t2b_int[2] += delta_time * quaternion_t2b[2];
        quaternion_t2b_int[3] += delta_time * quaternion_t2b[3];

        // Calculate the angular rate in the fixed frame
        quaternion_b2f[0] = + attitude_quaternion[0];
        quaternion_b2f[1] = - attitude_quaternion[1];
        quaternion_b2f[2] = - attitude_quaternion[2];
        quaternion_b2f[3] = - attitude_quaternion[3];

        angular_rate_f2b[0] = angular_rate_body[0];
        angular_rate_f2b[1] = angular_rate_body[1];
        angular_rate_f2b[2] = angular_rate_body[2];

        quaternion_rotate(quaternion_b2f, angular_rate_f2b, "passive");

        // Calculate the norm of the position vector
        position_f[0] = position_array[0][t];
        position_f[1] = position_array[1][t];
        position_f[2] = position_array[2][t];
        double position_norm = sqrt(pow(position_f[0], 2) + pow(position_f[1], 2) + pow(position_f[2], 2));

        // Calculate the norm of the velocity vector
        velocity_f[0] = velocity_array[0][t];
        velocity_f[1] = velocity_array[1][t];
        velocity_f[2] = velocity_array[2][t];
        double velocity_norm = sqrt(pow(velocity_f[0], 2) + pow(velocity_f[1], 2) + pow(velocity_f[2], 2));

        // Calculate the angular rate of the orbit frame
        cross_product_3(position_f, velocity_f, angular_rate_f2o);
        for(int i = 0; i < 3; i++) angular_rate_f2o[i] /= pow(position_norm, 2);

        // Calculate the relative angular rate from the orbit frame to the body frame in the fixed frame
        angular_rate_o2b[0] = angular_rate_f2b[0] - angular_rate_f2o[0];
        angular_rate_o2b[1] = angular_rate_f2b[1] - angular_rate_f2o[1];
        angular_rate_o2b[2] = angular_rate_f2b[2] - angular_rate_f2o[2];

        // Rotate the relative angular rate of the orbit from the fixed to the body frame
        angular_rate_o2b_body[0] = angular_rate_o2b[0];
        angular_rate_o2b_body[1] = angular_rate_o2b[1];
        angular_rate_o2b_body[2] = angular_rate_o2b[2];

        quaternion_rotate(attitude_quaternion, angular_rate_o2b_body, "passive");

        // Compute the magnetic field in body axes
        magnetic_field_body[0] = magnetic_field_array[0][t];
        magnetic_field_body[1] = magnetic_field_array[1][t];
        magnetic_field_body[2] = magnetic_field_array[2][t];

        quaternion_rotate(attitude_quaternion, magnetic_field_body, "passive");

        // Calculate the Gravitational disturbance torque in body axes
        gravitational_gradient_torque(total_inertia, attitude_quaternion, position_f, position_norm, gravity_disturbance_torque);

        // Calculate the Relative velocity of the oncoming air
        double atm_angular_rate[3] = {0.0, 0.0, omega_Earth};
        double atm_velocity[3] = {0.0, 0.0, 0.0};
        double air_velocity[3] = {0.0, 0.0, 0.0};

        cross_product_3(atm_angular_rate, position_f, atm_velocity);
        for(int i = 0; i < 3; i++)
            air_velocity[i] = atm_velocity[i] - velocity_f[i];

        double air_velocity_body[3] = {air_velocity[0], air_velocity[1], air_velocity[2]};
        quaternion_rotate(attitude_quaternion, air_velocity_body, "passive");

        // Calculate the Aerodynamic disturbance torque in body axes
        aerodynamic_torque(configuration, air_velocity_body, air_density[t], aerodynamic_disturbance_torque);

        // Calculate the SRP disturbance torque in body axes
        double light_ray[3] = {-sun_position[0][t], -sun_position[1][t], -sun_position[2][t]};
        quaternion_rotate(attitude_quaternion, light_ray, "passive");

        srp_torque(configuration, light_ray, sun_pressure * eclipse[t], srp_disturbance_torque);

        // Set the disturbance torques
        if (!gravitational_disturbance) for(int g=0; g < 3; g++)     gravity_disturbance_torque[g] = 0.0;
        if (!aerodynamic_disturbance)   for(int a=0; a < 3; a++) aerodynamic_disturbance_torque[a] = 0.0;
        if (!solarPressure_disturbance) for(int s=0; s < 3; s++)         srp_disturbance_torque[s] = 0.0;
        
        // Add up all disturbance torques
        for (int d = 0; d < 3; d++) disturbance_torque_body[d] = gravity_disturbance_torque[d] + aerodynamic_disturbance_torque[d] + srp_disturbance_torque[d];

        // Estimate disturbance torques from previous iteration for active disturbance rejection
        estimate_disturbance_torque(inertia, angular_acceleration_body, angular_rate_body, control_torque_body, disturbance_torques_estimate);

        // Perform attitude control according to the selected control mode
        if(control_mode == "no_torque"){
            // Assign null control torques
            control_torque_body[0] = 0.0;
            control_torque_body[1] = 0.0;
            control_torque_body[2] = 0.0;
        }

        if(control_mode == "ideal"){
            // Calculate the ideal control torque in body axes
            pd_control_torque(quaternion_t2b, Kp, angular_rate_o2b_body, Kd, target_torque_body);

            // Assign the ideal torques
            control_torque_body[0] = target_torque_body[0];
            control_torque_body[1] = target_torque_body[1];
            control_torque_body[2] = target_torque_body[2];
        }

        if(control_mode == "detumble"){
            // Eliminate satellite's initial angular rate

            // Compute modulus of angular rate
            double angular_rate = vector_norm(angular_rate_body, 3);

            // Compute modified bdot control gain
            double kw = 2 * omega_o * (1.0 + sin(geomagnetic_inclination[t])) * inertia_principal_min;
            average_kw += kw;

            // Compute target magnetic moment for detumbling
            magnetorquer_detumbling_moment(magnetic_field_body, angular_rate_body, kw, target_magnetic_moment);

             // Compute the actuator torque by modelling the reaction wheel
            magnetorquer_model(magnetorquers, magnetic_field_body, target_magnetic_moment, magnetorquer_torque_body, magnetorquer_current, magnetorquer_power, magnetorquer_duty_cycle);

            // Assign control torques to those created by the magnetorquers and the reaction wheel
            control_torque_body[0] = magnetorquer_torque_body[0];
            control_torque_body[1] = magnetorquer_torque_body[1];
            control_torque_body[2] = magnetorquer_torque_body[2];

            // Handle state logic
            if (angular_rate <= detumbling_threshold){
                if (detumbled == false){
                    printf("Detumbled in %.4f seconds\n", t * delta_time);
                    average_kw /= t;
                    printf("Kw constant %.12f\n", average_kw);
                }
                detumbled = true;
                control_mode = "mtq_afc";
            }
        }

        if(control_mode == "mtq_only"){
            // Calculate the actuation magnetic moment needed using the magnetic control law

            // Calculate the ideal control torque in body axes
            pd_control_torque(quaternion_t2b, Kp, angular_rate_o2b_body, Kd, target_torque_body);

            // Remove estimate disturbance torque
            if(active_force_control){
                for(int r = 0; r < 3; r++){
                    target_torque_body[r] -= disturbance_torques_estimate[r];
                }
            }

            // Compute the magnetic moment needed by the magnetorquers
            magnetorquer_torque_to_moment(magnetic_field_body, target_torque_body, target_magnetic_moment);            

            // Compute the actuator torque by modelling the reaction wheel
            magnetorquer_model(magnetorquers, magnetic_field_body, target_magnetic_moment, magnetorquer_torque_body, magnetorquer_current, magnetorquer_power, magnetorquer_duty_cycle);
            
            // Assign control torques to only those created by the magnetorquers
            control_torque_body[0] = magnetorquer_torque_body[0];
            control_torque_body[1] = magnetorquer_torque_body[1];
            control_torque_body[2] = magnetorquer_torque_body[2];
        }

        if(control_mode == "hybrid"){
            // Calculate the ideal control torque in body axes
            pd_control_torque(quaternion_t2b, Kp, angular_rate_o2b_body, Kd, target_torque_body);

            // Remove estimate disturbance torque
            if(active_force_control){
                for(int r = 0; r < 3; r++){
                    target_torque_body[r] -= disturbance_torques_estimate[r];
                }
            }

            // Compute the magnetic moment needed by the magnetorquers
            magnetorquer_torque_to_moment(magnetic_field_body, target_torque_body, target_magnetic_moment);

            // Perform momentum unloading
            double moment_unloading[3] = {0.0, 0.0, 0.0};
            if(RW_unloading){
                // Compute error between target and current RW angular rate
                double reactionWheel_error[3] = {0.0, 0.0, reactionWheel_target - reactionWheel_rate};

                // Compute magnetic moment according to proportional law
                cross_product_3(magnetic_field_body, reactionWheel_error, moment_unloading);
                for(int rw = 0; rw < 3; rw++){
                    moment_unloading[rw] *= Ku / pow(vector_norm(magnetic_field_body, 3), 2);

                    // Add unloading magnetic moment to target total moment
                    target_magnetic_moment[rw] += moment_unloading[rw];
                }
            }


            // Compute the actuator torque by modelling the reaction wheel
            magnetorquer_model(magnetorquers, magnetic_field_body, target_magnetic_moment, magnetorquer_torque_body, magnetorquer_current, magnetorquer_power, magnetorquer_duty_cycle);

            // Compute the torque that must be supplied by the reaction wheel
            double target_reaction_wheel_torque = (target_torque_body[2] - magnetorquer_torque_body[2]) * RW_relaxation;

            // Simulate the model of the reaction wheel
            reactionWheel_model(reactionWheel, delta_time, target_reaction_wheel_torque, reactionWheel_torque, reactionWheel_current, reactionWheel_power);
            reactionWheel_rate = reactionWheel.past_angular_rates[0];
            reactionWheel_acceleration = reactionWheel_torque / reactionWheel.mainInertia;

            // Assign control torques to those created by the magnetorquers and the reaction wheel
            control_torque_body[0] = magnetorquer_torque_body[0];
            control_torque_body[1] = magnetorquer_torque_body[1];
            control_torque_body[2] = magnetorquer_torque_body[2] + reactionWheel_torque;

        }

        // Compute the total torque in the body axes
        torque_body[0] = control_torque_body[0] + disturbance_torque_body[0];
        torque_body[1] = control_torque_body[1] + disturbance_torque_body[1];
        torque_body[2] = control_torque_body[2] + disturbance_torque_body[2];

        // Calculate the body angular acceleration using Euler's equation
        double angular_rate_wheels[3] = {0.0, 0.0, reactionWheel_rate};
        euler_equation_body_wheels(torque_body, angular_rate_body, angular_rate_wheels, inertia, inertia_wheels, total_inverse_inertia, angular_acceleration_body);
        //euler_equation_body(torque_body, angular_rate_body, inertia, total_inverse_inertia, angular_acceleration_body);

        // Propagate the fixed frame Angular rate with its acceleration and derivatives respectively
        numerical_propagation(angular_rate_body, angular_acceleration_body, 3, delta_time, angular_rate_body, integration_method);

        // Calculate the attitude quaternion derivative
        attitude_quaternion_differentiate(angular_rate_body, attitude_quaternion, attitude_quaternion_derivative);

        // Propagate the Attitude quaternion with its derivative respectively
        numerical_propagation(attitude_quaternion, attitude_quaternion_derivative, 4, delta_time, attitude_quaternion, integration_method);

        // Normalize the  next quaternion
        normalize_quaternion(attitude_quaternion);

        // Write the current time and data to the output datafile
        outFile << delta_time * t << ","
                << angular_rate_body[0] << "," << angular_rate_body[1] << "," << angular_rate_body[2] << ","
                //<< angular_rate_o2b_body[0] << "," << angular_rate_o2b_body[1] << "," << angular_rate_o2b_body[2] << ","
                << attitude_quaternion[0] << "," << attitude_quaternion[1] << "," << attitude_quaternion[2] << "," << attitude_quaternion[3] << ","
                << quaternion_f2o[0] << "," << quaternion_f2o[1] << "," << quaternion_f2o[2] << "," << quaternion_f2o[3] << ","
                << quaternion_t2b[0] << "," << quaternion_t2b[1] << "," << quaternion_t2b[2] << "," << quaternion_t2b[3] << "," 
                << magnetorquer_current[0] << "," << magnetorquer_current[1] << "," << magnetorquer_current[2] << ","
                << magnetorquer_power[0] << "," << magnetorquer_power[1] << "," << magnetorquer_power[2] << ","
                << magnetorquer_torque_body[0] << "," << magnetorquer_torque_body[1] << "," << magnetorquer_torque_body[2]  << ","
                << reactionWheel_torque << ","
                << reactionWheel_rate << ","
                << reactionWheel_current << ","
                << reactionWheel_power << ","
                << gravity_disturbance_torque[0] << "," << gravity_disturbance_torque[1] << "," << gravity_disturbance_torque[2] << ","
                << aerodynamic_disturbance_torque[0] << "," << aerodynamic_disturbance_torque[1] << "," << aerodynamic_disturbance_torque[2] << ","
                << srp_disturbance_torque[0] << "," << srp_disturbance_torque[1] << "," << srp_disturbance_torque[2] << "\n";


    }

    // Call external visualization scripts
    system("python ../visualization/control_visualizer.py");
    system("python ../visualization/3D_visualizer.py");
    system("python ../visualization/orbit_visualizer.py");

    return 0;
}

// Define the main-native functions

void euler_equation_body_wheels(double torque_body[3], double angular_rate[3], double wheels_angular_rate[3], double inertia_sat[3][3], double inertia_wheels[3][3], double inverse_inertia[3][3], double angular_acceleration[3]){
    // Function for computing the angular acceleration of the body on body axes for a given torque, angular rate and inertia tensor
    // The output angular acceleration should be a 1D array of length 3

    // TODO: Add reaction wheel components

    // Compute the matrix product between the total inertia tensor and angular rate + the wheels' inertia and their angular rate
    double IcWb[3] = {0};

    for(int m = 0; m < 3; m++){
        double result = 0;
        for(int i = 0; i < 3; i++)
            result += (inertia_sat[m][i] + inertia_wheels[m][i]) * angular_rate[i] + inertia_wheels[m][i] * wheels_angular_rate[i];
        
        IcWb[m] = result;
    }

    // Compute the cross product of the angular rate times the IcWb array
    double cross_product[3] = {0};
    cross_product_3(angular_rate, IcWb, cross_product);


    // Subtract torque by cross product
    double torque_minus_cross[3] = {0};
    torque_minus_cross[0] = torque_body[0] - cross_product[0];
    torque_minus_cross[1] = torque_body[1] - cross_product[1];
    torque_minus_cross[2] = torque_body[2] - cross_product[2];

    // Multiply inverse of inertia by 
    for(int m = 0; m < 3; m++){
        double result = 0;
        for(int i = 0; i < 3; i++)
            result += inverse_inertia[m][i] * torque_minus_cross[i];
        
        angular_acceleration[m] = result;
    }

}

void numerical_propagation(double state[], double derivative[], int length, const double delta_time, double* next_state, std::string method){
    // Function to numerically integrate a given state vector forwards in time

    // Apply simple Euler method
    if(method == "Euler"){
        for(int i = 0; i < length; i++){
            next_state[i] = state[i] + derivative[i] * delta_time;
        }
    }

    // Other methods have been removed but shall be re-implemented in the future for further acccuracy
}