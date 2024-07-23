// Include libraries

#include <iostream>
#include <fstream>
#include <math.h>

// Include the header files with the necessary functions

#include "../include/actuator_functions.hpp"

// Define functions

void magnetorquer_model(magnetorquersParameters magnetorquers, double magnetic_field_body[3], double desired_magnetic_moment[3], double resulting_torque[3], double resulting_current[3], double resulting_power[3], double resulting_duty_cycle[3]){
    // Function that emulates the behaviour of real magnetorquers for producing the final torque, with limited current draws

    // Find the unlimited values of the currents through each coil
    double mtq_current[3] = {desired_magnetic_moment[0] / magnetorquers.surface[0],
                             desired_magnetic_moment[1] / magnetorquers.surface[1],
                             desired_magnetic_moment[2] / magnetorquers.surface[2]};
    //printf("Desired current: %.4f, %.4f, %.4f (A) \n", mtq_current[0], mtq_current[1], mtq_current[2]);

    // Find the maximum unlimited current value
    int max_current_index = 0;
    double max_current = abs(mtq_current[0]);
    for(int i = 1; i < 3; i++){
        if(abs(mtq_current[i]) > max_current){
            max_current = abs(mtq_current[i]);
            max_current_index = i;
        }
    }

    // Limit the individual currents to the maximum possible individual value
    double k_individual = 1.0;
    if(max_current > magnetorquers.maximum_individual_current * magnetorquers.maximum_duty_cycle){
        k_individual = magnetorquers.maximum_individual_current * magnetorquers.maximum_duty_cycle / max_current;
    }

    resulting_current[0] = mtq_current[0] * k_individual;
    resulting_current[1] = mtq_current[1] * k_individual;
    resulting_current[2] = mtq_current[2] * k_individual;

    // Limit the currents to the maximum of the total current consumed
    double k_total = 1.0;
    double total_current = abs(resulting_current[0]) + abs(resulting_current[1]) + abs(resulting_current[2]);
    if(total_current > magnetorquers.maximum_total_current * magnetorquers.maximum_duty_cycle){
        k_total = magnetorquers.maximum_total_current * magnetorquers.maximum_duty_cycle / total_current;
    }
    //printf("Total current: %.4f (A),   ", total_current);
    //printf("Maximum total current: %.4f (A)\n", magnetorquers.maximum_total_current * magnetorquers.maximum_duty_cycle);
    //printf("K total: %.4f \n", k_total);

    resulting_current[0] = resulting_current[0] * k_total;
    resulting_current[1] = resulting_current[1] * k_total;
    resulting_current[2] = resulting_current[2] * k_total;

    //printf("Resulting current: %.4f, %.4f, %.4f (A) \n \n", resulting_current[0], resulting_current[1], resulting_current[2]);

    // Limit the total power draw
    double k_power = 1.0;
    double total_power = 0.0;
    for(int i = 0; i < 3; i++)
        total_power += magnetorquers.resistance[i] * pow(resulting_current[i], 2);
    if(total_power > magnetorquers.maximum_total_power)
        k_power = sqrt(magnetorquers.maximum_total_power / total_power);

    resulting_current[0] = resulting_current[0] * k_power;
    resulting_current[1] = resulting_current[1] * k_power;
    resulting_current[2] = resulting_current[2] * k_power;


    // Compute the resulting power drawn by the magnetorquers
    resulting_power[0] = magnetorquers.resistance[0] * pow(resulting_current[0], 2);
    resulting_power[1] = magnetorquers.resistance[1] * pow(resulting_current[1], 2);
    resulting_power[2] = magnetorquers.resistance[2] * pow(resulting_current[2], 2);

    // Compute the duty cycle being used on each magnetorquer
    resulting_duty_cycle[0] = resulting_current[0] * magnetorquers.maximum_individual_current;
    resulting_duty_cycle[1] = resulting_current[1] * magnetorquers.maximum_individual_current;
    resulting_duty_cycle[2] = resulting_current[2] * magnetorquers.maximum_individual_current;

    // Find the limited magnetic moment
    double limited_magnetic_moment[3] = {resulting_current[0] * magnetorquers.surface[0],
                                         resulting_current[1] * magnetorquers.surface[1],
                                         resulting_current[2] * magnetorquers.surface[2]};

    // Compute the resulting torque exerted by the magnetorquers
    cross_product_3(limited_magnetic_moment, magnetic_field_body, resulting_torque);

    //printf("Resulting torque: %.12f, %.12f, %.12f \n \n", resulting_torque[0], resulting_torque[1], resulting_torque[2]);
}

void reactionWheel_model(reactionWheelParameters& reactionWheel, double deltaTime, double& target_torque, double& resulting_torque, double& resulting_current, double& resulting_power){
    // Function to numerically propagate the angular rate of the reaction wheel using a difference equation of its discrete TF
    // This function also computes the current and power drawn, and updates the past references and rates
    // This function needs the 3 past angular rates and reference values for computing the next ones

    // Compute the desired angular acceleration for a given torque
    double target_angular_acceleration = - target_torque / reactionWheel.mainInertia;

    // Compute the next desired angular rate and update the reference values
    double target_next_rate = reactionWheel.past_angular_rates[0] + target_angular_acceleration * deltaTime;
    reactionWheel.past_references[3] = reactionWheel.past_references[2];
    reactionWheel.past_references[2] = reactionWheel.past_references[1];
    reactionWheel.past_references[1] = reactionWheel.past_references[0];
    reactionWheel.past_references[0] = target_next_rate;

    // Compute the next theoretical angular rate
    double next_rate = dot_product(reactionWheel.propagation_numerator_coefficients, reactionWheel.past_references, 4) - dot_product(reactionWheel.propagation_denominator_coefficients, reactionWheel.past_angular_rates, 3);

    // Limit the output to the maximum possible angular rate
    if(abs(next_rate) > reactionWheel.maximum_angular_rate)
        next_rate = sign(next_rate) * reactionWheel.maximum_angular_rate;

    // Limit the output to that corresponding to the maximum possible angular acceleration
    double angular_acceleration = (next_rate - reactionWheel.past_angular_rates[0]) / deltaTime;
    if(abs(angular_acceleration) > reactionWheel.maximum_angular_acceleration){
        angular_acceleration = sign(angular_acceleration) * reactionWheel.maximum_angular_acceleration;
        next_rate = reactionWheel.past_angular_rates[0] + angular_acceleration * deltaTime;
    }

    // printf("Target rate: %.8f ,   Next rate: %.8f \n", target_next_rate, next_rate);

    // Update the angular rate array
    reactionWheel.past_angular_rates[2] = reactionWheel.past_angular_rates[1];
    reactionWheel.past_angular_rates[1] = reactionWheel.past_angular_rates[0];
    reactionWheel.past_angular_rates[0] = next_rate;
    
    // Compute the resulting torque
    resulting_torque = - (reactionWheel.past_angular_rates[0] - reactionWheel.past_angular_rates[1]) / deltaTime * reactionWheel.mainInertia;

    // Compute the current being drawn for a given angular rate and acceleration at a given supply voltage
    double current_model_array[6] = {1.0, abs(reactionWheel.past_angular_rates[0]), reactionWheel.u_supply, pow(reactionWheel.past_angular_rates[0], 2), abs(reactionWheel.past_angular_rates[0]) * reactionWheel.u_supply, pow(reactionWheel.u_supply, 2)};
    resulting_current = abs(angular_acceleration) * reactionWheel.mainInertia / reactionWheel.current_kappa + dot_product(reactionWheel.current_coefficients, current_model_array, 6);
    
    // Compute the corresponding current being drawn
    resulting_power = reactionWheel.u_supply * resulting_current;
}