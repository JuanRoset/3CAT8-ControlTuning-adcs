// Declare the header guards

#ifndef ACTUATOR_FUNCTIONS_H
#define ACTUATOR_FUNCTIONS_H

// Include other source files

#include "algebraic_functions.h"
#include "../external/include/json.hpp"
using json = nlohmann::json;

// Create struct for saving the magnetorquers' parameters

struct magnetorquersParameters{
    // Physical parameters
    double surface[3] = {2.0108736, 2.0108736, 2.0108736};
    double resistance[3] = {39.42735, 39.42735, 39.42735};

    // Interface parameters
    double maximum_individual_current = 0.15;
    double maximum_total_current = 0.20;
    double maximum_total_power = 4.5;
    double maximum_duty_cycle = 0.5;
};

struct reactionWheelParameters{
    // Physical parameters
    double maximum_angular_rate = 8000 * 2 * pi / 60; // In 8000 rpm in radians per second
    double maximum_angular_acceleration = 113.35; // In radians per second square
    double mainInertia = 0.00000203;

    // Power draw parameters
    double current_kappa = 0.0035; // Current model constant for acceleration [A / Nm]
    double current_coefficients[6] = {0.009468, 1.994e-5, -0.0006313, 2.987e-8, -1.06e-6, 3.047e-5}; // Reaction wheel current model for maintaining certain angular rate
    double u_supply = 12.0; // Voltage of the power supply to the reaction wheel

    // Propagation parameters
    double propagation_numerator_coefficients[4]   = {0.194197817509306,   0.566275183225970,   0.564462042192420,   0.192384676475756};
    double propagation_denominator_coefficients[3] = {1.349253162033861,  -0.241033473812238,  -0.590899968818171};

    // State memory
    double past_references[4]    = {0.0, 0.0, 0.0, 0.0};
    double past_angular_rates[3] = {0.0, 0.0, 0.0};
};

// Define functions

void magnetorquer_model(magnetorquersParameters magnetorquers,
                        double magnetic_field_body[3],
                        double desired_magnetic_moment[3],
                        double resulting_torque[3],
                        double resulting_current[3],
                        double resulting_power[3],
                        double resulting_duty_cycle[3]);

void reactionWheel_model(reactionWheelParameters& reactionWheel,
                         double deltaTime,
                         double& target_torque,
                         double& resulting_torque,
                         double& resulting_current,
                         double& resulting_power);

#endif // ACTUATOR_FUNCTIONS_H