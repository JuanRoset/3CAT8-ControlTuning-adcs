#ifndef DISTURBANCE_TORQUES_H
#define DISTURBANCE_TORQUES_H


// Include libraries

#include <cstdio>
#include <cstdlib>
#include <iostream>

// Include source files

#include "surface_torques_model.hpp"
#include "algebraic_functions.h"

// Declare functions

void gravitational_gradient_torque(double inertia[3][3], double attitude_quaternion[4], double position[3], double position_norm, double *torque_body);
void aerodynamic_torque(int config, double air_velocity[3], double density, double * torque_body);
void srp_torque(int config, double ray[3], double pressure, double * torque_body);

// Functions based on the regression model

void aerodynamic_torque_regr(int config, double air_velocity[3], double density, double * torque_body);
void srp_torque_regr(int config, double ray[3], double pressure, double * torque_body);

#endif