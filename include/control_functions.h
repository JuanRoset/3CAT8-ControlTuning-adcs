// Declare the header guards

#ifndef CONTROL_H
#define CONTROL_H

// Include other source files

#include "algebraic_functions.h"

// Include libraries

#include <math.h>

// Define functions

#ifdef __cplusplus
extern "C" {
#endif

void pd_control_torque(double error_quaternion[4], double Kp[3], double angular_rate[3], double Kd[3], double torque_body[3]);
void pid_control_torque(double error_quaternion[4], double Kp[3], double int_quaternion[4], double Ki[3], double angular_rate[3], double Kd[3], double torque_body[3]);
void magnetorquer_torque_to_moment(double magnetic_field_body[3], double desired_torque[3], double resulting_moment[3]);
void magnetorquer_detumbling_moment(double magnetic_field_body[3], double angular_rate[3], double kw, double resulting_moment[3]);
void estimate_disturbance_torque(const double inertia[3][3], double angular_acceleration[3], double angular_rate[3], double control_torque[3], double *output);

#ifdef __cplusplus
}
#endif

#endif // CONTROL_H