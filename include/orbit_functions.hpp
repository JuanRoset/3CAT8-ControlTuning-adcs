// Declare the header guards

#ifndef ORBIT_H
#define ORBIT_H

// Include other source files

#include "constants.h"
#include "algebraic_functions.h"
#include "space_environment.h"
#include "../external/include/json.hpp"
using json = nlohmann::json;

// Include libraries

#include <iostream>
#include <fstream>
#include <cmath>

// Create struct for saving the orbital information

struct orbitalElements {

    // Values obtained from json
    double epoch = 0.0;
    double start_epoch = 0.0;
    double inclination = 0.0;
    double right_ascension = 0.0;
    double eccentricity = 0.0;
    double argument_of_perigee = 0.0;
    double mean_anomaly = 0.0;
    double mean_motion = 0.0;
    double eccentric_anomaly = 1.0;

    // Values derived from previous values
    double period = 0.0;
    double orbit_start_time = 0.0;
    double mean_angular_motion = 0.0;
    double semimajor_axis = 0.0;
    double angular_momentum = 0.0;
    double rotation_matrix[3][3] = {0.0};
    double normal_vector[3] = {0.0};

    // Start condition time variables
    int start_year = 0;
    int start_doy = 0;
    double start_day_fraction = 0.0;

};

struct sunActivity {

    // Solar activity parameters
    double f107A = 0.0; // f10.7 solar flux of 81 day average
    double f107d = 0.0; // f10.7 solar flux of daily reading
    double ap = 0.0; // ap magnetic index

};

// Define functions

double find_eccentric_anomaly(double eccentricity, double mean_anomaly, double initial_guess, double tolerance, int max_iterations);
void orbit_propagate(orbitalElements elements, double time, double* position, double* velocity, double* quaternion_rsw, double long_lat[2], double magnetic_field[3], double sun_position[3], double& eclipse, double& geomagnetic_inclination, double& density, std::string wmm_path);
void orbit_fetch(orbitalElements& elements, std::string filePath);
void sun_fetch(sunActivity& indices, std::string filePath);
double epoch_to_JD(int year,int day_in_year, double day_fraction);
double JD_to_GMST(double julian_date);

#endif // ORBIT_H