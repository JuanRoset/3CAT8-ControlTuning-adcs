#ifndef SIMULATION_UTILS_H
#define SIMULATION_UTILS_H

// Include other source files

#include "constants.h"
#include "orbit_functions.hpp"
#include "../external/include/json.hpp"
using json = nlohmann::json;


// Include libraries

#include <iostream>
#include <sstream>
#include <fstream>

// Declare the used functions

void satellite_fetch(double inertia[3][3], double & inertia_principal_min, std::string & control_mode, double Kp[3], double Kd[3], double & RW_relaxation, double & Ku, bool & active_force_control, bool & RW_unloading, int& configuration, std::string filePath);
void simulation_fetch(double & simulation_orbits, double & delta_time, bool & propagate_orbit, bool & gravitational_disturbance, bool & aerodynamic_disturbance, bool & solarPressure_disturbance, double euler_axis[3], double & theta, std::string filePath);
void write_orbit_state(double delta_time, orbitalElements elements, sunActivity indices, std::string orbitPath, std::string wmm_path, const double iteration_steps);
void read_orbit_state(double delta_time,
                      double simulation_time,
                      std::string orbitPath,
                      std::vector<std::vector<double>>& velocity_array,
                      std::vector<std::vector<double>>& position_array,
                      std::vector<std::vector<double>>& quaternion_rsw_array,
                      std::vector<double>& time_array,
                      std::vector<std::vector<double>>& longitude_latitude,
                      std::vector<std::vector<double>>& magnetic_field_array,
                      std::vector<std::vector<double>>& sun_position,
                      std::vector<double>& eclipse,
                      std::vector<double>& geomagnetic_inclination,
                      std::vector<double>& air_density);


#endif // SIMULATION_UTILS_H