// Declare the header guards

#ifndef SPACE_ENVIRONMENT_H
#define SPACE_ENVIRONMENT_H

// Include other source files

#include "constants.h"
#include "algebraic_functions.h"
#include "../external/include/geomag.h"
#include "../external/include/nrlmsise-00.h"

// Include libraries

#include <math.h>

// Define functions

#ifdef __cplusplus
extern "C" {
#endif

void dipole_magnetic_field(double position[3], double radius, double earth_rotation, double field[3]);
void igrf_magnetic_field(double long_lat[2], double position[3], double radius, double date, double field[3]);
double geomagneticInclination(double orbit_normal[3], double earth_rotation);
double calculate_density(int year, int doy, int seconds, double altitude, double g_lat, double g_long, double lst, double f107A, double f107, double ap);

#ifdef __cplusplus
}
#endif

#endif // SPACE_ENVIRONMENT_H