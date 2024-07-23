// Define the header guards

#ifndef GEOMAG_H
#define GEOMAG_H

// Include constants file

#include "../../include/constants.h"
#define NaN log(-1.0)

// Declare functions
void igrf_field(double dlat, double dlon, double altm, double time, double mag_field[3]);
static void E0000(int IENTRY, int *maxdeg, double alt,double glat,double glon, double time, double *dec, double *dip, double *ti, double *gv);
void geomag(int *maxdeg);
void geomg1(double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);
char geomag_introduction(double epochlowlim);

#endif