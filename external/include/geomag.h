// Define the header guards

#ifndef GEOMAG_H
#define GEOMAG_H

// Include constants file

#include "../../include/constants.h"

// Include libraries


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Define Nan number

#define NaN log(-1.0)

// Declare coefficient file path

//extern const char* wmm_path;

// Declare functions

void igrf_field(double dlat, double dlon, double altm, double time, double mag_field[3], const char *wmm_path);
static void E0000(int IENTRY, int *maxdeg, double alt,double glat,double glon, double time, double *dec, double *dip, double *ti, double *gv, const char *wmm_path);
void geomag(int *maxdeg, const char *wmm_path);
void geomg1(double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv, const char *wmm_path);
char geomag_introduction(double epochlowlim);

#endif