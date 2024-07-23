#ifndef NRLMSISE_00_H
#define NRLMSISE_00_H

// Include libraries

#include <math.h>          /* maths functions */
#include <stdio.h>         /* for error messages. TBD: remove this */
#include <stdlib.h>        /* for malloc/free */

// Define structs

struct nrlmsise_flags {
	int switches[24];
	double sw[24];
	double swc[24];
};

struct ap_array {
	double a[7];   
};

struct nrlmsise_input {
	int year;      /* year, currently ignored */
	int doy;       /* day of year */
	double sec;    /* seconds in day (UT) */
	double alt;    /* altitude in kilometers */
	double g_lat;  /* geodetic latitude */
	double g_long; /* geodetic longitude */
	double lst;    /* local apparent solar time (hours), see note below */
	double f107A;  /* 81 day average of F10.7 flux (centered on doy) */
	double f107;   /* daily F10.7 flux for previous day */
	double ap;     /* magnetic index(daily) */
	struct ap_array *ap_a; /* see above */
};

struct nrlmsise_output {
	double d[9];   /* densities */
	double t[2];   /* temperatures */
};


// Declare functions

void gtd7 (struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output);
void gtd7d(struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output);
void gts7 (struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output);
void ghp7 (struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output, double press);

#ifdef INLINE
#define __inline_double static inline double
#else
#define __inline_double double
#endif

#endif /* NRLMSISE_00_H */