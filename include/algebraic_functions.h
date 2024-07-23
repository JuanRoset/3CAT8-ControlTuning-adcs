// Declare the header guards

#ifndef ALGEBRAIC_FUNCTIONS_H
#define ALGEBRAIC_FUNCTIONS_H

// Include libraries

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

// Include other source files

#include "constants.h"

#ifdef __cplusplus
extern "C" {
#endif

// Declare matrix functions

void inverse_3x3_Matrix(double matrix[3][3], double result[3][3]);
void cross_product_3(double vec1[3], double vec2[3], double result[3]);
double dot_product(double vector_a[], double vector_b[], int size);
double vector_norm(double vector[], int size);
void normalize_vector(double vector[], int size);
double normalize_angle(double initial_angle);
double sign(double number);
void transpose_3x3_Matrix(double matrix[3][3], double result[3][3]);
void multiply_3x3_Matrices(double mat1[3][3], double mat2[3][3], double result[3][3]);
void matrix_multiply(double *mat1, double *mat2, int rows1, int cols1, int cols2, double *result);

// Declare quaternion functions

void quaternion_product(double input_p[4], double input_q[4], double* output);
void rotation_quaternion(double theta, double euler_axis[], double* output);
void quaternion_rotate(double quaternion[], double* output, const char *mode);
void attitude_quaternion_differentiate(double angular_rate[3], double quaternion[4], double* output);
void normalize_quaternion(double* output);
void quaternion_error(double input_p[4], double input_q[4], double* output);
void rotation_matrix_to_quaternion(double rotation[3][3], double quaternion[4]);

#ifdef __cplusplus
}
#endif

#endif // ALGEBRAIC_FUNCTIONS_H