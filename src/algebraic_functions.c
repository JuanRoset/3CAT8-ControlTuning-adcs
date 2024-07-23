// Include source file
#include "../include/algebraic_functions.h"

// Define the functions

void inverse_3x3_Matrix(double matrix[3][3], double result[3][3]){
    // Function to compute the inverse of a 3x3 matrix

    // Calculate the determinant of the input matrix
    double det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) - 
                 matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) + 
                 matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

    // Check if the determinant is close to zero, indicating singularity and inability to compute inverse
    if (fabs(det) < 1e-10) {
        fprintf(stderr, "Matrix is singular, cannot compute inverse.\n");
        return;
    }

    // Calculate the inverse determinant for efficiency
    double invDet = 1.0 / det;

    // Calculate the elements of the inverse matrix using the adjugate and inverse determinant
    result[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * invDet;
    result[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * invDet;
    result[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * invDet;

    result[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * invDet;
    result[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * invDet;
    result[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * invDet;

    result[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * invDet;
    result[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * invDet;
    result[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * invDet;
}

void cross_product_3(double vec1[3], double vec2[3], double result[3]){
    // Function for computing the cross product of two vectors of size 3

    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

double dot_product(double vector_a[], double vector_b[], int size){
    // Function for computing the dot product of two vectors of arbitrary (specified) size

    // Initialize result to zero
    double result = 0;
    for (int i = 0; i < size; ++i) {
        result += vector_a[i] * vector_b[i];
    }

    return result;
}

double vector_norm(double vector[], int size){
    // Function to find the norm of a given vector
    // Compute the array's norm
    double norm = 0;
    double norm_squared = 0;
    for(int i = 0; i < size; i++) {
        norm_squared += pow(vector[i], 2);
    }
    norm = sqrt(norm_squared);

    return norm;
}

void normalize_vector(double vector[], int size){
    // Function to normalize a given vector
    // The input-output array should be a 1 dimensional double array

    // Compute the array's norm
    double norm = vector_norm(vector, size);

    // Divide individual elements by the array's norm
    for(int i = 0; i < size; i++) {
        vector[i] /= norm;
    }
}

double normalize_angle(double initial_angle){
    // Function to normalize a given angle (in radians) to be within a range of 0, 2 * pi

    // Initialize solution angle
    double angle = 0.0;

    // Find if angle is outside of said range
    if(initial_angle < 0 || initial_angle > 2 * pi) {
        // Adjust angle to be within range by subtracting the nº of times it's outside it
        double n_rev = floor(initial_angle / (2 * pi));

        // Return corrected angle
        return initial_angle - n_rev * 2 * pi;
    }

    // If it isn't, return the same angle
    return initial_angle;
}

double sign(double number){
    // Function to return the sign of a double

    return (number < 0) ? - 1.0 : + 1.0;
}

void transpose_3x3_Matrix(double matrix[3][3], double result[3][3]){

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = matrix[j][i];
        }
    }

}

void multiply_3x3_Matrices(double mat1[3][3], double mat2[3][3], double result[3][3]) {
    // Initialize result matrix with zeros
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = 0;
        }
    }

    // Perform matrix multiplication
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) {
                sum += mat1[i][k] * mat2[k][j];
            }
            result[i][j] = sum;
        }
    }
}

void matrix_multiply(double *mat1, double *mat2, int rows1, int cols1, int cols2, double *result) {
    // Perform matrix multiplication

    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            double sum = 0.f;
            for (int k = 0; k < cols1; ++k) {
                sum += mat1[i * cols1 + k] * mat2[k * cols2 + j];
            }
            result[i * cols2 + j] = sum;
        }
    }
}

void quaternion_product(double input_p[4], double input_q[4], double* output){
    // Function for computing the product between two quaternions
    // Input arrays must have a size of 4

    // Extract the real parts of the quaternions
    double pw = input_p[0];
    double qw = input_q[0];

    // Extract the imaginary parts of the quaternions
    double pv[3] = {input_p[1], input_p[2], input_p[3]};
    double qv[3] = {input_q[1], input_q[2], input_q[3]};

    // Compute the dot product of the imaginary parts of the input quaternions
    double dot_product_pqv = 0;
    for(int i = 0; i < 3; i++)
        dot_product_pqv += pv[i] * qv[i];

    // Compute elements of the quaternion
    output[0] = pw * qw - dot_product_pqv;
    output[1] = pw * qv[0] + qw * pv[0] + pv[1] * qv[2] - pv[2] * qv[1];
    output[2] = pw * qv[1] + qw * pv[1] + pv[2] * qv[0] - pv[0] * qv[2];
    output[3] = pw * qv[2] + qw * pv[2] + pv[0] * qv[1] - pv[1] * qv[0];
}

void rotation_quaternion(double theta, double euler_axis[], double* output){
    // Function for computing the rotation quaternion for a given rotation angle (theta) and euler axis
    // Input euler axis must have a size of 3

    // Normalize rotation axis vector
    double norm = 0;
    for(int i = 0; i < 3; i++)
        norm += euler_axis[i] * euler_axis[i];
    norm = sqrt(norm);
    for(int i = 0; i < 3; i++)
        euler_axis[i] = euler_axis[i] / norm;

    // Declare quaternion of point coordinates
    double point_quaternion[4] = {0.0, output[0], output[1], output[2]};

    // Compute the rotation quaternion
    output[0] = cos(theta / 2);
    output[1] = euler_axis[0] * sin(theta / 2);
    output[2] = euler_axis[1] * sin(theta / 2);
    output[3] = euler_axis[2] * sin(theta / 2);
}

void quaternion_rotate(double quaternion[], double* output, const char *mode){
    // Function for performing the rotation of a point according to a rotation quaternion
    // Input quaternion must be of size 4 and represent a rotation quaternion
    // Input/output point must be an array of size 3 representing coordinates in 3D space

    // Declare quaternion for point coordinates
    double point_quaternion[4] = {0.0, output[0], output[1], output[2]};

    // Compute the conjugate of the rotation quaternion
    double quaternion_conjugate[4] = {quaternion[0], - quaternion[1], - quaternion[2], - quaternion[3]};

    // Declare arrays for intermediate and final results
    double intermediate[4];
    double final[4];

    // Passive rotation rotates the point about the fixed axes
    if (strcmp(mode, "active") == 0) {
        // Perform the necessary quaternion products for active rotation
        quaternion_product(quaternion, point_quaternion, intermediate);
        quaternion_product(intermediate, quaternion_conjugate, final);
    }
    // Active rotation expresses the point with respect to the rotated axes
    else if (strcmp(mode, "passive") == 0) {
        // Perform the necessary quaternion products for passive rotation
        quaternion_product(quaternion_conjugate, point_quaternion, intermediate);
        quaternion_product(intermediate, quaternion, final);
    }
    // Display error if no mode has been selected
    else {
        printf("Error in quaternion rotation: no mode was selected \n");
        return;
    }

    // Update the point with its new coordinates
    output[0] = final[1];
    output[1] = final[2];
    output[2] = final[3];
}

void attitude_quaternion_differentiate(double angular_rate[3], double quaternion[4], double* output){
    // Function for computing the derivative of the attitude quaternion for a given attitude quaternion and body axes angular rate
    // The angular rate must be a 1D array of length 3 and the attitude quaternion a 1D array of length 4
    // The output, or time derivative of the attitude quaternion, must enter as a 1D array of length 4

    // Declare and populate angular rate matrix
    double omega_matrix[4][4] = {
        {0, -angular_rate[0], -angular_rate[1], -angular_rate[2]},
        {angular_rate[0], 0, angular_rate[2], -angular_rate[1]},
        {angular_rate[1], -angular_rate[2], 0, angular_rate[0]},
        {angular_rate[2], angular_rate[1], -angular_rate[0], 0}
    };

    // Perform matrix multiplication and multiply by 1/2 (quaternion derivative formula)
    for(int m = 0; m < 4; m++){
        // Compute result for this position
        double result = 0;
        for(int n = 0; n < 4; n++){
            result += omega_matrix[m][n] * quaternion[n];
        }
        // Assign result for this position
        output[m] = result / 2.0;
    }
}


void normalize_quaternion(double* output){
    // Function for normalizing the attitude quaternion after numerical integration
    // The output must enter as a 1D array of length 4

    // If scalar part is negative, multiply by -1
    double sign = 1.0;    
    if(output[0] < 0.0)
        sign = - sign;

    for(int i = 0; i < 4; i++)
        output[i] *= sign;
    
    // Normalize the array (of size 4 <--> quaternion)
    normalize_vector(output, 4);
    
}

void quaternion_error(double q_desired[4], double q_actual[4], double* output){
    // Function for computing the error between two quaternions
    // Input arrays must have a size of 4

    // Create copies of the input quaternions to avoid modifying the originals
    double q_desired_copy[4] = {q_desired[0], q_desired[1], q_desired[2], q_desired[3]};
    double q_actual_copy[4] = {q_actual[0], q_actual[1], q_actual[2], q_actual[3]};

    // Normalize both quaternions
    normalize_quaternion(q_desired_copy);
    normalize_quaternion(q_actual_copy);

    // Extract the inverse of the desired quaternion
    double q_desired_inverse[4] = {q_desired_copy[0], -q_desired_copy[1], -q_desired_copy[2], -q_desired_copy[3]};

    // Apply the error formula
    quaternion_product(q_desired_inverse, q_actual_copy, output);

    // Normalize the resulting quaternion error
    normalize_quaternion(output);
}


void rotation_matrix_to_quaternion(double rotation[3][3], double quaternion[4]){
    // Function to find the rotation quaternion that corresponds to a given direction cosine Matrix
    // The function uses Shepperd's method for obtaining the quaternion)
    

    // Compute the numerical indicator
    double numerical_indicator[4] = {
        rotation[0][0] + rotation[1][1] + rotation[2][2],
        rotation[0][0],
        rotation[1][1],
        rotation[2][2]
    };

    // Initialize parameters of the solution
    int i_best = 0;
    double maximum = - DBL_MAX;

    // Find the best numerically suited solution
    for(int i = 0; i < 4; i++){
        if(numerical_indicator[i] > maximum){
            maximum = numerical_indicator[i];
            i_best = i;
        }
    }

    i_best = 3;
    
    // Compute the corresponding solution
    double temp_term;
    switch (i_best) {
        case 0:
            temp_term = sqrt(1 + rotation[0][0] + rotation[1][1] + rotation[2][2]);
            quaternion[0] = temp_term / 2;
            quaternion[1] = (rotation[2][1] - rotation[1][2]) / (2 * temp_term);
            quaternion[2] = (rotation[0][2] - rotation[2][0]) / (2 * temp_term);
            quaternion[3] = (rotation[1][0] - rotation[0][1]) / (2 * temp_term);
            break;
        case 1:
            temp_term = sqrt(1 + rotation[0][0] - rotation[1][1] - rotation[2][2]);
            quaternion[0] = (rotation[2][1] - rotation[1][2]) / (2 * temp_term);
            quaternion[1] = temp_term / 2;
            quaternion[2] = (rotation[0][1] + rotation[1][0]) / (2 * temp_term);
            quaternion[3] = (rotation[2][0] + rotation[0][2]) / (2 * temp_term);
            break;
        case 2:
            temp_term = sqrt(1 - rotation[0][0] + rotation[1][1] - rotation[2][2]);
            quaternion[0] = (rotation[0][2] - rotation[2][0]) / (2 * temp_term);
            quaternion[1] = (rotation[0][1] + rotation[1][0]) / (2 * temp_term);
            quaternion[2] = temp_term / 2;
            quaternion[3] = (rotation[1][2] + rotation[2][1]) / (2 * temp_term);
            break;
        case 3:
            temp_term = sqrt(1 - rotation[0][0] - rotation[1][1] + rotation[2][2]);
            quaternion[0] = (rotation[1][0] - rotation[0][1]) / (2 * temp_term);
            quaternion[1] = (rotation[2][0] + rotation[0][2]) / (2 * temp_term);
            quaternion[2] = (rotation[2][1] + rotation[1][2]) / (2 * temp_term);
            quaternion[3] = temp_term / 2;
            break;
    }


    // Normalize the output quaternion
    normalize_quaternion(quaternion);
}