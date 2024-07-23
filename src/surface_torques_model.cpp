// Include source files

#include "../include/surface_torques_model.hpp"

// Define functions

void ray_to_angles(double ray[3], double & pitch, double & yaw){
    // Apply own convention to convert ray unit vector to pitch and yaw angles

    // Compute projection of ray to yz plane
    double zy_projection[3] = {0.0, ray[1], ray[2]};
    double zy_norm = vector_norm(zy_projection, 3);
    double z_vector[3] = {0.0, 0.0, - 1.0};

    // Compute cosines of pitch and yaw
    double cos_pitch = dot_product(ray, zy_projection, 3) / zy_norm;
    double cos_yaw = - dot_product(z_vector, zy_projection, 3) / zy_norm;

    // Compute pitch and yaw angles
    pitch = - sign(ray[0]) * acos(cos_pitch);
    yaw = - sign(ray[1]) * acos(cos_yaw);
}

double Xaxis_torqueC(double pitch, double yaw, int config, const double coefs[], const int length){
    // Function for building the X-matrix for the regression model of X-axis torques

    double result = 0.0;
    double * X_array = (double *) malloc(length * sizeof(double));

    X_array[0] = 1.0;

    int order_i, order_j;
    if(config == 2){
        order_i = 3;
        order_j = 5;

        double yaw_proc = sign(sin(yaw)) * abs(sin(2 * yaw));
        X_array[1] = yaw_proc;
        
        int idxc = 2;
        for(int i = 0; i < order_i; i++){
            X_array[idxc] = pow(yaw, i + 1);
            idxc ++;
            for(int j = 0; j < order_j; j++){
                X_array[idxc] = pow(yaw_proc, i + 1) * pow(pitch, j + 1);
                idxc ++;
            }
        }
    }
    else if(config == 0 || config == 1){
        order_i = 13;
        order_j = 5;

        double yaw_temp = yaw;
        if(yaw < 0)
            yaw = - pi - yaw;
        
        int idxf, idx1, idx2;
        for(int i = 0; i < order_i; i++){
            for(int j = 0; j < order_j; j++){
                idxf = i * order_j + j;
                idx1 = 1 + 2 * idxf;
                idx2 = idx1 + 1;
                X_array[idx1] = abs(sin(yaw)) * pow(yaw, i + 1) * pow(cos(pitch), j + 1);
                X_array[idx2] = abs(sin(yaw)) * cos(pitch) * pow(yaw, i + 1) * pow(pitch, j + 1);
            }
        }
    }

    for(int k = 0; k < length; k++){
        result += coefs[k] * X_array[k];
    }

    free(X_array);

    return result;

}
double Yaxis_torqueC(double pitch, double yaw, int config, const double coefs[], const int length){
    // Function for building the X-matrix for the regression model of Y-axis torques

    double result = 0.0;
    double * X_array = (double *) malloc(length * sizeof(double));

    X_array[0] = 1.0;

    int order_i, order_j;
    if(config == 2){
        order_i = 7;
        order_j = 5;

        double yaw_proc = abs(cos(yaw));
        double pitch_proc = sin(2 * pitch);

        X_array[1] = pitch_proc;
        X_array[2] = cos(yaw / 4);

        int idxc = 2;
        for(int i = 0; i < order_i; i++){

            X_array[idxc + 1] = pow(yaw, i + 1) * yaw_proc * cos(pitch);
            idxc ++;
            
            for(int j = 0; j < order_j; j++){
                X_array[idxc + 1] = pow(yaw_proc, i + 1) * sign(pitch) * pow(abs(pitch), j + 1);
                idxc ++;
            }
        }
    }
    else if(config == 0 || config == 1){
        order_i = 7;
        order_j = 5;

        double yaw_proc = sin(2 * yaw) * sign(sin(yaw));
        double pitch_proc = cos(pitch);
        double pitch_slope = sin(pitch);

        X_array[1] = pow(cos(yaw / 2), 100) * pitch_slope;

        int idxc = 2;
        for(int i = 0; i < order_i; i++){

            if ((i + 1) % 2 != 0){
                X_array[idxc] = pow(pitch_slope, i + 1);
                idxc ++;
            }

            X_array[idxc] = pow(yaw, i + 1);
            idxc ++;
            
            for(int j = 0; j < order_j; j++){
                X_array[idxc] = pow(yaw, i + 1) * pow(pitch, j + 1) * cos(yaw / 2) * cos(pitch);
                X_array[idxc + 1] = pow(yaw_proc, i + 1) * pow(pitch_proc, j + 1);
                idxc += 2;
            }
        }
    }

    for(int k = 0; k < length; k++){
        result += coefs[k] * X_array[k];
    }

    free(X_array);

    return result;

}

double Zaxis_torqueC(double pitch, double yaw, int config, const double coefs[], const int length){
    // Function for building the X-matrix for the regression model of Z-axis torques

    double result = 0.0;
    double * X_array = (double *) malloc(length * sizeof(double));

    X_array[0] = 1.0;

    int order_i, order_j;
    if(config == 2){
        order_i = 7;
        order_j = 7;
    }
    else if(config == 0 || config == 1){
        order_i = 11;
        order_j = 7;
    }

    double yaw_proc = sin(yaw);
    double pitch_proc = cos(pitch);

    int idxc = 1;
    for (int i = 0; i < order_i; i++){
        if((i + 1) % 2 != 0){
            for(int j = 0; j < order_j; j++){
                X_array[idxc] = pow(yaw, i + 1) * pow(pitch, j + 1) * yaw_proc * pitch_proc;
                X_array[idxc + 1] = pow(yaw_proc, i + 1) * pow(pitch_proc, j + 1);
                idxc += 2;
            }
        }
    }

    for(int k = 0; k < length; k++){
        result += coefs[k] * X_array[k];
    }

    free(X_array);

    return result;

}