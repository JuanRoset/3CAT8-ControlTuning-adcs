// Include the header files with the necessary functions

#include "../include/space_environment.h"

// Define functions

void dipole_magnetic_field(double position[3], double radius, double earth_rotation, double field[3]){
    // Function to approximate Earth's magnetic field in LEO using the simplified dipole model
    // The position vector should be normalized to 1

    // Compute geomagnetic north pole unit vector (approximation as spherical Earth)
    double north_pole[3] = {0.0, 0.0, 1};

    double latitude_axis[3] = {0.0, 1.0, 0.0};
    double latitude_quaternion[4];
    rotation_quaternion((90 - geomagnetic_latitude) * pi / 180, latitude_axis, latitude_quaternion);
    quaternion_rotate(latitude_quaternion, north_pole, "active");
    
    double longitude_axis[3] = {0.0, 0.0, 1.0};
    double longitude_quaternion[4];
    rotation_quaternion((geomagnetic_longitude + earth_rotation) * pi / 180, longitude_axis, longitude_quaternion);
    quaternion_rotate(longitude_quaternion, north_pole, "active");

    // Compute position unit vector
    double normalized_position[3];
    for(int i = 0; i < 3; i++)
        normalized_position[i] = position[i] / radius;

    //printf("Geomagnetic North Pole: %.8f, %.8f, %.8f \n", north_pole[0], north_pole[1], north_pole[2]);
    //printf("Normalized Position: %.8f, %.8f, %.8f \n", normalized_position[0], normalized_position[1], normalized_position[2]);

    // Compute colaltitude as angle between the two unit vectors
    double colatitude = acos(dot_product(normalized_position, north_pole, 3));
    //printf("Latitude: %f \n", 180 / pi * colatitude);
    
    // Compute radial and normal magnetic field
    double B_radial = - 2 * magnetic_field_equator * pow(radius_Earth / radius, 3) * cos(colatitude);
    double B_normal = magnetic_field_equator * pow(radius_Earth / radius, 3) * sin(colatitude);

    // Compute the vector in the x-y plane, parpendicular to the position vector
    double perpendicular_vector[3] = {- normalized_position[1], normalized_position[0], 0};
    normalize_vector(perpendicular_vector, 3);

    // Find the normal or tangential unit vector (normal to sphere, points north)
    double normal_vector[3] = {0};
    cross_product_3(normalized_position, perpendicular_vector, normal_vector);

    // Calculate the magnetic field in the ECI frame
    field[0] = B_radial * normalized_position[0] + B_normal * normal_vector[0];
    field[1] = B_radial * normalized_position[1] + B_normal * normal_vector[1];
    field[2] = B_radial * normalized_position[2] + B_normal * normal_vector[2];
}

void igrf_magnetic_field(double long_lat[2], double position[3], double radius, double date, double field[3], const char *wmm_path){
    // Function to approximate Earth's magnetic field in LEO using the simplified dipole model
    // The position vector should be normalized to 1

    // Function to approximate Earth's magnetic field in LEO using the simplified dipole model
    // The position vector should be normalized to 1

    // Compute geomagnetic north pole unit vector (approximation as spherical Earth)
    double north_pole[3] = {0.0, 0.0, 1};

    // Compute position unit vector
    double normalized_position[3];
    for(int i = 0; i < 3; i++)
        normalized_position[i] = position[i] / radius;

    // Compute colaltitude as angle between the two unit vectors
    double colatitude = acos(dot_product(normalized_position, north_pole, 3));
    
    // Compute radial and normal magnetic field
    double field_NED[3];
    igrf_field(long_lat[1], long_lat[0], radius, date, field_NED, wmm_path);
    //double B_radial = - 2 * magnetic_field_equator * pow(radius_Earth / radius, 3) * cos(colatitude);
    //double B_normal = magnetic_field_equator * pow(radius_Earth / radius, 3) * sin(colatitude);

    // Compute the vector in the x-y plane, parpendicular to the position vector
    double east_vector[3] = {- normalized_position[1], normalized_position[0], 0};
    normalize_vector(east_vector, 3);

    // Find the normal or tangential unit vector (normal to sphere, points north)
    double north_vector[3] = {0};
    cross_product_3(normalized_position, east_vector, north_vector);

    // Find the downlooking vector
    double down_vector[3] = {-normalized_position[0], -normalized_position[1], -normalized_position[2]};

    // Calculate the magnetic field in the ECI frame
    for(int i = 0; i < 3; i++){
        field[i] = (field_NED[0] * north_vector[i] + field_NED[1] * east_vector[i] + field_NED[2] * down_vector[i]) * 10e-9;
    }

}

double geomagneticInclination(double orbit_normal[3], double earth_rotation){
    // Function to approximate Earth's magnetic field in LEO using the simplified dipole model
    // The position vector should be normalized to 1

    // Compute geomagnetic north pole unit vector (approximation as spherical Earth)
    double north_pole[3] = {0.0, 0.0, 1};

    double latitude_axis[3] = {0.0, 1.0, 0.0};
    double latitude_quaternion[4];
    rotation_quaternion((90 - geomagnetic_latitude) * pi / 180, latitude_axis, latitude_quaternion);
    quaternion_rotate(latitude_quaternion, north_pole, "active");
    
    double longitude_axis[3] = {0.0, 0.0, 1.0};
    double longitude_quaternion[4];
    rotation_quaternion((geomagnetic_longitude + earth_rotation) * pi / 180, longitude_axis, longitude_quaternion);
    quaternion_rotate(longitude_quaternion, north_pole, "active");

    // Compute angle between north pole and orbit plane normal
    double angle_cosine = dot_product(north_pole, orbit_normal, 3) / vector_norm(orbit_normal, 3) / vector_norm(north_pole, 3);
    return acos(angle_cosine);

}

double calculate_density(int year, int doy, int seconds, double altitude, double g_lat, double g_long, double lst, double f107A, double f107, double ap) {
    // Function to calculate density using NRLMSISE-00 model

    struct nrlmsise_input input;
    struct nrlmsise_flags flags;
    struct nrlmsise_output output;

    // Setting up input values
    input.doy = doy;
    input.year = year; /* without effect */
    input.sec = seconds;
    input.alt = altitude / 1000.0;
    input.g_lat = g_lat;
    input.g_long = g_long;
    input.lst = lst;
    input.f107A = f107A;
    input.f107 = f107;
    input.ap = ap;

    // Setting up flags
    flags.switches[0] = 0; // This flag is always 0
    for (int i = 1; i < 24; i++)
        flags.switches[i] = 1; // All other flags are set to 1

    // Calling gtd7 function
    gtd7(&input, &flags, &output);

    // Returning density
    return output.d[5] * 1000.0; // Density is at index 5 in output.d array
}