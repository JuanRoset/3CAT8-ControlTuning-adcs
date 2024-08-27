// Include the header files with the necessary functions

#include "../include/orbit_functions.hpp"

// Define functions

double find_eccentric_anomaly(orbitalElements elements, double mean_anomaly, double initialGuess, double tolerance = 1e-6, int maxIterations = 100){
    // Function to perform Newton-Raphson iteration for finding the mean eccentricity

    double eccentric_anomaly[2] = {initialGuess, initialGuess};
    double residue = 0;

    // Perform iterations with Newton-Raphson applied to the mean eccentric equation
    for (int i = 0; i < maxIterations; ++i) {
        
        eccentric_anomaly[1] = eccentric_anomaly[0] - (eccentric_anomaly[0] * (1 - elements.eccentricity) - mean_anomaly) / (1 - elements.eccentricity * cos(eccentric_anomaly[0]));
        residue = abs(eccentric_anomaly[1] - eccentric_anomaly[0]);

        // Check if the solution has converged
        if (fabs(residue) < tolerance) {
            // std::cout << "Converged in " << i + 1 << " iterations." << std::endl;
            return eccentric_anomaly[1];
        }

        // If residue is not small enough, reassign trial value
        eccentric_anomaly[0] = eccentric_anomaly[1];
    }

    // Return error if solution hasn't converged
    std::cerr << "Newton-Raphson method for finding Eccentric Anomaly did not converge within the specified number of iterations (" << maxIterations << ")." << std::endl;
    std::cerr << "Check orbital_functions.h code for solution." << std::endl;
    return NAN; // Not a Number
}

void orbit_propagate(orbitalElements elements, double time, double* position, double* velocity, double* quaternion_rsw, double long_lat[2], double magnetic_field[3], double sun_position[3], double& eclipse, double& geomagnetic_inclination, double& density, std::string wmm_path){
    // Function for finding the position and velocity vector for a given orbit  

    // Find mean and eccentric anomalies, M and E
    double real_time = time + elements.orbit_start_time;
    double mean_anomaly = elements.mean_angular_motion * real_time;
    double eccentric_anomaly = find_eccentric_anomaly(elements, mean_anomaly, 1.0, 1e-6, 10);

    // Find orbit anomaly
    double theta = 2 * atan(sqrt((1 + elements.eccentricity) / (1 - elements.eccentricity)) * tan(eccentric_anomaly / 2));

    // Find position vector norm (radius)
    double radius = pow(elements.angular_momentum, 2) / (mu_Earth * (1 + elements.eccentricity * cos(theta)));

    // Find position vector in PQW frame
    double position_PQW[3] = {0};
    position_PQW[0] = radius * cos(theta);
    position_PQW[1] = radius * sin(theta);
    position_PQW[2] = 0;

    // Find velocity vector in PQW frame
    double velocity_PQW[3] = {0};
    velocity_PQW[0] = mu_Earth / elements.angular_momentum * (- sin(theta));
    velocity_PQW[1] = mu_Earth / elements.angular_momentum * (elements.eccentricity + cos(theta));
    velocity_PQW[2] = 0;

    // Perform matrix multiplication to rotate position and velocity vectors
    for(int m = 0; m < 3; m++){
        // Compute result for position and velocity
        double position_result = 0;
        double velocity_result = 0;
        for(int n = 0; n < 3; n++){
            position_result += elements.rotation_matrix[m][n] * position_PQW[n];
            velocity_result += elements.rotation_matrix[m][n] * velocity_PQW[n];
        }
        // Assign results for position and velocity in ECI
        position[m] = position_result;
        velocity[m] = velocity_result;
    }

    // Find the vectors of the RSW (Radial, Along-track, Cross-track) basis

    // The Cross track unit vector is constant and normal to the orbital plane
    double crossTrack_RSW[3] = {elements.normal_vector[0], elements.normal_vector[1], elements.normal_vector[2]};

    // The Radial unit vector is the position vector divided by its norm
    double radial_RSW[3] = {0};
    double position_norm = sqrt(pow(position[0], 2) + pow(position[1], 2) + pow(position[2], 2));
    for(int i = 0; i < 3; i++)
        radial_RSW[i] = position[i] / position_norm;

    // The Along Track unit vector is the cross product between R and S
    double alongTrack_RSW[3] = {0};
    cross_product_3(crossTrack_RSW, radial_RSW, alongTrack_RSW);

    // Construct the Direct Cosine Matrix (DCM) for going from RSW orbit to ECI
    double DCM_RSW2ECI[3][3] = {{radial_RSW[0], alongTrack_RSW[0], crossTrack_RSW[0]},
                                {radial_RSW[1], alongTrack_RSW[1], crossTrack_RSW[1]},
                                {radial_RSW[2], alongTrack_RSW[2], crossTrack_RSW[2]}};

    // Find the RSW quaternion using the conjugate of the cayley method
    rotation_matrix_to_quaternion(DCM_RSW2ECI, quaternion_rsw);

    // Find lognitude of ground track
    double absolute_longitude = atan2(position[1], position[0]);

    // Find GMST angle of Earth
    double current_day_fraction = elements.start_day_fraction + time / (24 * 3600);
    double current_julian_date = epoch_to_JD(elements.start_year, elements.start_doy, current_day_fraction);
    double GMST = JD_to_GMST(current_julian_date);

    // Find latitude of ground track
    double longitude = absolute_longitude - GMST;

    // Normalize longitude to 0, 2 * pi
    longitude = normalize_angle(longitude);

    // Shift longitude to -pi, pi
    longitude -= pi;
    
    // Find latitude angle in ECI frame
    double projection_norm = sqrt(pow(position[0], 2) + pow(position[1], 2));
    double latitude_sign = position[2] / abs(position[2]);
    double latitude = latitude_sign * acos(projection_norm / position_norm);

    // Store latitude and longitude
    long_lat[0] = longitude * 180 / pi;
    long_lat[1] = latitude * 180 / pi;

    // Compute geomagnetic inclination
    geomagnetic_inclination = geomagneticInclination(elements.normal_vector,  normalize_angle(GMST));

    // Compute magnetic field in ECI frame
    double igrf_date = double(elements.start_year + (elements.start_doy + current_day_fraction / (24.0 * 3600.0)) / 365.2425);
    igrf_magnetic_field(long_lat, position, position_norm, igrf_date, magnetic_field, wmm_path.c_str());

    // Compute atmospheric density
    int day_seconds = int(current_day_fraction * 24 * 3600);
    double altitude_sphere = position_norm - radius_Earth;
    double LMST = normalize_angle(GMST + longitude) * 180 / pi;
    density = calculate_density(elements.start_year, elements.start_doy, day_seconds, altitude_sphere, long_lat[1], long_lat[0], LMST, 150, 175, 8);

    // Compute approximate sun position
    double mean_sun_anomaly = igrf_date * 2 * pi;
    double sun_position_barycentric[3] = {cos(mean_sun_anomaly), sin(mean_sun_anomaly), 0.0};
    double sun_to_eci_rotmat[3][3] = {{1, 0.0, 0.0}, {0.0, cos(axial_tilt_Earth), -sin(axial_tilt_Earth)}, {0.0, sin(axial_tilt_Earth), cos(axial_tilt_Earth)}};
    matrix_multiply(&sun_to_eci_rotmat[0][0], &sun_position_barycentric[0], 3, 3, 1, sun_position);

    // Compute shadow and penumbra alpha angles
    double alpha_s = pi / 2 - asin((radius_Sun - radius_Earth) / AU);
    //double alpha_p = pi / 2 + asin((radius_Sun + radius_Earth) / AU);
    double theta_s = pi - alpha_s + acos(radius_Earth / position_norm);
    //double theta_p = pi - alpha_p + acos(radius_Earth / position_norm);

    // Compute the theta angle between the sun's position and the satellite's position
    double theta_satellite = acos(dot_product(position, sun_position, 3) / position_norm);
    eclipse = (theta_satellite > theta_s) ? 0.0 : 1.0;

}

void orbit_fetch(orbitalElements& elements, std::string filePath){
    // Function for reading the JSON file containing the satellite's parameters and retreiving the orbital elements

    // Open the JSON file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening file\n";
    }

    // Read the JSON data from the file
    json root;
    file >> root; // Read directly into the JSON object

    file.close();

    // Access specific elements from the JSON data
    elements.epoch = root["OrbitalElements"]["epoch"].get<double>();
    elements.start_epoch = root["SimulationParameters"]["start_epoch"].get<double>();
    elements.inclination = root["OrbitalElements"]["inclination"].get<double>();
    elements.right_ascension = root["OrbitalElements"]["right_ascension"].get<double>();
    elements.eccentricity = root["OrbitalElements"]["eccentricity"].get<double>();
    elements.argument_of_perigee = root["OrbitalElements"]["argument of perigee"].get<double>();
    elements.mean_anomaly = root["OrbitalElements"]["mean_anomaly"].get<double>();
    elements.mean_motion = root["OrbitalElements"]["mean_motion"].get<double>();

    // Compute other relevant information
    elements.period = 24 * 3600 / elements.mean_motion;
    elements.mean_angular_motion = 2 * pi / elements.period;
    elements.semimajor_axis = pow(mu_Earth / pow(elements.mean_angular_motion, 2), 1.0 / 3.0);
    elements.angular_momentum = sqrt(mu_Earth * elements.semimajor_axis * (1 - pow(elements.eccentricity, 2)));

    // Compute the year
    elements.start_year = 2000 + int(elements.start_epoch / 1000);
    double year_difference = int(elements.start_epoch / 1000) - int(elements.epoch / 1000);

    // Compute the day of year
    double epoch_doy = int(elements.epoch) - int(elements.epoch / 1000) * 1000;
    elements.start_doy = int(elements.start_epoch) - int(elements.start_epoch / 1000) * 1000;
    double day_difference = elements.start_doy - epoch_doy;

    // COmpute the fraction of day
    double fraction_epoch = elements.epoch - int(elements.epoch);
    elements.start_day_fraction = elements.start_epoch - int(elements.start_epoch);
    double fraction_difference = elements.start_day_fraction - fraction_epoch;

    // Compute the time elapsed within the satellite's orbit
    double total_delta_time = (year_difference * 365 + (day_difference + fraction_difference)) * 24 * 3600;
    double start_revolutions = int(total_delta_time / elements.period);

    elements.orbit_start_time = total_delta_time - start_revolutions * elements.period;

    // Compute the orbital rotation matrix
    double omega = elements.argument_of_perigee * pi / 180.0;
    double inclination = elements.inclination * pi / 180.0;
    double Omega = elements.right_ascension * pi / 180.0;

    elements.rotation_matrix[0][0] = + cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(inclination);
    elements.rotation_matrix[0][1] = - cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(inclination);
    elements.rotation_matrix[0][2] = + sin(Omega) * sin(inclination);
    elements.rotation_matrix[1][0] = + sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(inclination);
    elements.rotation_matrix[1][1] = - sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(inclination);
    elements.rotation_matrix[1][2] = - cos(Omega) * sin(inclination);
    elements.rotation_matrix[2][0] = + sin(omega) * sin(inclination);
    elements.rotation_matrix[2][1] = + cos(omega) * sin(inclination);
    elements.rotation_matrix[2][2] = + cos(inclination);

    // Find the unit vector normal to the orbit plane
    double normal_vector_PQW[3] = {0.0, 0.0, 1.0};
    for(int m = 0; m < 3; m++){
        // Compute result for position and velocity
        double result = 0;
        for(int n = 0; n < 3; n++){
            result += elements.rotation_matrix[m][n] * normal_vector_PQW[n];
        }
        // Assign results for position and velocity in ECI
        elements.normal_vector[m] = result;
    }

    // Normalize the normal vector
    normalize_vector(elements.normal_vector, 3);
}

void sun_fetch(sunActivity& indices, std::string filePath){
    // Function for reading the JSON file containing the satellite's parameters and retreiving the orbital elements

    // Open the JSON file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening file\n";
    }

    // Read the JSON data from the file
    json root;
    file >> root; // Read directly into the JSON object

    file.close();

    // Access specific elements from the JSON data
    indices.f107A = root["SunActivity"]["f107A"].get<double>();
    indices.f107d = root["SunActivity"]["f107"].get<double>();
    indices.ap = root["SunActivity"]["ap"].get<double>();

}

double epoch_to_JD(int year,int day_in_year, double day_fraction){
    // Function used for computing the Julian Date from a given year , doy and day fraction
    // The function returns the JD date with a double
    // References:
    // - Navipedia: https://gssc.esa.int/navipedia/index.php/Julian_Date

    // Initialize day and month
    int day;
    int month;

    // Load days in each month
    int days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Check if the year is a leap year and update days february accordingly
    // Leap year if divisible by 4 and not divisible by 100, unless divisible by 400
    if((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))
        days_in_month[1] = 29;

    // Find day and month
    for(int i = 0; i < 12; i++){
        if(day_in_year > days_in_month[i])
            // Subtract month days to day in year
            day_in_year -= days_in_month[i];
        else{
            // Assign day and month
            day = day_in_year;
            month = i + 1;

            // Break out of the loop
            break;
        }
    }

    // Compute Julian Day number (valid from 1900 - 2100)
    double julian_date = 0;
    int mod_year = year;
    int mod_month = month;
    if(month <= 2){
        mod_year = year - 1;
        mod_month = month + 12;
    }
    julian_date = int(365.25 * mod_year) + int(30.6001 * (mod_month + 1)) + day + day_fraction + 1720981.5;
    return julian_date;

}

double JD_to_GMST(double julian_date){
    // Function to convert from a given Julian Date to the corresponding GMST in radians
    // This function has been translated from a MATLAB function found in:
    // - https://es.mathworks.com/matlabcentral/fileexchange/28176-julian-date-to-greenwich-mean-sidereal-time
    // The theory behind this function can be found in:
    // - https://aa.usno.navy.mil/faq/GAST

    // Variables declaration
    double JD0, JDmin, JDmax, H, D, D0, T, GMST;
    
    // Initialize variables
    JD0 = julian_date;
    JDmin = floor(julian_date) - 0.5;
    JDmax = floor(julian_date) + 0.5;
    
    // Adjust JD0 based on julian_date to ensure proper calculation
    if (julian_date > JDmin)
        JD0 = JDmin;
    if (julian_date > JDmax)
        JD0 = JDmax;
    
    // Compute time in hours past previous midnight
    H = (julian_date - JD0) * 24;
    
    // Compute the number of days since J2000 for JD and JD0
    D = julian_date - 2451545.0;
    D0 = JD0 - 2451545.0;
    
    // Compute the number of centuries since J2000
    T = D / 36525;
    
    // Calculate GMST in hours (0h to 24h) and then convert to degrees
    GMST = fmod(6.697374558 + 0.06570982441908 * D0  + 1.00273790935 * H +
                0.000026 * (pow(T, 2)), 24) * 15;
    
    // Convert GMST to radians and return
    return GMST / 180 * pi;
}