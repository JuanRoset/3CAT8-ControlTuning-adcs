// Include the header files with the necessary functions

#include "../include/simulation_utils.hpp"

// Define functions

void satellite_fetch(double inertia[3][3], double & inertia_principal_min, std::string & control_mode, double Kp[3], double Kd[3], double & RW_relaxation, double & Ku, bool & active_force_control, bool & RW_unloading, int& configuration, std::string filePath){

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

    // Fetch the desired satellite configuration
    configuration = root["SatelliteParameters"]["configuration"].get<int>();

    // Fetch the correct inertia and control gain values
    inertia_principal_min = root["SatelliteParameters"]["Ipmin"].get<double>();
    switch(configuration){
        case 0:
        inertia[0][0] = root["PhysicalParameters0"]["Ixx"].get<double>();
        inertia[1][1] = root["PhysicalParameters0"]["Iyy"].get<double>();
        inertia[2][2] = root["PhysicalParameters0"]["Izz"].get<double>();
        inertia[1][0] = root["PhysicalParameters0"]["Ixy"].get<double>();
        inertia[2][0] = root["PhysicalParameters0"]["Ixz"].get<double>();
        inertia[2][1] = root["PhysicalParameters0"]["Iyz"].get<double>();
        inertia[0][1] = root["PhysicalParameters0"]["Ixy"].get<double>();
        inertia[0][2] = root["PhysicalParameters0"]["Ixz"].get<double>();
        inertia[1][2] = root["PhysicalParameters0"]["Iyz"].get<double>();
        Kp[0] = root["ControlParameters"]["Kpx0"].get<double>();
        Kp[1] = root["ControlParameters"]["Kpy0"].get<double>();
        Kp[2] = root["ControlParameters"]["Kpz0"].get<double>();
        Kd[0] = root["ControlParameters"]["Kdx0"].get<double>();
        Kd[1] = root["ControlParameters"]["Kdy0"].get<double>();
        Kd[2] = root["ControlParameters"]["Kdz0"].get<double>();
        break;
        case 1:
        inertia[0][0] = root["PhysicalParameters1"]["Ixx"].get<double>();
        inertia[1][1] = root["PhysicalParameters1"]["Iyy"].get<double>();
        inertia[2][2] = root["PhysicalParameters1"]["Izz"].get<double>();
        inertia[1][0] = root["PhysicalParameters1"]["Ixy"].get<double>();
        inertia[2][0] = root["PhysicalParameters1"]["Ixz"].get<double>();
        inertia[2][1] = root["PhysicalParameters1"]["Iyz"].get<double>();
        inertia[0][1] = root["PhysicalParameters1"]["Ixy"].get<double>();
        inertia[0][2] = root["PhysicalParameters1"]["Ixz"].get<double>();
        inertia[1][2] = root["PhysicalParameters1"]["Iyz"].get<double>();
        Kp[0] = root["ControlParameters"]["Kpx1"].get<double>();
        Kp[1] = root["ControlParameters"]["Kpy1"].get<double>();
        Kp[2] = root["ControlParameters"]["Kpz1"].get<double>();
        Kd[0] = root["ControlParameters"]["Kdx1"].get<double>();
        Kd[1] = root["ControlParameters"]["Kdy1"].get<double>();
        Kd[2] = root["ControlParameters"]["Kdz1"].get<double>();
        break;
        case 2:
        inertia[0][0] = root["PhysicalParameters2"]["Ixx"].get<double>();
        inertia[1][1] = root["PhysicalParameters2"]["Iyy"].get<double>();
        inertia[2][2] = root["PhysicalParameters2"]["Izz"].get<double>();
        inertia[1][0] = root["PhysicalParameters2"]["Ixy"].get<double>();
        inertia[2][0] = root["PhysicalParameters2"]["Ixz"].get<double>();
        inertia[2][1] = root["PhysicalParameters2"]["Iyz"].get<double>();
        inertia[0][1] = root["PhysicalParameters2"]["Ixy"].get<double>();
        inertia[0][2] = root["PhysicalParameters2"]["Ixz"].get<double>();
        inertia[1][2] = root["PhysicalParameters2"]["Iyz"].get<double>();
        Kp[0] = root["ControlParameters"]["Kpx2"].get<double>();
        Kp[1] = root["ControlParameters"]["Kpy2"].get<double>();
        Kp[2] = root["ControlParameters"]["Kpz2"].get<double>();
        Kd[0] = root["ControlParameters"]["Kdx2"].get<double>();
        Kd[1] = root["ControlParameters"]["Kdy2"].get<double>();
        Kd[2] = root["ControlParameters"]["Kdz2"].get<double>();
        break;
    }

    // Fetch the general control parameters
    RW_unloading = root["SatelliteParameters"]["unloading"].get<bool>();
    active_force_control = root["SatelliteParameters"]["afc"].get<bool>();
    Ku = root["ControlParameters"]["Ku"].get<double>();
    RW_relaxation = root["ControlParameters"]["RW_relax"].get<double>();

    // Switch the control modes
    int control_idx = root["SatelliteParameters"]["mission"].get<int>();
    switch(control_idx){
        case 0: control_mode = "no_control";break;
        case 1: control_mode = "ideal"; break;
        case 2: control_mode = "detumble"; break;
        case 3: control_mode = "mtq_only"; break;
        case 4: control_mode = "hybrid"; break;
        default:
            std::cout << "Unexpected mission mode: " << control_idx << std::endl;
            break;
    }

}

void simulation_fetch(double & simulation_orbits, double & delta_time, bool & propagate_orbit, bool & gravitational_disturbance, bool & aerodynamic_disturbance, bool & solarPressure_disturbance, double euler_axis[3], double & theta, std::string filePath){
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

    // Fetch the desired satellite configuration
    simulation_orbits = root["SimulationParameters"]["orbit_number"].get<double>();
    delta_time = root["SimulationParameters"]["delta_time"].get<double>();
    propagate_orbit = root["SimulationParameters"]["propagate_orbit"].get<bool>();

    // Fetch the disturbance configurations
    gravitational_disturbance = root["SimulationParameters"]["gravitational_disturbance"].get<bool>();
    aerodynamic_disturbance = root["SimulationParameters"]["aerodynamic_disturbance"].get<bool>();
    solarPressure_disturbance = root["SimulationParameters"]["solarPressure_disturbance"].get<bool>();

    // Fetch the initial disturbances
    euler_axis[0] = root["SimulationParameters"]["axis_disturbance_x"].get<double>();
    euler_axis[1] = root["SimulationParameters"]["axis_disturbance_y"].get<double>();
    euler_axis[2] = root["SimulationParameters"]["axis_disturbance_z"].get<double>();
    theta = root["SimulationParameters"]["angle_disturbance_deg"].get<double>() * pi / 180.0;
}


void write_orbit_state(double delta_time, orbitalElements elements, sunActivity indices, std::string orbitPath, std::string wmm_path, const double iteration_steps){
    // Function for calling the orbit propagation at all times and writing the results on the orbit file

    // Preallocate and store propagation of orbit positions for simulation time
    double velocity[3] = {0.0};
    double position[3] = {0.0};
    double quaternion_rsw[4] = {0.0};
    double current_time = 0.0;
    double long_lat[2] = {0.0};
    double magnetic_field[3] = {0.0};
    double density = 0.0;
    double magnetic_inclination = 0.0;
    double sun_position[3] = {0.0};
    double eclipse = 0.0;
    
    // Open the orbit file for writing position and orbit quaternion
    std::ofstream orbitFileOutput(orbitPath);

    // Check if the file is opened successfully
    if (!orbitFileOutput.is_open()) {
        std::cerr << "Error opening the file for writing: " << orbitPath << std::endl;
    }

    // Write output data header
    orbitFileOutput << "Time" << ","
                    << "VelocityX,VelocityY,VelocityZ" << ","
                    << "PositionX,PositionY,PositionZ" << ","
                    << "QuaternionW,QuaternionX,QuaternionY,QuaternionZ" << "," 
                    << "Longitude,Latitude" << ","
                    << "MagX,MagY,MagZ" << ","
                    << "SunX,SunY,SunZ" << ","
                    << "Eclipse" << ","
                    << "MagneticInclination" << ","
                    << "rho" << "\n";

    // Print status
    std::cout << "Propagating orbit" << std::endl;

    for(int t = 0; t < iteration_steps; t++){
        // Update time, position and velocity
        current_time = t * delta_time;
        orbit_propagate(elements, current_time, position, velocity, quaternion_rsw, long_lat, magnetic_field, sun_position, eclipse, magnetic_inclination, density, wmm_path);


        // Write the current time and data to the output datafile
        orbitFileOutput << current_time << "," 
                    << velocity[0] << "," << velocity[1] << "," << velocity[2] << "," 
                    << position[0] << "," << position[1] << "," << position[2] << "," 
                    << quaternion_rsw[0] << "," << quaternion_rsw[1] << "," << quaternion_rsw[2] << "," << quaternion_rsw[3] << ","
                    << long_lat[0] << "," << long_lat[1] << ","
                    << magnetic_field[0] << "," << magnetic_field[1] << "," << magnetic_field[2] << ","
                    << sun_position[0] << "," << sun_position[1] << "," << sun_position[2] << ","
                    << eclipse << ","
                    << magnetic_inclination << ","
                    << density << "\n";
    }

    // Close output file
    orbitFileOutput.close();
}

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
                      std::vector<double>& air_density){
    
    // Function reading the orbit state from the orbit file and storing it on the given vectors

    // Open the orbit file for writing position and orbit quaternion
    std::ifstream orbitFileInput(orbitPath);

    // Check if the file is opened successfully
    if (!orbitFileInput.is_open()) {
        std::cerr << "Error opening the file for reading: " << orbitPath << std::endl;
    }

    // Read the header line to ignore it
    std::string header;
    std::getline(orbitFileInput, header);

    // Read and process the rest of the lines
    std::string line;

    // Initialize counter
    int t = 0;

    // Read through the csv and store data in vector
    while (std::getline(orbitFileInput, line) && (t - 1) * delta_time < simulation_time) {
        std::istringstream iss(line);
        std::vector<double> tokens;
        double token;
        char delimiter;
        while (iss >> token >> delimiter) {
            tokens.push_back(token);
        }

        // If there's a last token without a delimiter
        iss >> token;
        tokens.push_back(token);

        // Store the data into the corresponding vectors
        time_array[t] = tokens[0];
        velocity_array[0][t] = tokens[1];
        velocity_array[1][t] = tokens[2];
        velocity_array[2][t] = tokens[3];
        position_array[0][t] = tokens[4];
        position_array[1][t] = tokens[5];
        position_array[2][t] = tokens[6];
        quaternion_rsw_array[0][t] = tokens[7];
        quaternion_rsw_array[1][t] = tokens[8];
        quaternion_rsw_array[2][t] = tokens[9];
        quaternion_rsw_array[3][t] = tokens[10];
        longitude_latitude[0][t] = tokens[11];
        longitude_latitude[1][t] = tokens[12];
        magnetic_field_array[0][t] = tokens[13];
        magnetic_field_array[1][t] = tokens[14];
        magnetic_field_array[2][t] = tokens[15];
        sun_position[0][t] = tokens[16];
        sun_position[1][t] = tokens[17];
        sun_position[2][t] = tokens[18];
        eclipse[t] = tokens[19];
        geomagnetic_inclination[t] = tokens[20];
        air_density[t] = tokens[21];

        // Update counter
        t++;
    }
}