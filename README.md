# ³cat-8 SITL-ADCS simulator

## Overview

This C/C++ code simulates the dynamics and environment of a satellite in Low Earth Orbit. The simulation uses quaternions and Euler's equations of rotational dynamics as the backbone for the propagation of the satellite's state. Magnetorquers are used along with a single reaction wheel about the Z-axis to apply the control torques. The disturbances that act on the satellite are also modeled in detail.

The core of the code has been programmed using C/C++ both for speed as well as for direct translation to the OBC software. This allows for the core functionalities of the final code to be tested in simulation.

## Dependencies

- C++ compiler (supporting C++11 or later)
- Python (for visualization using visualization scripts)
    - pygame
    - sys
    - pandas
    - pyquaternion
- [nlohmann/json](https://github.com/nlohmann/json): JSON library for C++
- [WMM 2020](https://www.ncei.noaa.gov/products/world-magnetic-model): An implementation of the World Magnetic Model in C.
- [NRLMSISE-00](https://es.mathworks.com/matlabcentral/fileexchange/56253-nrlmsise-00-atmosphere-model): Atmosphere Model

## Setup

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/nanosatlab/3cat-8-adcs-SITL.git
   cd 3cat-8-adcs-SITL
   ```

2. **Compile the Code:**

   ```bash
   mkdir build
   cd build
   cmake ..
   cmake --build .
   ```

The _parameters/satellite_parameters.json_ file can be used for changing the physical parameters of the satellite, its orbital elements, simulation parameters, control gains, and more before executing the code.

## Usage

Run the compiled executable `main`:

   ```bash
   cd build/Debug
   ./main.exe
   ```

If the _propagate_orbit_ parameter is set to true in the JSON file, the code will propagate the orbit with the specified orbital elements and store the orbital state in the `data/orbit_file.csv` file. If this is set to false, the orbital parameters from previous code executions will be read from this file. Afterwards, the simulation results will be stored in the `data/simulation_output.csv` file for posterior analysis.

## Configuration

The configuration file (_parameters/satellite_parameters.json_) includes various parameters grouped under different categories. Some of the key parameters are:

- **Simulation Parameters:**
  - `start_epoch`: Simulation start time
  - `delta_time`: Time step for the simulation
  - `orbit_number`: Number of orbits to simulate
  - `propagate_orbit`: Boolean to indicate if orbit propagation is needed
  - `gravitational_disturbance`: Enable gravitational disturbance
  - `aerodynamic_disturbance`: Enable aerodynamic disturbance
  - `solarPressure_disturbance`: Enable solar pressure disturbance

- **Satellite Parameters:**
  - `configuration`: Satellite configuration index
  - `mission`: Mission profile index
  - `unloading`: Boolean to enable reaction wheel unloading
  - `afc`: Boolean to enable automatic frequency control
  - `Ipmin`: Minimum moment of inertia parameter

- **Control Parameters:**
  - `Kpxi`, `Kpyi`, `Kpzi`: Proportional gains for control for the $i^{th}$ mission profile
  - `Kdxi`, `Kdyi`, `Kdzi`: Derivative gains for control for the $i^{th}$ mission profile
  - `Ku`: Control gain for reaction wheel unloading
  - `RW_relax`: Relaxation factor for reaction wheel

- **Orbital Elements:**
  - `epoch`, `inclination`, `right_ascension`, `eccentricity`, `argument of perigee`, `mean_anomaly`, `mean_motion`: Keplerian orbital elements

- **Physical Parameters:**
  - `Ixx`, `Iyy`, `Izz`, `Ixy`, `Ixz`, `Iyz`: Moments of inertia

## Output

The simulation results are stored in a CSV file (`data/simulation_output.csv`). The columns include time, angular rates, and quaternion components.

## Visualization

After running the simulation, Python scripts in the `visualization` folder can be used to visualize the results. Each script is designed to provide different insights into the simulation data:

1. **Attitude 3D Visualizer**

   Animates the 3D orientation of the satellite across time in 3D.

   ```bash
   python visualization/3D_visualizer.py
   ```

2. **Control Visualizer**

   Plots all of the data related to the performance of the ADCS of the satellite.

   ```bash
   python visualization/control_visualizer.py
   ```

3. **Orbit Visualizer**

   Plots the 3D orbit of the satellite and its corresponding ground track.

   ```bash
   python visualization/orbit_visualizer.py
   ```

### Dependencies for Visualization

Ensure you have the following Python libraries installed for the visualization scripts:

```bash
pip install pygame pandas pyquaternion
```

## Code Structure

- `main.cpp`: Main program file
- `actuator_functions.cpp`: Contains functions related to actuator modeling
- `algebraic_functions.c`: Contains basic algebraic and quaternion operations
- `control_functions.c`: Functions for control algorithms
- `disturbance_torques.cpp`: Models disturbance torques
- `orbit_functions.cpp`: Functions related to orbit propagation
- `simulation_utils.cpp`: Utility functions for the simulation
- `space_environment.c`: Models the space environment
- `surface_torques_model.cpp`: Models surface torques on the satellite
- `visualization/3D_visualizer.py`: Python script for 3D visualization of the simulation results
- `visualization/control_visualizer.py`: Python script for visualizing control-related data
- `visualization/orbit_visualizer.py`: Python script for visualizing orbital data

## Development

### Tasks to Complete

[ ] **Enhance Visualization**
   - Improve the Python visualization scripts to provide more insightful visualizations of the simulation results.

[ ] **Leave code as function to perform control optimization**
   - Develop a function to quantify and simplify the perturbance torques experienced in the LEO space environment.

[ ] **Improve disturbance torques computation**
   - Add efficient functions to compute the disturbance torques for simple geometries without relying on the machine-learning methods currently implemented.

[ ] **Performance Optimization**
   - Optimize the code for better performance, considering larger simulation times and more complex scenarios.

[ ] **Documentation Update**
   - Keep the code documentation up-to-date with any changes or additions.

[ ] **Unit Testing**
   - Implement unit tests for critical functions to ensure the correctness of the code.

[ ] **User Interface**
   - Develop a simple user interface to enable users to interact with the simulation parameters easily.

Feel free to pick any task or propose new ones based on the project's needs. Contributions are welcome!