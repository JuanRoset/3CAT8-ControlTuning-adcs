cmake_minimum_required(VERSION 3.0)
project(AttitudeSimulator)

# Set C++ standard (e.g., C++11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add executable and list source files
add_executable(main
    src/main.cpp
    src/actuator_functions.cpp
    src/algebraic_functions.c
    src/control_functions.c
    src/disturbance_torques.cpp
    src/orbit_functions.cpp
    src/simulation_utils.cpp
    src/space_environment.c
    src/surface_torques_model.cpp
    external/src/geomag.c
    external/src/nrlmsise-00.c
)

# Include directories for header files
target_include_directories(main
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/include
)

# Example: Link math library (libm on Unix-like systems)
# target_link_libraries(main PRIVATE m)

# Generate configuration files if needed
configure_file(parameters/astronomical_parameters.json ${CMAKE_BINARY_DIR}/parameters/astronomical_parameters.json COPYONLY)
configure_file(parameters/satellite_parameters.json ${CMAKE_BINARY_DIR}/parameters/satellite_parameters.json COPYONLY)