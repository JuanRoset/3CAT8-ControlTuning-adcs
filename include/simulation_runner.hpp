#ifndef SIMULATION_RUNNER_HPP
#define SIMULATION_RUNNER_HPP

struct Gains {
    double Kp[3];
    double Kd[3];
};

struct SimulationResult {
    double cost;
    double final_attitude_error;
    double max_angular_rate;
    double total_power;
    bool stable;
};

SimulationResult run_simulation(const Gains& gains,
                                bool write_output,
                                bool run_visualization);

#endif
