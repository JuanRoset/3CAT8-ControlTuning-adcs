#ifndef SIMULATION_RUNNER_CASCADE_HPP
#define SIMULATION_RUNNER_CASCADE_HPP

struct CascadeGains {
    double Kp_phi[3];
    double Kd_phi[3];
    double Kp_omega[3];
    double Kd_omega[3];
};

struct SimulationResult {
    double cost;
    double final_attitude_error;
    double max_angular_rate;
    double total_power;
    bool stable;
};

SimulationResult run_simulation_cascade(const CascadeGains& gains,
                                        bool write_output,
                                        bool run_visualization);

#endif
