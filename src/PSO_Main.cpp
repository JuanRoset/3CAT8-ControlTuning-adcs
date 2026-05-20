#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>

// This file uses the same simulation interface as MonteCarlo_Main.cpp.
// The actual run_simulation(...) function is implemented in Main_Simulation.cpp
// and declared in simulation_runner.hpp.

#include "../include/simulation_runner.hpp"

// ================================================================
// PSO configuration
// ================================================================

constexpr int DIM = 6;  // [Kp_x, Kp_y, Kp_z, Kd_x, Kd_y, Kd_z]

// Gain ranges.
// The PSO works internally in log10-space because the gains span
// several orders of magnitude.
constexpr double KP_MIN = 1e-8;
constexpr double KP_MAX = 1e-5;
constexpr double KD_MIN = 1e-5;
constexpr double KD_MAX = 1e-1;

constexpr double LOG_KP_MIN = -8.0;
constexpr double LOG_KP_MAX = -5.0;
constexpr double LOG_KD_MIN = -5.0;
constexpr double LOG_KD_MAX = -1.0;

// Recommended starting values for this problem.
// w decreases from W_MAX to W_MIN during the optimization.
constexpr double W_MAX = 0.8;
constexpr double W_MIN = 0.4;
constexpr double C1 = 1.5;
constexpr double C2 = 1.5;

// Maximum velocity as a fraction of each variable range.
// Since we are in log10-space, this limits jumps in orders of magnitude.
constexpr double VMAX_FRACTION = 0.25;

struct Particle {
    double position[DIM];
    double velocity[DIM];

    double best_position[DIM];
    double best_cost;

    bool has_personal_best;
};

double uniform_real(double a, double b, std::mt19937& rng)
{
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

double clamp_value(double x, double x_min, double x_max)
{
    if (x < x_min) return x_min;
    if (x > x_max) return x_max;
    return x;
}

double get_lower_bound(int d)
{
    if (d < 3) return LOG_KP_MIN;
    return LOG_KD_MIN;
}

double get_upper_bound(int d)
{
    if (d < 3) return LOG_KP_MAX;
    return LOG_KD_MAX;
}

double get_vmax(int d)
{
    const double range = get_upper_bound(d) - get_lower_bound(d);
    return VMAX_FRACTION * range;
}

Gains position_to_gains(const double position[DIM])
{
    Gains gains;

    for (int i = 0; i < 3; i++) {
        gains.Kp[i] = std::pow(10.0, position[i]);
        gains.Kd[i] = std::pow(10.0, position[i + 3]);
    }

    return gains;
}

void initialize_particle(Particle& p, std::mt19937& rng)
{
    p.best_cost = std::numeric_limits<double>::infinity();
    p.has_personal_best = false;

    for (int d = 0; d < DIM; d++) {
        const double x_min = get_lower_bound(d);
        const double x_max = get_upper_bound(d);
        const double vmax = get_vmax(d);

        p.position[d] = uniform_real(x_min, x_max, rng);

        // Small random initial velocity.
        p.velocity[d] = uniform_real(-0.2 * vmax, 0.2 * vmax, rng);

        p.best_position[d] = p.position[d];
    }
}

void write_result_header(std::ofstream& file)
{
    file << "Iteration,Particle,"
         << "Kp_x,Kp_y,Kp_z,"
         << "Kd_x,Kd_y,Kd_z,"
         << "Cost,FinalAttitudeError,MaxAngularRate,TotalPower,Stable,"
         << "IsPersonalBest,IsGlobalBest\n";
}

void write_result_row(std::ofstream& file,
                      int iteration,
                      int particle_id,
                      const Gains& gains,
                      const SimulationResult& result,
                      bool is_personal_best,
                      bool is_global_best)
{
    file << iteration << ","
         << particle_id << ","
         << gains.Kp[0] << "," << gains.Kp[1] << "," << gains.Kp[2] << ","
         << gains.Kd[0] << "," << gains.Kd[1] << "," << gains.Kd[2] << ","
         << result.cost << ","
         << result.final_attitude_error << ","
         << result.max_angular_rate << ","
         << result.total_power << ","
         << result.stable << ","
         << is_personal_best << ","
         << is_global_best
         << "\n";
}

void write_best_history_header(std::ofstream& file)
{
    file << "Iteration,BestCost,"
         << "Kp_x,Kp_y,Kp_z,"
         << "Kd_x,Kd_y,Kd_z\n";
}

void write_best_history_row(std::ofstream& file,
                            int iteration,
                            double global_best_cost,
                            const double global_best_position[DIM])
{
    Gains best_gains = position_to_gains(global_best_position);

    file << iteration << ","
         << global_best_cost << ","
         << best_gains.Kp[0] << "," << best_gains.Kp[1] << "," << best_gains.Kp[2] << ","
         << best_gains.Kd[0] << "," << best_gains.Kd[1] << "," << best_gains.Kd[2]
         << "\n";
}

int main(int argc, char* argv[])
{
    int n_particles = 20;
    int n_iterations = 30;

    if (argc > 1) {
        n_particles = std::stoi(argv[1]);
    }

    if (argc > 2) {
        n_iterations = std::stoi(argv[2]);
    }

    if (n_particles <= 0) {
        std::cerr << "Error: n_particles must be positive." << std::endl;
        return 1;
    }

    if (n_iterations <= 0) {
        std::cerr << "Error: n_iterations must be positive." << std::endl;
        return 1;
    }

    std::random_device rd;
    std::mt19937 rng(rd());

    std::ofstream results_file("../data/pso_results.csv");
    if (!results_file.is_open()) {
        std::cerr << "Error opening ../data/pso_results.csv" << std::endl;
        return 1;
    }

    std::ofstream best_history_file("../data/pso_best_history.csv");
    if (!best_history_file.is_open()) {
        std::cerr << "Error opening ../data/pso_best_history.csv" << std::endl;
        return 1;
    }

    write_result_header(results_file);
    write_best_history_header(best_history_file);

    Particle* swarm = new Particle[n_particles];

    for (int i = 0; i < n_particles; i++) {
        initialize_particle(swarm[i], rng);
    }

    double global_best_position[DIM] = {0.0};
    double global_best_cost = std::numeric_limits<double>::infinity();
    bool has_global_best = false;

    SimulationResult global_best_result = {
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
        false
    };

    // ================================================================
    // Main PSO loop
    // ================================================================
    for (int iter = 0; iter < n_iterations; iter++) {

        double w = W_MAX;
        if (n_iterations > 1) {
            w = W_MAX - (W_MAX - W_MIN) * double(iter) / double(n_iterations - 1);
        }

        std::cout << "\n===== PSO ITERATION "
                  << iter + 1 << "/" << n_iterations
                  << " | w = " << w
                  << " =====" << std::endl;

        // ------------------------------------------------------------
        // Evaluate all particles
        // ------------------------------------------------------------
        for (int i = 0; i < n_particles; i++) {

            Gains gains = position_to_gains(swarm[i].position);

            // Do not write full time histories or launch visualizers during PSO.
            SimulationResult result = run_simulation(gains, false, false);

            bool valid_result = result.stable && std::isfinite(result.cost);
            bool is_personal_best = false;
            bool is_global_best = false;

            if (valid_result && result.cost < swarm[i].best_cost) {
                swarm[i].best_cost = result.cost;
                swarm[i].has_personal_best = true;

                for (int d = 0; d < DIM; d++) {
                    swarm[i].best_position[d] = swarm[i].position[d];
                }

                is_personal_best = true;
            }

            if (valid_result && result.cost < global_best_cost) {
                global_best_cost = result.cost;
                global_best_result = result;
                has_global_best = true;

                for (int d = 0; d < DIM; d++) {
                    global_best_position[d] = swarm[i].position[d];
                }

                is_global_best = true;
            }

            write_result_row(results_file,
                             iter + 1,
                             i + 1,
                             gains,
                             result,
                             is_personal_best,
                             is_global_best);

            std::cout << "[Particle " << i + 1 << "/" << n_particles << "] "
                      << "Cost = " << result.cost
                      << " | Stable = " << result.stable;

            if (is_global_best) {
                std::cout << "  <-- new global best";
            }
            else if (is_personal_best) {
                std::cout << "  <-- new personal best";
            }

            std::cout << std::endl;
        }

        if (has_global_best) {
            write_best_history_row(best_history_file,
                                   iter + 1,
                                   global_best_cost,
                                   global_best_position);

            std::cout << "Current best cost: " << global_best_cost << std::endl;
        }
        else {
            std::cout << "No stable particle found yet." << std::endl;
        }

        // ------------------------------------------------------------
        // Update all particles
        // ------------------------------------------------------------
        for (int i = 0; i < n_particles; i++) {

            for (int d = 0; d < DIM; d++) {
                const double r1 = uniform_real(0.0, 1.0, rng);
                const double r2 = uniform_real(0.0, 1.0, rng);

                const double x = swarm[i].position[d];
                const double v = swarm[i].velocity[d];

                double cognitive_term = 0.0;
                if (swarm[i].has_personal_best) {
                    cognitive_term = C1 * r1 * (swarm[i].best_position[d] - x);
                }

                double social_term = 0.0;
                if (has_global_best) {
                    social_term = C2 * r2 * (global_best_position[d] - x);
                }

                double new_velocity = w * v + cognitive_term + social_term;

                const double vmax = get_vmax(d);
                new_velocity = clamp_value(new_velocity, -vmax, vmax);

                double new_position = x + new_velocity;

                const double x_min = get_lower_bound(d);
                const double x_max = get_upper_bound(d);

                // If a particle hits a boundary, clamp it and damp the velocity
                // so it does not keep bouncing too aggressively.
                if (new_position < x_min) {
                    new_position = x_min;
                    new_velocity *= -0.5;
                }
                else if (new_position > x_max) {
                    new_position = x_max;
                    new_velocity *= -0.5;
                }

                swarm[i].position[d] = new_position;
                swarm[i].velocity[d] = new_velocity;
            }
        }
    }

    results_file.close();
    best_history_file.close();

    std::cout << "\n===== PSO BEST RESULT =====" << std::endl;

    if (has_global_best) {
        Gains best_gains = position_to_gains(global_best_position);

        std::cout << "Cost: " << global_best_cost << std::endl;
        std::cout << "Stable: " << global_best_result.stable << std::endl;
        std::cout << "Final attitude error: " << global_best_result.final_attitude_error << std::endl;
        std::cout << "Max angular rate: " << global_best_result.max_angular_rate << std::endl;
        std::cout << "Total power: " << global_best_result.total_power << std::endl;

        std::cout << "Kp = ["
                  << best_gains.Kp[0] << ", "
                  << best_gains.Kp[1] << ", "
                  << best_gains.Kp[2] << "]" << std::endl;

        std::cout << "Kd = ["
                  << best_gains.Kd[0] << ", "
                  << best_gains.Kd[1] << ", "
                  << best_gains.Kd[2] << "]" << std::endl;

        // Rerun the best case once, now saving the full output and launching visualizers.
        run_simulation(best_gains, true, true);
    }
    else {
        std::cerr << "No stable gains found. Final simulation skipped." << std::endl;
    }

    delete[] swarm;

    return 0;
}
