#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <vector>

// This header must declare:
//   struct CascadeGains {
//       double Kp_phi[3];
//       double Kd_phi[3];
//       double Kp_omega[3];
//       double Kd_omega[3];
//   };
//
//   struct SimulationResult {
//       double cost;
//       double final_attitude_error;
//       double max_angular_rate;
//       double total_power;
//       bool stable;
//   };
//
//   SimulationResult run_simulation_cascade(const CascadeGains& gains,
//                                          bool write_output,
//                                          bool run_visualization);
//
// If you decide to keep a single header, replace this include by your actual one.
#include "../include/simulation_runner_cascade.hpp"

constexpr double KP_MIN = 1e-8;
constexpr double KP_MAX = 1e-5;
constexpr double KD_MIN = 1e-5;
constexpr double KD_MAX = 1e-1;

constexpr double SIGMA_LOCAL = 0.4;

// Fraction of simulations used for the first global exploration phase.
// The rest are used for local refinement around the current best gains.
constexpr double EXPLORATION_FRACTION = 0.8;
constexpr int N_BEST_FOR_LOCAL = 3;

double sample_log_uniform(double k_min, double k_max, std::mt19937& rng)
{
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);

    const double x = uniform01(rng);
    const double log_k = std::log10(k_min) + x * (std::log10(k_max) - std::log10(k_min));

    return std::pow(10.0, log_k);
}

double perturb_positive(double current_best, double k_min, double k_max, std::mt19937& rng)
{
    std::normal_distribution<double> normal(1.0, SIGMA_LOCAL);

    double g = normal(rng);

    // The Gaussian may very rarely return a negative multiplier.
    // Negative gains are not physically useful here, so resample.
    int attempts = 0;
    while (g <= 0.0 && attempts < 20) {
        g = normal(rng);
        attempts++;
    }

    if (g <= 0.0) {
        g = 1.0;
    }

    const double candidate = current_best * g;

    if (candidate < k_min) {
        return k_min;
    }
    if (candidate > k_max) {
        return k_max;
    }
    return candidate;
}

CascadeGains sample_global_gains(std::mt19937& rng)
{
    CascadeGains gains = {};

    for (int i = 0; i < 3; i++) {
        // First cascade step: attitude/angle-error loop
        gains.Kp_phi[i] = sample_log_uniform(KP_MIN, KP_MAX, rng);
        gains.Kd_phi[i] = sample_log_uniform(KD_MIN, KD_MAX, rng);

        // Second cascade step: angular-rate loop
        gains.Kp_omega[i] = sample_log_uniform(KP_MIN, KP_MAX, rng);
        gains.Kd_omega[i] = sample_log_uniform(KD_MIN, KD_MAX, rng);
    }

    return gains;
}

CascadeGains sample_local_gains(const CascadeGains& best_gains, std::mt19937& rng)
{
    CascadeGains gains = {};

    for (int i = 0; i < 3; i++) {
        // First cascade step: attitude/angle-error loop
        gains.Kp_phi[i] = perturb_positive(best_gains.Kp_phi[i], KP_MIN, KP_MAX, rng);
        gains.Kd_phi[i] = perturb_positive(best_gains.Kd_phi[i], KD_MIN, KD_MAX, rng);

        // Second cascade step: angular-rate loop
        gains.Kp_omega[i] = perturb_positive(best_gains.Kp_omega[i], KP_MIN, KP_MAX, rng);
        gains.Kd_omega[i] = perturb_positive(best_gains.Kd_omega[i], KD_MIN, KD_MAX, rng);
    }

    return gains;
}

void write_result_row(std::ofstream& file,
                      int simulation_id,
                      const std::string& phase,
                      const CascadeGains& gains,
                      const SimulationResult& result,
                      bool is_best)
{
    file << simulation_id << ","
         << phase << ","

         << gains.Kp_phi[0] << "," << gains.Kp_phi[1] << "," << gains.Kp_phi[2] << ","
         << gains.Kd_phi[0] << "," << gains.Kd_phi[1] << "," << gains.Kd_phi[2] << ","

         << gains.Kp_omega[0] << "," << gains.Kp_omega[1] << "," << gains.Kp_omega[2] << ","
         << gains.Kd_omega[0] << "," << gains.Kd_omega[1] << "," << gains.Kd_omega[2] << ","

         << result.cost << ","
         << result.final_attitude_error << ","
         << result.max_angular_rate << ","
         << result.total_power << ","
         << result.stable << ","
         << is_best
         << "\n";
}

void print_gain_vector(const std::string& name, const double gains[3])
{
    std::cout << name << " = ["
              << gains[0] << ", "
              << gains[1] << ", "
              << gains[2] << "]" << std::endl;
}

struct CascadeCandidate {
    CascadeGains gains;
    SimulationResult result;
};

void update_top_candidates(std::vector<CascadeCandidate>& top_candidates,
                           const CascadeGains& gains,
                           const SimulationResult& result,
                           int max_candidates)
{
    if (!result.stable) {
        return;
    }

    CascadeCandidate candidate;
    candidate.gains = gains;
    candidate.result = result;

    top_candidates.push_back(candidate);

    std::sort(top_candidates.begin(), top_candidates.end(),
              [](const CascadeCandidate& a, const CascadeCandidate& b) {
                  return a.result.cost < b.result.cost;
              });

    if (static_cast<int>(top_candidates.size()) > max_candidates) {
        top_candidates.pop_back();
    }
}

int main(int argc, char* argv[])
{
    int n_simulations = 100;

    if (argc > 1) {
        n_simulations = std::stoi(argv[1]);
    }

    if (n_simulations <= 0) {
        std::cerr << "Error: n_simulations must be positive." << std::endl;
        return 1;
    }

    std::random_device rd;
    std::mt19937 rng(rd());

    const int n_exploration = std::max(1, static_cast<int>(EXPLORATION_FRACTION * n_simulations));
    const int n_local = n_simulations - n_exploration;

    const int n_best_for_local = std::min(N_BEST_FOR_LOCAL, n_local);

    if (n_local > 0 && n_best_for_local <= 0) {
        std::cerr << "Error: N_BEST_FOR_LOCAL must be positive when local simulations exist." << std::endl;
        return 1;
    }

    std::ofstream results_file("../data/montecarlo_cascade_results.csv");

    if (!results_file.is_open()) {
        std::cerr << "Error opening ../data/montecarlo_cascade_results.csv" << std::endl;
        return 1;
    }

    results_file << "Simulation,Phase,"
                 << "Kp_phi_x,Kp_phi_y,Kp_phi_z,"
                 << "Kd_phi_x,Kd_phi_y,Kd_phi_z,"
                 << "Kp_omega_x,Kp_omega_y,Kp_omega_z,"
                 << "Kd_omega_x,Kd_omega_y,Kd_omega_z,"
                 << "Cost,FinalAttitudeError,MaxAngularRate,TotalPower,Stable,IsBest\n";

    std::vector<CascadeCandidate> top_candidates;
    
    CascadeGains best_gains = {};
    SimulationResult best_result = {
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
        false
    };

    int simulation_id = 0;

    // ================================================================
    // Phase 1: global log-uniform exploration
    // ================================================================
    for (int i = 0; i < n_exploration; i++) {
        simulation_id++;

        CascadeGains test_gains = sample_global_gains(rng);

        // During Monte Carlo, do not write the full time history and do not run visualization.
        SimulationResult result = run_simulation_cascade(test_gains, false, false);

        update_top_candidates(top_candidates, test_gains, result, n_best_for_local);

        bool is_best = false;
        if (result.stable && result.cost < best_result.cost) {
            best_gains = test_gains;
            best_result = result;
            is_best = true;
        }

        write_result_row(results_file, simulation_id, "exploration", test_gains, result, is_best);

        std::cout << "[Exploration " << simulation_id << "/" << n_simulations << "] "
                  << "Cost = " << result.cost;

        if (is_best) {
            std::cout << "  <-- new best";
        }

        std::cout << std::endl;
    }

    // ================================================================
    // Phase 2: local Gaussian refinement around the current best gains
    // ================================================================
    if (!top_candidates.empty() && n_local > 0) {

        const int active_candidates = std::min(static_cast<int>(top_candidates.size()), n_best_for_local);

        for (int i = 0; i < n_local; i++) {
            simulation_id++;

            const int candidate_index = i % active_candidates;
            const CascadeGains& center_gains = top_candidates[candidate_index].gains;

            CascadeGains test_gains = sample_local_gains(center_gains, rng);

            SimulationResult result = run_simulation_cascade(test_gains, false, false);

            update_top_candidates(top_candidates, test_gains, result, n_best_for_local);

            bool is_best = false;
            if (result.stable && result.cost < best_result.cost) {
                best_gains = test_gains;
                best_result = result;
                is_best = true;
            }

            std::string phase_name = "local_from_best_" + std::to_string(candidate_index + 1);

            write_result_row(results_file, simulation_id, phase_name, test_gains, result, is_best);

            std::cout << "[Local " << simulation_id << "/" << n_simulations << "] "
                    << "Center = best " << candidate_index + 1
                    << ", Cost = " << result.cost;

            if (is_best) {
                std::cout << "  <-- new best";
            }

            std::cout << std::endl;
        }
    }
    else {
        std::cerr << "No stable solution found during exploration, or no local simulations available. Local phase skipped." << std::endl;
    }

    results_file.close();

    std::cout << "\n===== BEST CASCADE RESULT =====" << std::endl;
    std::cout << "Cost: " << best_result.cost << std::endl;
    std::cout << "Stable: " << best_result.stable << std::endl;

    print_gain_vector("Kp_phi", best_gains.Kp_phi);
    print_gain_vector("Kd_phi", best_gains.Kd_phi);
    print_gain_vector("Kp_omega", best_gains.Kp_omega);
    print_gain_vector("Kd_omega", best_gains.Kd_omega);

    // Optional: rerun the best case once, now saving the full output and launching the visualizers.
    if (best_result.stable) {
        run_simulation_cascade(best_gains, true, true);
    }
    else {
        std::cerr << "No stable gains found. Final simulation skipped." << std::endl;
    }

    return 0;
}
