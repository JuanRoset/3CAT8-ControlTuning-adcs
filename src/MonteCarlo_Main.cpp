#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <vector>
#include <chrono>

// IMPORTANT:
// This file assumes that Gains, SimulationResult and run_simulation(...)
// are declared in a small header, for example:
//
//   #include "../include/simulation_runner.hpp"
//
// If you have not created that header yet, see the chat response.
// Replace the include below with your actual path.

#include "../include/simulation_runner.hpp"

constexpr double KP_MIN = 1e-8;
constexpr double KP_MAX = 1e-5;
constexpr double KD_MIN = 1e-5;
constexpr double KD_MAX = 1e-1;

constexpr double SIGMA_LOCAL = 0.4;

// Fraction of simulations used for the first global exploration phase and number of best simulations for local refinement.
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

Gains sample_global_gains(std::mt19937& rng)
{
    Gains gains;

    for (int i = 0; i < 3; i++) {
        gains.Kp[i] = sample_log_uniform(KP_MIN, KP_MAX, rng);
        gains.Kd[i] = sample_log_uniform(KD_MIN, KD_MAX, rng);
    }

    return gains;
}

Gains sample_local_gains(const Gains& best_gains, std::mt19937& rng)
{
    Gains gains;

    for (int i = 0; i < 3; i++) {
        gains.Kp[i] = perturb_positive(best_gains.Kp[i], KP_MIN, KP_MAX, rng);
        gains.Kd[i] = perturb_positive(best_gains.Kd[i], KD_MIN, KD_MAX, rng);
    }

    return gains;
}

struct Candidate {
    Gains gains;
    SimulationResult result;
};

void update_top_candidates(std::vector<Candidate>& top_candidates,
                           const Gains& gains,
                           const SimulationResult& result,
                           int max_candidates)
{
    if (!result.stable) {
        return;
    }

    Candidate candidate;
    candidate.gains = gains;
    candidate.result = result;

    top_candidates.push_back(candidate);

    std::sort(top_candidates.begin(), top_candidates.end(),
              [](const Candidate& a, const Candidate& b) {
                  return a.result.cost < b.result.cost;
              });

    if (static_cast<int>(top_candidates.size()) > max_candidates) {
        top_candidates.pop_back();
    }
}

void write_result_row(std::ofstream& file,
                      int simulation_id,
                      const std::string& phase,
                      const Gains& gains,
                      const SimulationResult& result,
                      bool is_best)
{
    file << simulation_id << ","
         << phase << ","
         << gains.Kp[0] << "," << gains.Kp[1] << "," << gains.Kp[2] << ","
         << gains.Kd[0] << "," << gains.Kd[1] << "," << gains.Kd[2] << ","
         << result.cost << ","
         << result.final_attitude_error << ","
         << result.max_angular_rate << ","
         << result.total_power << ","
         << result.stable << ","
         << is_best
         << "\n";
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

    auto now = std::chrono::high_resolution_clock::now()
               .time_since_epoch()
               .count();

    std::seed_seq seed_seq{
        static_cast<unsigned int>(now),
        static_cast<unsigned int>(now >> 32),
        rd(),
        rd(),
        rd(),
        rd()
    };

    std::mt19937 rng(seed_seq);

    std::cout << "Monte Carlo random seed generated from time and random_device." << std::endl;

    const int n_exploration = std::max(1, static_cast<int>(EXPLORATION_FRACTION * n_simulations));
    const int n_local = n_simulations - n_exploration;

    const int n_best_for_local = std::min(N_BEST_FOR_LOCAL, n_local);

    if (n_local > 0 && n_best_for_local <= 0) {
        std::cerr << "Error: N_BEST_FOR_LOCAL must be positive when local simulations exist." << std::endl;
        return 1;
    }

    std::ofstream results_file("../data/montecarlo_results.csv");

    if (!results_file.is_open()) {
        std::cerr << "Error opening ../data/montecarlo_results.csv" << std::endl;
        return 1;
    }

    results_file << "Simulation,Phase,"
                 << "Kp_x,Kp_y,Kp_z,"
                 << "Kd_x,Kd_y,Kd_z,"
                 << "Cost,FinalAttitudeError,MaxAngularRate,TotalPower,Stable,IsBest\n";

    std::vector<Candidate> top_candidates;

    Gains best_gains = {};
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

        Gains test_gains = sample_global_gains(rng);

        // During Monte Carlo, do not write the full time history and do not run visualization.
        SimulationResult result = run_simulation(test_gains, false, false);

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
            const Gains& center_gains = top_candidates[candidate_index].gains;

            Gains test_gains = sample_local_gains(center_gains, rng);

            SimulationResult result = run_simulation(test_gains, false, false);

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

    std::cout << "\n===== BEST RESULT =====" << std::endl;
    std::cout << "Cost: " << best_result.cost << std::endl;
    std::cout << "Stable: " << best_result.stable << std::endl;

    std::cout << "Kp = ["
              << best_gains.Kp[0] << ", "
              << best_gains.Kp[1] << ", "
              << best_gains.Kp[2] << "]" << std::endl;

    std::cout << "Kd = ["
              << best_gains.Kd[0] << ", "
              << best_gains.Kd[1] << ", "
              << best_gains.Kd[2] << "]" << std::endl;

    // Optional: rerun the best case once, now saving the full output and launching the visualizers.
    if (best_result.stable) {
        run_simulation(best_gains, true, true);
    }
    else {
        std::cerr << "No stable gains found. Final simulation skipped." << std::endl;
    }

    return 0;
}
