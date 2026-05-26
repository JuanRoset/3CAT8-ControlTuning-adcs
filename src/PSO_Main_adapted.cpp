#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>

// This file uses the same simulation interface as MonteCarlo_Main.cpp.
// The actual run_simulation(...) function is implemented in Main_Simulation.cpp
// and declared in simulation_runner.hpp.

#include "../include/simulation_runner.hpp"

// ================================================================
// PSO configuration
// ================================================================

constexpr int DIM = 6;  // [Kp_x, Kp_y, Kp_z, Kd_x, Kd_y, Kd_z]

// Fixed PSO size requested for the MC500-seeded run.
constexpr int N_PARTICLES = 20;
constexpr int N_ITERATIONS = 25;

// Initial swarm composition:
// - 18 particles from the best exploration cases in montecarlo_results.csv
// - 1 particle from the best final/local exploitation case in montecarlo_results.csv
// - 1 particle from Pau Climent's TFG constants
constexpr int N_EXPLORATION_SEEDS = 18;

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

// PSO coefficients.
// w decreases from W_MAX to W_MIN during the optimization.
constexpr double W_MAX = 0.8;
constexpr double W_MIN = 0.4;
constexpr double C1 = 1.5;
constexpr double C2 = 1.5;

// Maximum velocity in log10-space.
// VMAX = 0.30 means that a particle cannot change more than roughly
// a factor 10^0.30 ≈ 2 per iteration in any gain.
constexpr double VMAX_LOG = 0.30;

// Initial velocity scale, in log10-space.
constexpr double INITIAL_VELOCITY_FRACTION = 0.20;

// Monte Carlo results generated from the 500-simulation run.
const std::string MONTECARLO_RESULTS_PATH = "../data/montecarlo_results.csv";

// Pau Climent's TFG constants from the reference table.
constexpr double PAU_TFG_KP[3] = {
    6.0e-8,   // X-axis Kp
    5.0e-6,   // Y-axis Kp
    4.0e-8    // Z-axis Kp
};

constexpr double PAU_TFG_KD[3] = {
    8.0e-5,   // X-axis Kd
    1.0e-3,   // Y-axis Kd
    4.0e-5    // Z-axis Kd
};

struct Particle {
    double position[DIM];
    double velocity[DIM];

    double best_position[DIM];
    double best_cost;

    bool has_personal_best;
};

struct MonteCarloCandidate {
    int simulation_id;
    std::string phase;
    Gains gains;
    double cost;
    double final_attitude_error;
    double max_angular_rate;
    double total_power;
    bool stable;
};

struct InitialSeed {
    std::string source;
    int simulation_id;
    double csv_cost;
    Gains gains;
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

std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;

    while (std::getline(ss, token, ',')) {
        tokens.push_back(token);
    }

    return tokens;
}

double safe_stod(const std::string& value)
{
    if (value == "inf" || value == "INF" || value == "Infinity" || value == "infinity") {
        return std::numeric_limits<double>::infinity();
    }
    return std::stod(value);
}

std::vector<MonteCarloCandidate> read_montecarlo_results(const std::string& path)
{
    std::vector<MonteCarloCandidate> candidates;

    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening Monte Carlo results file: " << path << std::endl;
        return candidates;
    }

    std::string line;
    std::getline(file, line); // header

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line(line);
        if (cols.size() < 13) continue;

        MonteCarloCandidate c;
        c.simulation_id = std::stoi(cols[0]);
        c.phase = cols[1];

        c.gains.Kp[0] = safe_stod(cols[2]);
        c.gains.Kp[1] = safe_stod(cols[3]);
        c.gains.Kp[2] = safe_stod(cols[4]);
        c.gains.Kd[0] = safe_stod(cols[5]);
        c.gains.Kd[1] = safe_stod(cols[6]);
        c.gains.Kd[2] = safe_stod(cols[7]);

        c.cost = safe_stod(cols[8]);
        c.final_attitude_error = safe_stod(cols[9]);
        c.max_angular_rate = safe_stod(cols[10]);
        c.total_power = safe_stod(cols[11]);
        c.stable = (std::stoi(cols[12]) != 0);

        if (c.stable && std::isfinite(c.cost)) {
            candidates.push_back(c);
        }
    }

    return candidates;
}

bool compare_candidate_cost(const MonteCarloCandidate& a, const MonteCarloCandidate& b)
{
    return a.cost < b.cost;
}

std::vector<InitialSeed> build_initial_seeds(const std::string& montecarlo_path)
{
    std::vector<InitialSeed> seeds;
    std::vector<MonteCarloCandidate> all_candidates = read_montecarlo_results(montecarlo_path);

    std::vector<MonteCarloCandidate> exploration_candidates;
    std::vector<MonteCarloCandidate> exploitation_candidates;

    for (const auto& c : all_candidates) {
        if (c.phase == "exploration") {
            exploration_candidates.push_back(c);
        }
        else {
            // In the current Monte Carlo code the local refinement phases are named
            // local_from_best_1, local_from_best_2, etc. Any non-exploration phase
            // is treated as an exploitation/local-refinement result.
            exploitation_candidates.push_back(c);
        }
    }

    std::sort(exploration_candidates.begin(), exploration_candidates.end(), compare_candidate_cost);
    std::sort(exploitation_candidates.begin(), exploitation_candidates.end(), compare_candidate_cost);
    std::sort(all_candidates.begin(), all_candidates.end(), compare_candidate_cost);

    if (static_cast<int>(exploration_candidates.size()) < N_EXPLORATION_SEEDS) {
        std::cerr << "Error: not enough stable exploration candidates in "
                  << montecarlo_path << ". Found " << exploration_candidates.size()
                  << ", required " << N_EXPLORATION_SEEDS << "." << std::endl;
        return seeds;
    }

    // 18 best exploration particles.
    for (int i = 0; i < N_EXPLORATION_SEEDS; i++) {
        InitialSeed s;
        s.source = "MC_exploration_best_" + std::to_string(i + 1);
        s.simulation_id = exploration_candidates[i].simulation_id;
        s.csv_cost = exploration_candidates[i].cost;
        s.gains = exploration_candidates[i].gains;
        seeds.push_back(s);
    }

    // 1 best final/local exploitation particle.
    // If no local phase exists, fall back to the best overall stable candidate.
    MonteCarloCandidate final_candidate;
    if (!exploitation_candidates.empty()) {
        final_candidate = exploitation_candidates.front();
    }
    else if (!all_candidates.empty()) {
        final_candidate = all_candidates.front();
    }
    else {
        std::cerr << "Error: no stable Monte Carlo candidates found in "
                  << montecarlo_path << "." << std::endl;
        return std::vector<InitialSeed>();
    }

    InitialSeed final_seed;
    final_seed.source = "MC_best_final_exploitation";
    final_seed.simulation_id = final_candidate.simulation_id;
    final_seed.csv_cost = final_candidate.cost;
    final_seed.gains = final_candidate.gains;
    seeds.push_back(final_seed);

    // 1 particle from Pau Climent's TFG constants.
    InitialSeed pau_seed;
    pau_seed.source = "Pau_TFG_reference";
    pau_seed.simulation_id = -1;
    pau_seed.csv_cost = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < 3; i++) {
        pau_seed.gains.Kp[i] = PAU_TFG_KP[i];
        pau_seed.gains.Kd[i] = PAU_TFG_KD[i];
    }
    seeds.push_back(pau_seed);

    if (static_cast<int>(seeds.size()) != N_PARTICLES) {
        std::cerr << "Error: initial seed count is " << seeds.size()
                  << " but N_PARTICLES is " << N_PARTICLES << "." << std::endl;
        return std::vector<InitialSeed>();
    }

    return seeds;
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
    (void)d;
    return VMAX_LOG;
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

void gains_to_position(const Gains& gains, double position[DIM])
{
    for (int i = 0; i < 3; i++) {
        const double kp = clamp_value(gains.Kp[i], KP_MIN, KP_MAX);
        const double kd = clamp_value(gains.Kd[i], KD_MIN, KD_MAX);

        position[i] = std::log10(kp);
        position[i + 3] = std::log10(kd);
    }
}

void initialize_particle_from_seed(Particle& p,
                                   const InitialSeed& seed,
                                   std::mt19937& rng)
{
    p.best_cost = std::numeric_limits<double>::infinity();
    p.has_personal_best = false;

    double seed_position[DIM];
    gains_to_position(seed.gains, seed_position);

    for (int d = 0; d < DIM; d++) {
        const double x_min = get_lower_bound(d);
        const double x_max = get_upper_bound(d);
        const double vmax = get_vmax(d);

        // The initial position is exactly the selected Monte Carlo / Pau case.
        p.position[d] = clamp_value(seed_position[d], x_min, x_max);

        // A small initial velocity lets the particle start moving after the first
        // evaluation while preserving the exact seed at iteration 1.
        p.velocity[d] = uniform_real(-INITIAL_VELOCITY_FRACTION * vmax,
                                      INITIAL_VELOCITY_FRACTION * vmax,
                                      rng);

        p.best_position[d] = p.position[d];
    }
}

void write_initial_swarm_header(std::ofstream& file)
{
    file << "Particle,SeedSource,SeedSimulation,SeedCSVCost,"
         << "Kp_x,Kp_y,Kp_z,"
         << "Kd_x,Kd_y,Kd_z\n";
}

void write_initial_swarm_row(std::ofstream& file,
                             int particle_id,
                             const InitialSeed& seed)
{
    file << particle_id << ","
         << seed.source << ","
         << seed.simulation_id << ","
         << seed.csv_cost << ","
         << seed.gains.Kp[0] << "," << seed.gains.Kp[1] << "," << seed.gains.Kp[2] << ","
         << seed.gains.Kd[0] << "," << seed.gains.Kd[1] << "," << seed.gains.Kd[2]
         << "\n";
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

int main()
{
    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<InitialSeed> initial_seeds = build_initial_seeds(MONTECARLO_RESULTS_PATH);
    if (static_cast<int>(initial_seeds.size()) != N_PARTICLES) {
        std::cerr << "Could not initialize the requested 20-particle swarm." << std::endl;
        return 1;
    }

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

    std::ofstream initial_swarm_file("../data/pso_initial_swarm.csv");
    if (!initial_swarm_file.is_open()) {
        std::cerr << "Error opening ../data/pso_initial_swarm.csv" << std::endl;
        return 1;
    }

    write_result_header(results_file);
    write_best_history_header(best_history_file);
    write_initial_swarm_header(initial_swarm_file);

    Particle swarm[N_PARTICLES];

    std::cout << "PSO configuration: " << N_PARTICLES
              << " particles x " << N_ITERATIONS
              << " iterations = " << N_PARTICLES * N_ITERATIONS
              << " simulations." << std::endl;

    std::cout << "Initial swarm composition:" << std::endl;
    std::cout << "  - 18 particles: best stable exploration cases from Monte Carlo CSV" << std::endl;
    std::cout << "  -  1 particle : best stable exploitation/local case from Monte Carlo CSV" << std::endl;
    std::cout << "  -  1 particle : Pau Climent TFG reference constants" << std::endl;

    for (int i = 0; i < N_PARTICLES; i++) {
        initialize_particle_from_seed(swarm[i], initial_seeds[i], rng);
        write_initial_swarm_row(initial_swarm_file, i + 1, initial_seeds[i]);

        std::cout << "Particle " << i + 1
                  << " initialized from " << initial_seeds[i].source;
        if (initial_seeds[i].simulation_id > 0) {
            std::cout << " | MC simulation " << initial_seeds[i].simulation_id
                      << " | CSV cost " << initial_seeds[i].csv_cost;
        }
        std::cout << std::endl;
    }

    initial_swarm_file.close();

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
    for (int iter = 0; iter < N_ITERATIONS; iter++) {

        double w = W_MAX;
        if (N_ITERATIONS > 1) {
            w = W_MAX - (W_MAX - W_MIN) * double(iter) / double(N_ITERATIONS - 1);
        }

        std::cout << "\n===== PSO ITERATION "
                  << iter + 1 << "/" << N_ITERATIONS
                  << " | w = " << w
                  << " =====" << std::endl;

        // ------------------------------------------------------------
        // Evaluate all particles
        // ------------------------------------------------------------
        for (int i = 0; i < N_PARTICLES; i++) {

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

            std::cout << "[Particle " << i + 1 << "/" << N_PARTICLES << "] "
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
        for (int i = 0; i < N_PARTICLES; i++) {

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

    return 0;
}
