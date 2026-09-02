#include "helixtracking/benchmark.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace {

constexpr std::uint64_t kMasterSeed = 20260902;
constexpr std::size_t kTrialsPerCondition = 500;

}

int main() {

    std::filesystem::create_directories(
        "outputs/benchmark"
    );

    std::ofstream csv(
        "outputs/benchmark/monte_carlo.csv"
    );

    csv
        << helixtracking::benchmark_csv_header()
        << '\n';

    const helixtracking::SimulationConfig baseline_sim{};
    const helixtracking::ReconstructionConfig baseline_reco{};

    const std::vector<double> sigma_xy_values{
        0.0005,
        0.0010,
        0.0015,
        0.0020,
        0.0030,
        0.0040,
        0.0060,
        0.0080,
        0.0120
    };

    for (const double sigma_xy : sigma_xy_values) {

        helixtracking::SimulationConfig sim =
            baseline_sim;

        helixtracking::ReconstructionConfig reco =
            baseline_reco;

        sim.transverse_sigma_m = sigma_xy;

        for (
            std::size_t trial = 0;
            trial < kTrialsPerCondition;
            ++trial
        ) {

            const std::uint64_t seed =
                kMasterSeed + trial;

            const auto row =
                helixtracking::run_benchmark_event(
                    "sigma_xy",
                    trial,
                    seed,
                    sim,
                    reco
                );

            csv
                << helixtracking::benchmark_csv_row(row)
                << '\n';
        }

        std::cout
            << "Completed sigma_xy = "
            << 1000.0 * sigma_xy
            << " mm\n";
    }

    std::cout
        << "Monte Carlo benchmark complete.\n";

    return 0;
}