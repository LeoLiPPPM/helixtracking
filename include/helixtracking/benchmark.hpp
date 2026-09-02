#pragma once

#include "helixtracking/track.hpp"

#include <cstddef>
#include <cstdint>
#include <string>

namespace helixtracking {

struct BenchmarkRow {
    std::string experiment;
    std::size_t trial{};
    std::uint64_t seed{};

    std::size_t measurement_count{};
    std::size_t outlier_count{};

    double sigma_xy_m{};
    double sigma_z_m{};

    double magnetic_field_t{};
    double true_radius_m{};
    double true_pt_gev_c{};

    double radial_threshold_m{};
    std::size_t ransac_iterations{};

    bool success{};

    double efficiency{};
    double false_positive_rate{};
    double purity{};

    double pt_relative_error{};
    double pt_signed_relative_error{};

    double radial_rmse_m{};
    double longitudinal_rmse_m{};

    double runtime_ms{};
};

BenchmarkRow run_benchmark_event(
    const std::string& experiment,
    std::size_t trial,
    std::uint64_t seed,
    const SimulationConfig& simulation,
    const ReconstructionConfig& reconstruction
);

std::string benchmark_csv_header();

std::string benchmark_csv_row(const BenchmarkRow& row);

} // namespace helixtracking
