#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace helixtracking {

struct Hit {
    int layer{};
    double x_m{};
    double y_m{};
    double z_m{};
    bool injected_outlier{};
};

struct Circle {
    double center_x_m{};
    double center_y_m{};
    double radius_m{};
};

struct SimulationConfig {
    std::size_t measurement_count{56};
    std::size_t outlier_count{14};
    double radius_m{0.82};
    double center_x_m{0.12};
    double center_y_m{-0.08};
    double start_angle_rad{2.4};
    double angle_span_rad{5.0};
    double pitch_m_per_rad{0.14};
    double transverse_sigma_m{0.0015};
    double longitudinal_sigma_m{0.0025};
};

struct ReconstructionConfig {
    std::size_t ransac_iterations{900};
    double radial_threshold_m{0.006};
    double longitudinal_threshold_m{0.012};
    double magnetic_field_t{1.5};
    int charge_number{1};
};

struct Reconstruction {
    Circle circle{};
    double pitch_m_per_rad{};
    double z_intercept_m{};
    double transverse_momentum_gev_c{};
    double radial_rmse_m{};
    double longitudinal_rmse_m{};
    std::vector<bool> inlier_mask;
};

std::vector<Hit> simulate_helix(const SimulationConfig& config, std::uint64_t seed);
Circle circle_through_three(const Hit& a, const Hit& b, const Hit& c);
Reconstruction reconstruct_track(
    const std::vector<Hit>& hits,
    const ReconstructionConfig& config,
    std::uint64_t seed
);
void write_outputs(
    const std::filesystem::path& directory,
    const std::vector<Hit>& hits,
    const Reconstruction& reconstruction,
    const SimulationConfig& truth
);

}  // namespace helixtracking

