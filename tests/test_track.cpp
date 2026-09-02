#include "helixtracking/track.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

void require(const bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}

}  // namespace

int main() {
    const helixtracking::Hit a{0, 3.0, -1.0, 0.0, false};
    const helixtracking::Hit b{1, 1.0, 1.0, 0.0, false};
    const helixtracking::Hit c{2, -1.0, -1.0, 0.0, false};
    const auto exact = helixtracking::circle_through_three(a, b, c);
    require(std::abs(exact.center_x_m - 1.0) < 1e-12, "exact circle x center");
    require(std::abs(exact.center_y_m + 1.0) < 1e-12, "exact circle y center");
    require(std::abs(exact.radius_m - 2.0) < 1e-12, "exact circle radius");

    const helixtracking::SimulationConfig truth{};
    const helixtracking::ReconstructionConfig config{};
    const auto hits = helixtracking::simulate_helix(truth, 17);
    const auto fit = helixtracking::reconstruct_track(hits, config, 19);
    require(std::abs(fit.circle.radius_m - truth.radius_m) < 0.004, "noisy radius recovery");
    require(std::abs(fit.pitch_m_per_rad - truth.pitch_m_per_rad) < 0.003, "noisy pitch recovery");
    require(fit.radial_rmse_m < 0.003, "radial residual quality");

    std::size_t accepted_outliers = 0;
    for (std::size_t i = 0; i < hits.size(); ++i) {
        if (hits[i].injected_outlier && fit.inlier_mask[i]) {
            ++accepted_outliers;
        }
    }
    require(accepted_outliers <= 1, "RANSAC outlier rejection");
    std::cout << "All HelixTracking tests passed.\n";
    return 0;
}

