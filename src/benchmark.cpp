#include "helixtracking/benchmark.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace helixtracking {

namespace {

constexpr double kMomentumFactor = 0.299792458;

double nan_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

} // namespace

BenchmarkRow run_benchmark_event(
    const std::string& experiment,
    const std::size_t trial,
    const std::uint64_t seed,
    const SimulationConfig& simulation,
    const ReconstructionConfig& reconstruction
) {
    BenchmarkRow row{};

    row.experiment = experiment;
    row.trial = trial;
    row.seed = seed;

    row.measurement_count = simulation.measurement_count;
    row.outlier_count = simulation.outlier_count;

    row.sigma_xy_m = simulation.transverse_sigma_m;
    row.sigma_z_m = simulation.longitudinal_sigma_m;

    row.magnetic_field_t = reconstruction.magnetic_field_t;
    row.true_radius_m = simulation.radius_m;

    row.true_pt_gev_c =
        kMomentumFactor
        * std::abs(static_cast<double>(reconstruction.charge_number))
        * reconstruction.magnetic_field_t
        * simulation.radius_m;

    row.radial_threshold_m = reconstruction.radial_threshold_m;
    row.ransac_iterations = reconstruction.ransac_iterations;

    const auto hits = simulate_helix(simulation, seed);

    try {
        const auto start = std::chrono::steady_clock::now();

        const Reconstruction result =
            reconstruct_track(
                hits,
                reconstruction,
                seed ^ 0x9E3779B97F4A7C15ULL
            );

        const auto stop = std::chrono::steady_clock::now();

        row.runtime_ms =
            std::chrono::duration<double, std::milli>(
                stop - start
            ).count();

        std::size_t true_hits = 0;
        std::size_t outlier_hits = 0;

        std::size_t true_positives = 0;
        std::size_t false_positives = 0;

        for (std::size_t i = 0; i < hits.size(); ++i) {

            if (hits[i].injected_outlier) {
                ++outlier_hits;

                if (result.inlier_mask[i]) {
                    ++false_positives;
                }

            } else {
                ++true_hits;

                if (result.inlier_mask[i]) {
                    ++true_positives;
                }
            }
        }

        row.efficiency =
            static_cast<double>(true_positives)
            / static_cast<double>(true_hits);

        row.false_positive_rate =
            outlier_hits == 0
            ? 0.0
            : static_cast<double>(false_positives)
                / static_cast<double>(outlier_hits);

        const std::size_t accepted =
            true_positives + false_positives;

        row.purity =
            accepted == 0
            ? nan_value()
            : static_cast<double>(true_positives)
                / static_cast<double>(accepted);

        row.pt_signed_relative_error =
            (
                result.transverse_momentum_gev_c
                - row.true_pt_gev_c
            )
            / row.true_pt_gev_c;

        row.pt_relative_error =
            std::abs(row.pt_signed_relative_error);

        row.radial_rmse_m =
            result.radial_rmse_m;

        row.longitudinal_rmse_m =
            result.longitudinal_rmse_m;

        row.success = true;

    } catch (const std::exception&) {

        row.success = false;

        row.efficiency = nan_value();
        row.false_positive_rate = nan_value();
        row.purity = nan_value();

        row.pt_relative_error = nan_value();
        row.pt_signed_relative_error = nan_value();

        row.radial_rmse_m = nan_value();
        row.longitudinal_rmse_m = nan_value();
        row.runtime_ms = nan_value();
    }

    return row;
}


std::string benchmark_csv_header() {
    return
        "experiment,"
        "trial,"
        "seed,"
        "measurement_count,"
        "outlier_count,"
        "sigma_xy_m,"
        "sigma_z_m,"
        "magnetic_field_t,"
        "true_radius_m,"
        "true_pt_gev_c,"
        "radial_threshold_m,"
        "ransac_iterations,"
        "success,"
        "efficiency,"
        "false_positive_rate,"
        "purity,"
        "pt_relative_error,"
        "pt_signed_relative_error,"
        "radial_rmse_m,"
        "longitudinal_rmse_m,"
        "runtime_ms";
}


std::string benchmark_csv_row(
    const BenchmarkRow& row
) {
    std::ostringstream stream;

    stream << std::setprecision(12)

        << row.experiment << ','
        << row.trial << ','
        << row.seed << ','

        << row.measurement_count << ','
        << row.outlier_count << ','

        << row.sigma_xy_m << ','
        << row.sigma_z_m << ','

        << row.magnetic_field_t << ','
        << row.true_radius_m << ','
        << row.true_pt_gev_c << ','

        << row.radial_threshold_m << ','
        << row.ransac_iterations << ','

        << static_cast<int>(row.success) << ','

        << row.efficiency << ','
        << row.false_positive_rate << ','
        << row.purity << ','

        << row.pt_relative_error << ','
        << row.pt_signed_relative_error << ','

        << row.radial_rmse_m << ','
        << row.longitudinal_rmse_m << ','

        << row.runtime_ms;

    return stream.str();
}

} // namespace helixtracking
