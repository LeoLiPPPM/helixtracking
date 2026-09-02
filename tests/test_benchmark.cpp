#include "helixtracking/benchmark.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

namespace {

void require(const bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}

bool nearly_equal(
    const double a,
    const double b,
    const double tolerance = 1e-12
) {
    return std::abs(a - b) < tolerance;
}

std::size_t count_csv_columns(const std::string& line) {
    if (line.empty()) {
        return 0;
    }

    std::size_t columns = 1;

    for (const char c : line) {
        if (c == ',') {
            ++columns;
        }
    }

    return columns;
}

}  // namespace


int main() {

    // ------------------------------------------------------------
    // 1. Baseline event should reconstruct successfully.
    // ------------------------------------------------------------

    const helixtracking::SimulationConfig simulation{};
    const helixtracking::ReconstructionConfig reconstruction{};

    const auto baseline =
        helixtracking::run_benchmark_event(
            "baseline_test",
            0,
            20260902,
            simulation,
            reconstruction
        );

    require(
        baseline.success,
        "baseline benchmark reconstruction succeeds"
    );


    // ------------------------------------------------------------
    // 2. Configuration information must be copied correctly.
    // ------------------------------------------------------------

    require(
        baseline.experiment == "baseline_test",
        "experiment label preserved"
    );

    require(
        baseline.trial == 0,
        "trial number preserved"
    );

    require(
        baseline.seed == 20260902,
        "seed preserved"
    );

    require(
        baseline.measurement_count
            == simulation.measurement_count,
        "measurement count preserved"
    );

    require(
        baseline.outlier_count
            == simulation.outlier_count,
        "outlier count preserved"
    );

    require(
        nearly_equal(
            baseline.sigma_xy_m,
            simulation.transverse_sigma_m
        ),
        "transverse resolution preserved"
    );

    require(
        nearly_equal(
            baseline.sigma_z_m,
            simulation.longitudinal_sigma_m
        ),
        "longitudinal resolution preserved"
    );

    require(
        nearly_equal(
            baseline.magnetic_field_t,
            reconstruction.magnetic_field_t
        ),
        "magnetic field preserved"
    );

    require(
        nearly_equal(
            baseline.true_radius_m,
            simulation.radius_m
        ),
        "true radius preserved"
    );

    require(
        baseline.ransac_iterations
            == reconstruction.ransac_iterations,
        "RANSAC iteration count preserved"
    );

    require(
        nearly_equal(
            baseline.radial_threshold_m,
            reconstruction.radial_threshold_m
        ),
        "radial threshold preserved"
    );


    // ------------------------------------------------------------
    // 3. True pT bookkeeping must obey pT = 0.299792458 |q| B R.
    // ------------------------------------------------------------

    constexpr double kMomentumFactor = 0.299792458;

    const double expected_true_pt =
        kMomentumFactor
        * std::abs(
            static_cast<double>(
                reconstruction.charge_number
            )
        )
        * reconstruction.magnetic_field_t
        * simulation.radius_m;

    require(
        nearly_equal(
            baseline.true_pt_gev_c,
            expected_true_pt
        ),
        "true transverse momentum calculation"
    );


    // ------------------------------------------------------------
    // 4. Successful-event classification metrics must be valid.
    // ------------------------------------------------------------

    require(
        baseline.efficiency >= 0.0
            && baseline.efficiency <= 1.0,
        "efficiency lies in [0, 1]"
    );

    require(
        baseline.false_positive_rate >= 0.0
            && baseline.false_positive_rate <= 1.0,
        "false-positive rate lies in [0, 1]"
    );

    require(
        baseline.purity >= 0.0
            && baseline.purity <= 1.0,
        "purity lies in [0, 1]"
    );


    // ------------------------------------------------------------
    // 5. Reconstruction error quantities must be physically valid.
    // ------------------------------------------------------------

    require(
        baseline.pt_relative_error >= 0.0,
        "absolute pT relative error is nonnegative"
    );

    require(
        baseline.radial_rmse_m >= 0.0,
        "radial RMSE is nonnegative"
    );

    require(
        baseline.longitudinal_rmse_m >= 0.0,
        "longitudinal RMSE is nonnegative"
    );

    require(
        baseline.runtime_ms >= 0.0,
        "runtime is nonnegative"
    );


    // ------------------------------------------------------------
    // 6. Signed and absolute pT error must be consistent.
    // ------------------------------------------------------------

    require(
        nearly_equal(
            baseline.pt_relative_error,
            std::abs(
                baseline.pt_signed_relative_error
            )
        ),
        "absolute pT error equals absolute signed pT error"
    );


    // ------------------------------------------------------------
    // 7. Known baseline should remain reasonably high quality.
    //
    // These intentionally use broad physical bounds rather than
    // demanding one exact floating-point reconstruction.
    // ------------------------------------------------------------

    require(
        baseline.efficiency > 0.95,
        "baseline tracking efficiency"
    );

    require(
        baseline.false_positive_rate < 0.10,
        "baseline false-positive rejection"
    );

    require(
        baseline.purity > 0.95,
        "baseline track purity"
    );

    require(
        baseline.pt_relative_error < 0.01,
        "baseline pT reconstruction within 1 percent"
    );

    require(
        baseline.radial_rmse_m < 0.005,
        "baseline radial RMSE below 5 mm"
    );


    // ------------------------------------------------------------
    // 8. Same seed + same configuration must reproduce the same
    //    physics result.
    //
    // runtime_ms is intentionally NOT compared.
    // ------------------------------------------------------------

    const auto repeated =
        helixtracking::run_benchmark_event(
            "baseline_test",
            0,
            20260902,
            simulation,
            reconstruction
        );

    require(
        repeated.success == baseline.success,
        "deterministic success state"
    );

    require(
        nearly_equal(
            repeated.efficiency,
            baseline.efficiency
        ),
        "deterministic efficiency"
    );

    require(
        nearly_equal(
            repeated.false_positive_rate,
            baseline.false_positive_rate
        ),
        "deterministic false-positive rate"
    );

    require(
        nearly_equal(
            repeated.purity,
            baseline.purity
        ),
        "deterministic purity"
    );

    require(
        nearly_equal(
            repeated.pt_relative_error,
            baseline.pt_relative_error
        ),
        "deterministic pT error"
    );

    require(
        nearly_equal(
            repeated.pt_signed_relative_error,
            baseline.pt_signed_relative_error
        ),
        "deterministic signed pT error"
    );

    require(
        nearly_equal(
            repeated.radial_rmse_m,
            baseline.radial_rmse_m
        ),
        "deterministic radial RMSE"
    );

    require(
        nearly_equal(
            repeated.longitudinal_rmse_m,
            baseline.longitudinal_rmse_m
        ),
        "deterministic longitudinal RMSE"
    );


    // ------------------------------------------------------------
    // 9. Changing seed should actually produce a different event.
    // ------------------------------------------------------------

    const auto different_seed =
        helixtracking::run_benchmark_event(
            "different_seed",
            1,
            20260903,
            simulation,
            reconstruction
        );

    require(
        different_seed.success,
        "second deterministic seed reconstructs successfully"
    );

    const bool physics_result_changed =
        !nearly_equal(
            different_seed.radial_rmse_m,
            baseline.radial_rmse_m
        )
        ||
        !nearly_equal(
            different_seed.pt_relative_error,
            baseline.pt_relative_error
        );

    require(
        physics_result_changed,
        "different seed changes simulated event result"
    );


    // ------------------------------------------------------------
    // 10. Zero-outlier configuration should report zero FPR.
    // ------------------------------------------------------------

    auto no_outlier_simulation = simulation;
    no_outlier_simulation.outlier_count = 0;

    const auto no_outliers =
        helixtracking::run_benchmark_event(
            "no_outliers",
            0,
            424242,
            no_outlier_simulation,
            reconstruction
        );

    require(
        no_outliers.success,
        "zero-outlier event reconstructs successfully"
    );

    require(
        nearly_equal(
            no_outliers.false_positive_rate,
            0.0
        ),
        "zero injected outliers gives zero false-positive rate"
    );


    // ------------------------------------------------------------
    // 11. CSV schema should contain the expected number of columns.
    // ------------------------------------------------------------

    const std::string header =
        helixtracking::benchmark_csv_header();

    const std::string row =
        helixtracking::benchmark_csv_row(
            baseline
        );

    constexpr std::size_t expected_columns = 21;

    require(
        count_csv_columns(header)
            == expected_columns,
        "CSV header has 21 columns"
    );

    require(
        count_csv_columns(row)
            == expected_columns,
        "CSV row has 21 columns"
    );


    // ------------------------------------------------------------
    // 12. Important CSV fields should actually be present.
    // ------------------------------------------------------------

    require(
        header.find("experiment")
            != std::string::npos,
        "CSV contains experiment column"
    );

    require(
        header.find("seed")
            != std::string::npos,
        "CSV contains seed column"
    );

    require(
        header.find("efficiency")
            != std::string::npos,
        "CSV contains efficiency column"
    );

    require(
        header.find("false_positive_rate")
            != std::string::npos,
        "CSV contains false-positive-rate column"
    );

    require(
        header.find("pt_relative_error")
            != std::string::npos,
        "CSV contains pT error column"
    );

    require(
        header.find("runtime_ms")
            != std::string::npos,
        "CSV contains runtime column"
    );


    std::cout
        << "All HelixTracking benchmark tests passed.\n";

    return 0;
}