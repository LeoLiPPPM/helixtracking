#include "helixtracking/track.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct Arguments {
    std::filesystem::path output{"outputs"};
    std::uint64_t seed{20260902};
};

Arguments parse_arguments(const int argc, char** argv) {
    Arguments args;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--output" && i + 1 < argc) {
            args.output = argv[++i];
        } else if (option == "--seed" && i + 1 < argc) {
            args.seed = std::stoull(argv[++i]);
        } else {
            throw std::invalid_argument("usage: helixtracking [--output PATH] [--seed INTEGER]");
        }
    }
    return args;
}

}  // namespace

int main(const int argc, char** argv) {
    try {
        const Arguments args = parse_arguments(argc, argv);
        const helixtracking::SimulationConfig truth{};
        const helixtracking::ReconstructionConfig reconstruction_config{};
        const auto hits = helixtracking::simulate_helix(truth, args.seed);
        const auto result = helixtracking::reconstruct_track(hits, reconstruction_config, args.seed + 1);
        helixtracking::write_outputs(args.output, hits, result, truth);

        std::size_t accepted = 0;
        std::size_t true_positives = 0;
        std::size_t false_positives = 0;
        for (std::size_t i = 0; i < hits.size(); ++i) {
            if (result.inlier_mask[i]) {
                ++accepted;
                if (hits[i].injected_outlier) {
                    ++false_positives;
                } else {
                    ++true_positives;
                }
            }
        }
        const double radius_error_percent = 100.0 * (result.circle.radius_m - truth.radius_m) / truth.radius_m;
        const double pitch_error_percent = 100.0 * (result.pitch_m_per_rad - truth.pitch_m_per_rad) / truth.pitch_m_per_rad;
        const double truth_momentum = 0.299792458 * reconstruction_config.magnetic_field_t * truth.radius_m;
        const double momentum_error_percent = 100.0 * (result.transverse_momentum_gev_c - truth_momentum) / truth_momentum;

        std::filesystem::create_directories(args.output);
        std::ofstream summary(args.output / "summary.json");
        summary << std::fixed << std::setprecision(8)
                << "{\n"
                << "  \"total_hits\": " << hits.size() << ",\n"
                << "  \"accepted_hits\": " << accepted << ",\n"
                << "  \"true_positive_hits\": " << true_positives << ",\n"
                << "  \"false_positive_hits\": " << false_positives << ",\n"
                << "  \"estimated_radius_m\": " << result.circle.radius_m << ",\n"
                << "  \"radius_error_percent\": " << radius_error_percent << ",\n"
                << "  \"estimated_pitch_m_per_rad\": " << result.pitch_m_per_rad << ",\n"
                << "  \"pitch_error_percent\": " << pitch_error_percent << ",\n"
                << "  \"estimated_pt_gev_c\": " << result.transverse_momentum_gev_c << ",\n"
                << "  \"momentum_error_percent\": " << momentum_error_percent << ",\n"
                << "  \"radial_rmse_mm\": " << 1000.0 * result.radial_rmse_m << ",\n"
                << "  \"longitudinal_rmse_mm\": " << 1000.0 * result.longitudinal_rmse_m << "\n"
                << "}\n";

        std::cout << "HelixTracking reconstructed " << accepted << '/' << hits.size()
                  << " hits; radius error = " << radius_error_percent
                  << "%, pT error = " << momentum_error_percent << "%\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "helixtracking: " << error.what() << '\n';
        return 1;
    }
}

