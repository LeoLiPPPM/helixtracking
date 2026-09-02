#include "helixtracking/track.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>

namespace helixtracking {
namespace {

constexpr double kMomentumFactor = 0.299792458;  // pT[GeV/c] = factor * q * B[T] * R[m]
constexpr double kTwoPi = 2.0 * std::numbers::pi;

double radial_residual(const Hit& hit, const Circle& circle) {
    return std::abs(std::hypot(hit.x_m - circle.center_x_m, hit.y_m - circle.center_y_m) - circle.radius_m);
}

std::array<double, 3> solve_three_by_three(
    std::array<std::array<double, 3>, 3> matrix,
    std::array<double, 3> rhs
) {
    for (std::size_t pivot = 0; pivot < 3; ++pivot) {
        std::size_t best = pivot;
        for (std::size_t row = pivot + 1; row < 3; ++row) {
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot])) {
                best = row;
            }
        }
        if (std::abs(matrix[best][pivot]) < 1e-14) {
            throw std::runtime_error("singular normal equation in circle refinement");
        }
        std::swap(matrix[pivot], matrix[best]);
        std::swap(rhs[pivot], rhs[best]);
        for (std::size_t row = pivot + 1; row < 3; ++row) {
            const double factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (std::size_t col = pivot; col < 3; ++col) {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    std::array<double, 3> answer{};
    for (int row = 2; row >= 0; --row) {
        double value = rhs[static_cast<std::size_t>(row)];
        for (std::size_t col = static_cast<std::size_t>(row) + 1; col < 3; ++col) {
            value -= matrix[static_cast<std::size_t>(row)][col] * answer[col];
        }
        answer[static_cast<std::size_t>(row)] = value / matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }
    return answer;
}

Circle refine_circle(const std::vector<Hit>& hits, const std::vector<std::size_t>& indices) {
    if (indices.size() < 3) {
        throw std::runtime_error("at least three inliers are required to refine a circle");
    }
    std::array<std::array<double, 3>, 3> normal{};
    std::array<double, 3> rhs{};
    for (const std::size_t index : indices) {
        const Hit& hit = hits[index];
        const std::array<double, 3> row{hit.x_m, hit.y_m, 1.0};
        const double target = -(hit.x_m * hit.x_m + hit.y_m * hit.y_m);
        for (std::size_t i = 0; i < 3; ++i) {
            rhs[i] += row[i] * target;
            for (std::size_t j = 0; j < 3; ++j) {
                normal[i][j] += row[i] * row[j];
            }
        }
    }
    const auto [d, e, f] = solve_three_by_three(normal, rhs);
    const double center_x = -0.5 * d;
    const double center_y = -0.5 * e;
    const double radius_squared = center_x * center_x + center_y * center_y - f;
    if (radius_squared <= 0.0) {
        throw std::runtime_error("circle refinement produced a non-positive radius");
    }
    return Circle{center_x, center_y, std::sqrt(radius_squared)};
}

std::vector<double> unwrapped_angles(
    const std::vector<Hit>& hits,
    const std::vector<std::size_t>& ordered_indices,
    const Circle& circle
) {
    std::vector<double> angles;
    angles.reserve(ordered_indices.size());
    for (const std::size_t index : ordered_indices) {
        double angle = std::atan2(hits[index].y_m - circle.center_y_m, hits[index].x_m - circle.center_x_m);
        if (!angles.empty()) {
            while (angle - angles.back() > std::numbers::pi) {
                angle -= kTwoPi;
            }
            while (angle - angles.back() < -std::numbers::pi) {
                angle += kTwoPi;
            }
        }
        angles.push_back(angle);
    }
    return angles;
}

std::pair<double, double> fit_line(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::runtime_error("line fit requires paired samples");
    }
    const double x_mean = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
    const double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / static_cast<double>(y.size());
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        numerator += (x[i] - x_mean) * (y[i] - y_mean);
        denominator += (x[i] - x_mean) * (x[i] - x_mean);
    }
    if (denominator <= std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error("line fit has no angular extent");
    }
    const double slope = numerator / denominator;
    return {slope, y_mean - slope * x_mean};
}

std::vector<std::size_t> sorted_by_layer(const std::vector<Hit>& hits, const std::vector<std::size_t>& indices) {
    std::vector<std::size_t> sorted = indices;
    std::sort(sorted.begin(), sorted.end(), [&hits](const std::size_t left, const std::size_t right) {
        return hits[left].layer < hits[right].layer;
    });
    return sorted;
}

}  // namespace

std::vector<Hit> simulate_helix(const SimulationConfig& config, const std::uint64_t seed) {
    if (config.measurement_count < 4 || config.radius_m <= 0.0 || config.angle_span_rad == 0.0) {
        throw std::invalid_argument("invalid helix simulation configuration");
    }
    std::mt19937_64 generator(seed);
    std::normal_distribution<double> transverse_noise(0.0, config.transverse_sigma_m);
    std::normal_distribution<double> longitudinal_noise(0.0, config.longitudinal_sigma_m);
    std::uniform_real_distribution<double> outlier_xy(-1.25, 1.25);
    std::uniform_real_distribution<double> outlier_z(-0.35, 1.05);
    std::uniform_int_distribution<int> outlier_layer(0, static_cast<int>(config.measurement_count - 1));

    std::vector<Hit> hits;
    hits.reserve(config.measurement_count + config.outlier_count);
    for (std::size_t i = 0; i < config.measurement_count; ++i) {
        const double fraction = static_cast<double>(i) / static_cast<double>(config.measurement_count - 1);
        const double angle = config.start_angle_rad + fraction * config.angle_span_rad;
        hits.push_back(Hit{
            static_cast<int>(i),
            config.center_x_m + config.radius_m * std::cos(angle) + transverse_noise(generator),
            config.center_y_m + config.radius_m * std::sin(angle) + transverse_noise(generator),
            config.pitch_m_per_rad * angle + longitudinal_noise(generator),
            false,
        });
    }
    for (std::size_t i = 0; i < config.outlier_count; ++i) {
        hits.push_back(Hit{
            outlier_layer(generator),
            outlier_xy(generator),
            outlier_xy(generator),
            outlier_z(generator),
            true,
        });
    }
    return hits;
}

Circle circle_through_three(const Hit& a, const Hit& b, const Hit& c) {
    const double determinant = 2.0 * (
        a.x_m * (b.y_m - c.y_m) + b.x_m * (c.y_m - a.y_m) + c.x_m * (a.y_m - b.y_m)
    );
    if (std::abs(determinant) < 1e-12) {
        throw std::runtime_error("three collinear hits do not define a circle");
    }
    const double a2 = a.x_m * a.x_m + a.y_m * a.y_m;
    const double b2 = b.x_m * b.x_m + b.y_m * b.y_m;
    const double c2 = c.x_m * c.x_m + c.y_m * c.y_m;
    const double center_x = (a2 * (b.y_m - c.y_m) + b2 * (c.y_m - a.y_m) + c2 * (a.y_m - b.y_m)) / determinant;
    const double center_y = (a2 * (c.x_m - b.x_m) + b2 * (a.x_m - c.x_m) + c2 * (b.x_m - a.x_m)) / determinant;
    return Circle{center_x, center_y, std::hypot(a.x_m - center_x, a.y_m - center_y)};
}

Reconstruction reconstruct_track(
    const std::vector<Hit>& hits,
    const ReconstructionConfig& config,
    const std::uint64_t seed
) {
    if (hits.size() < 4 || config.ransac_iterations == 0 || config.radial_threshold_m <= 0.0) {
        throw std::invalid_argument("invalid reconstruction input");
    }
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<std::size_t> choose(0, hits.size() - 1);
    std::vector<std::size_t> best_indices;
    double best_error = std::numeric_limits<double>::infinity();

    for (std::size_t iteration = 0; iteration < config.ransac_iterations; ++iteration) {
        const std::size_t i = choose(generator);
        const std::size_t j = choose(generator);
        const std::size_t k = choose(generator);
        if (i == j || i == k || j == k) {
            continue;
        }
        Circle candidate{};
        try {
            candidate = circle_through_three(hits[i], hits[j], hits[k]);
        } catch (const std::runtime_error&) {
            continue;
        }
        if (candidate.radius_m < 0.05 || candidate.radius_m > 5.0) {
            continue;
        }
        std::vector<std::size_t> inliers;
        double error = 0.0;
        for (std::size_t index = 0; index < hits.size(); ++index) {
            const double residual = radial_residual(hits[index], candidate);
            if (residual <= config.radial_threshold_m) {
                inliers.push_back(index);
                error += residual * residual;
            }
        }
        if (inliers.size() > best_indices.size() || (inliers.size() == best_indices.size() && error < best_error)) {
            best_indices = std::move(inliers);
            best_error = error;
        }
    }
    if (best_indices.size() < 4) {
        throw std::runtime_error("RANSAC could not find a valid track");
    }

    Circle circle = refine_circle(hits, best_indices);
    auto ordered = sorted_by_layer(hits, best_indices);
    auto angles = unwrapped_angles(hits, ordered, circle);
    std::vector<double> z_values;
    z_values.reserve(ordered.size());
    for (const std::size_t index : ordered) {
        z_values.push_back(hits[index].z_m);
    }
    auto [pitch, intercept] = fit_line(angles, z_values);

    std::vector<std::size_t> final_indices;
    for (std::size_t position = 0; position < ordered.size(); ++position) {
        const std::size_t index = ordered[position];
        const double z_residual = std::abs(hits[index].z_m - (pitch * angles[position] + intercept));
        if (radial_residual(hits[index], circle) <= config.radial_threshold_m &&
            z_residual <= config.longitudinal_threshold_m) {
            final_indices.push_back(index);
        }
    }
    if (final_indices.size() < 4) {
        throw std::runtime_error("longitudinal consistency rejected too many hits");
    }

    circle = refine_circle(hits, final_indices);
    ordered = sorted_by_layer(hits, final_indices);
    angles = unwrapped_angles(hits, ordered, circle);
    z_values.clear();
    for (const std::size_t index : ordered) {
        z_values.push_back(hits[index].z_m);
    }
    std::tie(pitch, intercept) = fit_line(angles, z_values);

    double radial_squared_error = 0.0;
    double longitudinal_squared_error = 0.0;
    std::vector<bool> mask(hits.size(), false);
    for (std::size_t position = 0; position < ordered.size(); ++position) {
        const std::size_t index = ordered[position];
        mask[index] = true;
        const double radial = radial_residual(hits[index], circle);
        const double longitudinal = hits[index].z_m - (pitch * angles[position] + intercept);
        radial_squared_error += radial * radial;
        longitudinal_squared_error += longitudinal * longitudinal;
    }

    return Reconstruction{
        circle,
        pitch,
        intercept,
        kMomentumFactor * std::abs(static_cast<double>(config.charge_number)) * config.magnetic_field_t * circle.radius_m,
        std::sqrt(radial_squared_error / static_cast<double>(ordered.size())),
        std::sqrt(longitudinal_squared_error / static_cast<double>(ordered.size())),
        std::move(mask),
    };
}

void write_outputs(
    const std::filesystem::path& directory,
    const std::vector<Hit>& hits,
    const Reconstruction& reconstruction,
    const SimulationConfig& truth
) {
    std::filesystem::create_directories(directory);
    std::ofstream hit_file(directory / "hits.csv");
    hit_file << "layer,x_m,y_m,z_m,injected_outlier,reconstructed_inlier\n" << std::setprecision(12);
    for (std::size_t i = 0; i < hits.size(); ++i) {
        const Hit& hit = hits[i];
        hit_file << hit.layer << ',' << hit.x_m << ',' << hit.y_m << ',' << hit.z_m << ','
                 << static_cast<int>(hit.injected_outlier) << ',' << static_cast<int>(reconstruction.inlier_mask[i]) << '\n';
    }

    std::ofstream curve_file(directory / "fitted_curve.csv");
    curve_file << "angle_rad,x_m,y_m,z_m\n" << std::setprecision(12);
    constexpr std::size_t curve_points = 240;
    for (std::size_t i = 0; i < curve_points; ++i) {
        const double fraction = static_cast<double>(i) / static_cast<double>(curve_points - 1);
        const double angle = truth.start_angle_rad + fraction * truth.angle_span_rad;
        curve_file << angle << ','
                   << reconstruction.circle.center_x_m + reconstruction.circle.radius_m * std::cos(angle) << ','
                   << reconstruction.circle.center_y_m + reconstruction.circle.radius_m * std::sin(angle) << ','
                   << reconstruction.pitch_m_per_rad * angle + reconstruction.z_intercept_m << '\n';
    }
}

}  // namespace helixtracking
