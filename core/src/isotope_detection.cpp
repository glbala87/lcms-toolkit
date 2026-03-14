#include "lcms/algorithms/isotope_detection.hpp"
#include <numeric>
#include <cmath>
#include <algorithm>

namespace lcms {
namespace algorithms {

std::vector<double> averagineDistribution(double mass, int num_peaks) {
    std::vector<double> distribution(num_peaks, 0.0);
    if (num_peaks <= 0) return distribution;

    if (mass <= 0) {
        distribution[0] = 1.0;
        return distribution;
    }

    // Number of averagine units
    double n_units = mass / AVERAGINE_MASS;

    // Expected heavy isotope count (Poisson lambda)
    // C13: 1.1%, N15: 0.366%, O18: 0.205%, S34: 4.25%, H2: 0.012%
    double lambda = n_units * (
        4.9384 * 0.0111 +   // C13
        1.3577 * 0.00366 +  // N15
        1.4773 * 0.00205 +  // O18
        0.0417 * 0.0425 +   // S34
        7.7583 * 0.00012    // H2
    );

    // Poisson distribution
    for (int i = 0; i < num_peaks; ++i) {
        double log_prob = -lambda + i * std::log(lambda);
        // log(i!)
        double log_fact = 0.0;
        for (int j = 2; j <= i; ++j) {
            log_fact += std::log(static_cast<double>(j));
        }
        distribution[i] = std::exp(log_prob - log_fact);
    }

    // Normalize to max = 1
    double max_val = *std::max_element(distribution.begin(), distribution.end());
    if (max_val > 0) {
        for (auto& v : distribution) {
            v /= max_val;
        }
    }

    return distribution;
}

std::pair<int, double> assignChargeState(
    const std::vector<MZ>& mz_peaks,
    double mz_tolerance,
    int max_charge) {

    if (mz_peaks.size() < 2) return {0, 0.0};

    std::vector<MZ> sorted_mz = mz_peaks;
    std::sort(sorted_mz.begin(), sorted_mz.end());

    std::vector<double> spacings;
    for (size_t i = 1; i < sorted_mz.size(); ++i) {
        spacings.push_back(sorted_mz[i] - sorted_mz[i - 1]);
    }

    int best_charge = 0;
    double best_confidence = 0.0;

    for (int charge = 1; charge <= max_charge; ++charge) {
        double expected = NEUTRON_MASS / charge;
        int matches = 0;
        double total_error = 0.0;

        for (double spacing : spacings) {
            double error = std::abs(spacing - expected);
            if (error <= mz_tolerance) {
                matches++;
                total_error += error;
            }
        }

        if (matches > 0) {
            double confidence = static_cast<double>(matches) / spacings.size();
            double avg_error = total_error / matches;
            confidence *= (1.0 - avg_error / mz_tolerance);
            if (confidence > best_confidence) {
                best_confidence = confidence;
                best_charge = charge;
            }
        }
    }

    return {best_charge, best_confidence};
}

namespace {

double cosineSimilarity(const std::vector<double>& a, const std::vector<double>& b) {
    size_t n = std::min(a.size(), b.size());
    if (n == 0) return 0.0;

    double dot = 0.0, norm_a = 0.0, norm_b = 0.0;
    for (size_t i = 0; i < n; ++i) {
        dot += a[i] * b[i];
        norm_a += a[i] * a[i];
        norm_b += b[i] * b[i];
    }

    norm_a = std::sqrt(norm_a);
    norm_b = std::sqrt(norm_b);

    if (norm_a == 0.0 || norm_b == 0.0) return 0.0;
    return dot / (norm_a * norm_b);
}

} // anonymous namespace

std::vector<IsotopePattern> detectIsotopePatterns(
    const Spectrum& spectrum,
    const IsotopeDetectionOptions& options) {

    size_t n = spectrum.size();
    if (n < static_cast<size_t>(options.min_peaks)) return {};

    const auto& mz_data = spectrum.mzData();
    const auto& int_data = spectrum.intensityData();

    // Create sorted indices by m/z
    std::vector<Index> mz_order(n);
    std::iota(mz_order.begin(), mz_order.end(), 0);
    std::sort(mz_order.begin(), mz_order.end(),
              [&](Index a, Index b) { return mz_data[a] < mz_data[b]; });

    // Create intensity-sorted indices for seed selection
    std::vector<Index> int_order(n);
    std::iota(int_order.begin(), int_order.end(), 0);
    std::sort(int_order.begin(), int_order.end(),
              [&](Index a, Index b) { return int_data[a] > int_data[b]; });

    std::vector<bool> used(n, false);
    std::vector<IsotopePattern> patterns;

    for (Index seed_idx : int_order) {
        if (used[seed_idx]) continue;
        if (int_data[seed_idx] < options.min_intensity) continue;

        MZ seed_mz = mz_data[seed_idx];
        Intensity seed_int = int_data[seed_idx];

        IsotopePattern best_pattern;
        double best_score = -1.0;
        std::vector<Index> best_indices;

        for (int charge = options.min_charge; charge <= options.max_charge; ++charge) {
            double expected_spacing = NEUTRON_MASS / charge;

            std::vector<std::pair<MZ, Intensity>> pattern_peaks = {{seed_mz, seed_int}};
            std::vector<Index> pattern_indices = {seed_idx};

            // Forward isotope peaks
            for (int iso = 1; iso <= 10; ++iso) {
                MZ expected_mz = seed_mz + iso * expected_spacing;
                Index best_cand = n;
                double best_diff = options.mz_tolerance + 1.0;

                for (Index idx : mz_order) {
                    if (used[idx]) continue;
                    double diff = std::abs(mz_data[idx] - expected_mz);
                    if (diff <= options.mz_tolerance && diff < best_diff) {
                        best_diff = diff;
                        best_cand = idx;
                    }
                }

                if (best_cand >= n) break;
                pattern_peaks.push_back({mz_data[best_cand], int_data[best_cand]});
                pattern_indices.push_back(best_cand);
            }

            if (static_cast<int>(pattern_peaks.size()) < options.min_peaks) continue;

            // Score against theoretical
            MZ mono_mz = pattern_peaks[0].first;
            double neutral_mass = (mono_mz - PROTON_MASS) * charge;
            auto theoretical = averagineDistribution(neutral_mass,
                                                      static_cast<int>(pattern_peaks.size()));

            std::vector<double> observed;
            double max_obs = 0.0;
            for (const auto& [mz, intensity] : pattern_peaks) {
                observed.push_back(intensity);
                if (intensity > max_obs) max_obs = intensity;
            }
            if (max_obs > 0) {
                for (auto& v : observed) v /= max_obs;
            }

            double score = cosineSimilarity(theoretical, observed);

            if (score > best_score && score >= options.min_score) {
                best_score = score;
                best_pattern.monoisotopic_mz = mono_mz;
                best_pattern.charge = static_cast<ChargeState>(charge);
                best_pattern.peaks = pattern_peaks;
                best_pattern.score = score;
                best_pattern.rt = spectrum.rt();
                best_indices = pattern_indices;
            }
        }

        if (best_score >= options.min_score) {
            patterns.push_back(best_pattern);
            for (Index idx : best_indices) {
                used[idx] = true;
            }
        }
    }

    std::sort(patterns.begin(), patterns.end(),
              [](const IsotopePattern& a, const IsotopePattern& b) {
                  return a.score > b.score;
              });

    return patterns;
}

std::vector<DeconvolutedMass> deconvoluteSpectrum(
    const Spectrum& spectrum,
    const IsotopeDetectionOptions& options) {

    auto patterns = detectIsotopePatterns(spectrum, options);
    if (patterns.empty()) return {};

    std::vector<DeconvolutedMass> masses;

    for (const auto& pattern : patterns) {
        double neutral_mass = pattern.neutralMass();
        bool merged = false;

        for (auto& dm : masses) {
            if (std::abs(dm.neutral_mass - neutral_mass) <= options.mass_tolerance) {
                dm.intensity += pattern.totalIntensity();
                int charge = static_cast<int>(pattern.charge);
                bool found = false;
                for (auto cs : dm.charge_states) {
                    if (cs == pattern.charge) { found = true; break; }
                }
                if (!found) dm.charge_states.push_back(pattern.charge);
                dm.patterns.push_back(pattern);

                // Weighted average mass
                double total_int = 0.0;
                double weighted_mass = 0.0;
                for (const auto& p : dm.patterns) {
                    double ti = p.totalIntensity();
                    total_int += ti;
                    weighted_mass += p.neutralMass() * ti;
                }
                if (total_int > 0) dm.neutral_mass = weighted_mass / total_int;

                merged = true;
                break;
            }
        }

        if (!merged) {
            DeconvolutedMass dm;
            dm.neutral_mass = neutral_mass;
            dm.intensity = pattern.totalIntensity();
            dm.charge_states = {pattern.charge};
            dm.patterns = {pattern};
            dm.quality_score = pattern.score;
            masses.push_back(dm);
        }
    }

    // Update quality scores
    for (auto& dm : masses) {
        double avg_score = 0.0;
        for (const auto& p : dm.patterns) avg_score += p.score;
        avg_score /= dm.patterns.size();
        double charge_bonus = std::min(0.2, 0.1 * (dm.numChargeStates() - 1.0));
        dm.quality_score = std::min(1.0, avg_score + charge_bonus);
    }

    std::sort(masses.begin(), masses.end(),
              [](const DeconvolutedMass& a, const DeconvolutedMass& b) {
                  return a.intensity > b.intensity;
              });

    return masses;
}

} // namespace algorithms
} // namespace lcms
