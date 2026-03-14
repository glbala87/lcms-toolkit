#pragma once

#include "../spectrum.hpp"
#include "../peak.hpp"
#include "../types.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace lcms {
namespace algorithms {

/**
 * @brief Represents a detected isotope envelope.
 */
struct IsotopePattern {
    MZ monoisotopic_mz = 0.0;
    ChargeState charge = 1;
    std::vector<std::pair<MZ, Intensity>> peaks;
    double score = 0.0;
    RetentionTime rt = 0.0;

    /// Calculate neutral mass
    [[nodiscard]] double neutralMass() const {
        constexpr double PROTON_MASS = 1.007276;
        return (monoisotopic_mz - PROTON_MASS) * std::abs(static_cast<int>(charge));
    }

    /// Total intensity across all peaks
    [[nodiscard]] Intensity totalIntensity() const {
        Intensity total = 0.0;
        for (const auto& [mz, intensity] : peaks) {
            total += intensity;
        }
        return total;
    }

    /// Number of isotope peaks
    [[nodiscard]] size_t numPeaks() const { return peaks.size(); }
};

/**
 * @brief Represents a deconvoluted neutral mass.
 */
struct DeconvolutedMass {
    double neutral_mass = 0.0;
    Intensity intensity = 0.0;
    std::vector<ChargeState> charge_states;
    std::vector<IsotopePattern> patterns;
    double quality_score = 0.0;

    [[nodiscard]] size_t numChargeStates() const { return charge_states.size(); }
};

/**
 * @brief Options for isotope pattern detection.
 */
struct IsotopeDetectionOptions {
    /// Minimum charge state
    int min_charge = 1;

    /// Maximum charge state
    int max_charge = 6;

    /// m/z tolerance for matching isotope peaks (Da)
    double mz_tolerance = 0.01;

    /// Minimum number of isotope peaks
    int min_peaks = 2;

    /// Minimum pattern score (cosine similarity)
    double min_score = 0.5;

    /// Minimum intensity for seed peaks
    Intensity min_intensity = 0.0;

    /// Tolerance for grouping neutral masses (Da)
    double mass_tolerance = 0.5;
};

/// Neutron mass difference (C13 - C12)
constexpr double NEUTRON_MASS = 1.003355;

/// Proton mass
constexpr double PROTON_MASS = 1.007276;

/// Averagine residue mass
constexpr double AVERAGINE_MASS = 111.1254;

/**
 * @brief Generate theoretical isotope distribution using averagine model.
 *
 * @param mass Neutral mass
 * @param num_peaks Number of peaks to generate
 * @return Normalized intensity distribution
 */
std::vector<double> averagineDistribution(double mass, int num_peaks = 6);

/**
 * @brief Assign charge state from m/z spacing.
 *
 * @param mz_peaks Sorted m/z values
 * @param mz_tolerance Tolerance for matching
 * @param max_charge Maximum charge to consider
 * @return Pair of (charge_state, confidence)
 */
std::pair<int, double> assignChargeState(
    const std::vector<MZ>& mz_peaks,
    double mz_tolerance = 0.01,
    int max_charge = 6);

/**
 * @brief Detect isotope patterns in a spectrum.
 *
 * @param spectrum Input spectrum
 * @param options Detection options
 * @return Detected isotope patterns sorted by score
 */
std::vector<IsotopePattern> detectIsotopePatterns(
    const Spectrum& spectrum,
    const IsotopeDetectionOptions& options = {});

/**
 * @brief Deconvolute spectrum to neutral masses.
 *
 * @param spectrum Input spectrum
 * @param options Detection options
 * @return Deconvoluted masses sorted by intensity
 */
std::vector<DeconvolutedMass> deconvoluteSpectrum(
    const Spectrum& spectrum,
    const IsotopeDetectionOptions& options = {});

} // namespace algorithms
} // namespace lcms
