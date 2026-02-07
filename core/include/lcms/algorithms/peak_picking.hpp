#pragma once

#include "../spectrum.hpp"
#include "../peak.hpp"
#include <vector>
#include <functional>
#include <cmath>

namespace lcms {
namespace algorithms {

/**
 * @brief Peak picking method
 */
enum class PeakPickingMethod : std::uint8_t {
    LOCAL_MAXIMUM,      // Simple local maximum detection
    CENTROID,           // Weighted centroid calculation
    GAUSSIAN_FIT,       // Gaussian curve fitting
    WAVELET             // Continuous wavelet transform
};

/**
 * @brief Options for peak picking
 */
struct PeakPickingOptions {
    /// Picking method
    PeakPickingMethod method = PeakPickingMethod::CENTROID;

    /// Minimum signal-to-noise ratio
    double min_snr = 3.0;

    /// Minimum absolute intensity
    Intensity min_intensity = 0.0;

    /// Minimum peak width (in m/z units)
    MZ min_width = 0.0;

    /// Maximum peak width (in m/z units, 0 = no limit)
    MZ max_width = 0.0;

    /// Number of points for local maximum detection
    int window_size = 3;

    /// Noise estimation method ("mad", "percentile", "fixed")
    std::string noise_method = "mad";

    /// Fixed noise level (when noise_method = "fixed")
    double fixed_noise = 0.0;

    /// Percentile for noise estimation (when noise_method = "percentile")
    double noise_percentile = 5.0;

    /// Fit Gaussian to each peak for better parameters
    bool fit_peaks = true;

    /// Maximum number of iterations for Gaussian fitting
    int max_fit_iterations = 20;
};

/**
 * @brief Peak picking algorithm for mass spectra.
 *
 * Detects peaks in profile mass spectra and converts them to
 * centroided data with peak properties (area, width, SNR).
 */
class PeakPicker {
public:
    PeakPicker() = default;
    explicit PeakPicker(const PeakPickingOptions& options)
        : options_(options) {}

    /**
     * @brief Pick peaks from a spectrum.
     *
     * @param spectrum Input profile spectrum
     * @return Detected peaks
     */
    PeakList pick(const Spectrum& spectrum);

    /**
     * @brief Pick peaks and create centroided spectrum.
     *
     * @param spectrum Input profile spectrum
     * @return Centroided spectrum
     */
    Spectrum centroid(const Spectrum& spectrum);

    /**
     * @brief Estimate noise level in a spectrum.
     *
     * @param spectrum Input spectrum
     * @return Estimated noise level
     */
    double estimateNoise(const Spectrum& spectrum) const;

    /**
     * @brief Get/set options
     */
    const PeakPickingOptions& options() const { return options_; }
    void setOptions(const PeakPickingOptions& opt) { options_ = opt; }

private:
    PeakPickingOptions options_;

    /// Find local maxima indices
    std::vector<Index> findLocalMaxima(const Spectrum& spectrum) const;

    /// Determine peak boundaries for a local maximum
    std::pair<Index, Index> findPeakBoundaries(const Spectrum& spectrum,
                                                Index apex) const;

    /// Calculate centroid position using weighted average
    MZ calculateCentroid(const Spectrum& spectrum, Index left,
                         Index right) const;

    /// Calculate peak area (trapezoidal integration)
    double calculateArea(const Spectrum& spectrum, Index left,
                         Index right) const;

    /// Fit Gaussian to peak region
    bool fitGaussian(const Spectrum& spectrum, Index left, Index right,
                     double& center, double& sigma, double& amplitude) const;

    /// Calculate FWHM from peak region
    MZ calculateFWHM(const Spectrum& spectrum, Index apex, Index left,
                     Index right) const;
};

/**
 * @brief Convenience function for peak picking.
 */
inline PeakList pickPeaks(const Spectrum& spectrum,
                          const PeakPickingOptions& options = {}) {
    PeakPicker picker(options);
    return picker.pick(spectrum);
}

/**
 * @brief Convenience function for centroiding.
 */
inline Spectrum centroidSpectrum(const Spectrum& spectrum,
                                 const PeakPickingOptions& options = {}) {
    PeakPicker picker(options);
    return picker.centroid(spectrum);
}

} // namespace algorithms
} // namespace lcms
