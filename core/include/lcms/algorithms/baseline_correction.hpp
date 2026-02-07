#pragma once

#include "../spectrum.hpp"
#include "../chromatogram.hpp"
#include <vector>
#include <cmath>

namespace lcms {
namespace algorithms {

/**
 * @brief Baseline correction method
 */
enum class BaselineMethod : std::uint8_t {
    TOPHAT,         // Morphological top-hat filter
    SNIP,           // Statistics-sensitive Non-linear Iterative Peak-clipping
    ROLLINGBALL,    // Rolling ball algorithm
    POLYNOMIAL,     // Polynomial fitting
    ASYMMETRIC_LS   // Asymmetric least squares
};

/**
 * @brief Options for baseline correction
 */
struct BaselineCorrectionOptions {
    /// Correction method
    BaselineMethod method = BaselineMethod::SNIP;

    /// Structuring element half-width for morphological methods (m/z or RT)
    double half_window = 50.0;

    /// Number of iterations for SNIP
    int snip_iterations = 40;

    /// Smoothing factor for SNIP
    bool snip_smooth = true;

    /// Ball radius for rolling ball
    double ball_radius = 100.0;

    /// Polynomial degree for polynomial fitting
    int polynomial_degree = 3;

    /// Lambda parameter for asymmetric LS
    double als_lambda = 1e6;

    /// Asymmetry parameter for asymmetric LS (0-1)
    double als_p = 0.01;

    /// Number of iterations for asymmetric LS
    int als_iterations = 10;

    /// Return baseline instead of corrected signal
    bool return_baseline = false;
};

/**
 * @brief Baseline correction algorithms.
 *
 * Removes baseline drift from spectra and chromatograms using
 * various morphological and statistical methods.
 */
class BaselineCorrector {
public:
    BaselineCorrector() = default;
    explicit BaselineCorrector(const BaselineCorrectionOptions& options)
        : options_(options) {}

    /**
     * @brief Correct baseline of a spectrum.
     *
     * @param spectrum Input spectrum
     * @return Baseline-corrected spectrum (or baseline if return_baseline=true)
     */
    Spectrum correct(const Spectrum& spectrum);

    /**
     * @brief Correct baseline of a chromatogram.
     *
     * @param chromatogram Input chromatogram
     * @return Baseline-corrected chromatogram
     */
    Chromatogram correct(const Chromatogram& chromatogram);

    /**
     * @brief Estimate baseline of intensity values.
     *
     * @param intensity Intensity values
     * @param spacing Average spacing between points
     * @return Estimated baseline values
     */
    std::vector<double> estimateBaseline(const std::vector<double>& intensity,
                                         double spacing = 1.0);

    /**
     * @brief Get/set options
     */
    const BaselineCorrectionOptions& options() const { return options_; }
    void setOptions(const BaselineCorrectionOptions& opt) { options_ = opt; }

private:
    BaselineCorrectionOptions options_;

    /// Top-hat morphological baseline estimation
    std::vector<double> topHatBaseline(const std::vector<double>& intensity,
                                        int half_window) const;

    /// SNIP baseline estimation
    std::vector<double> snipBaseline(const std::vector<double>& intensity,
                                      int iterations) const;

    /// Rolling ball baseline estimation
    std::vector<double> rollingBallBaseline(const std::vector<double>& intensity,
                                             int radius) const;

    /// Morphological erosion
    std::vector<double> erode(const std::vector<double>& data,
                               int half_window) const;

    /// Morphological dilation
    std::vector<double> dilate(const std::vector<double>& data,
                                int half_window) const;
};

/**
 * @brief Convenience function for baseline correction of spectrum.
 */
inline Spectrum correctBaseline(const Spectrum& spectrum,
                                const BaselineCorrectionOptions& options = {}) {
    BaselineCorrector corrector(options);
    return corrector.correct(spectrum);
}

/**
 * @brief Convenience function for baseline correction of chromatogram.
 */
inline Chromatogram correctBaseline(const Chromatogram& chromatogram,
                                    const BaselineCorrectionOptions& options = {}) {
    BaselineCorrector corrector(options);
    return corrector.correct(chromatogram);
}

} // namespace algorithms
} // namespace lcms
