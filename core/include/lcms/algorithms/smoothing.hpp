#pragma once

#include "../spectrum.hpp"
#include "../chromatogram.hpp"
#include <vector>
#include <cmath>

namespace lcms {
namespace algorithms {

/**
 * @brief Smoothing method
 */
enum class SmoothingMethod : std::uint8_t {
    GAUSSIAN,           // Gaussian kernel smoothing
    SAVITZKY_GOLAY,     // Savitzky-Golay filter
    MOVING_AVERAGE,     // Simple moving average
    LOESS,              // Local polynomial regression
    WHITTAKER           // Whittaker smoother
};

/**
 * @brief Options for smoothing
 */
struct SmoothingOptions {
    /// Smoothing method
    SmoothingMethod method = SmoothingMethod::GAUSSIAN;

    /// Window size (number of points, must be odd for SG)
    int window_size = 5;

    /// Gaussian sigma (in points) for Gaussian smoothing
    double gaussian_sigma = 1.5;

    /// Polynomial order for Savitzky-Golay (must be < window_size)
    int sg_order = 2;

    /// Lambda parameter for Whittaker smoother
    double whittaker_lambda = 2.0;

    /// Number of passes for iterative methods
    int passes = 1;

    /// Preserve peak heights (scale after smoothing)
    bool preserve_height = false;
};

/**
 * @brief Smoothing algorithms for spectra and chromatograms.
 *
 * Provides various smoothing/filtering algorithms to reduce noise
 * while preserving signal features.
 */
class Smoother {
public:
    Smoother() = default;
    explicit Smoother(const SmoothingOptions& options)
        : options_(options) {}

    /**
     * @brief Smooth a spectrum.
     *
     * @param spectrum Input spectrum
     * @return Smoothed spectrum
     */
    Spectrum smooth(const Spectrum& spectrum);

    /**
     * @brief Smooth a chromatogram.
     *
     * @param chromatogram Input chromatogram
     * @return Smoothed chromatogram
     */
    Chromatogram smooth(const Chromatogram& chromatogram);

    /**
     * @brief Smooth intensity values.
     *
     * @param intensity Input intensities
     * @return Smoothed intensities
     */
    std::vector<double> smooth(const std::vector<double>& intensity);

    /**
     * @brief Get/set options
     */
    const SmoothingOptions& options() const { return options_; }
    void setOptions(const SmoothingOptions& opt) { options_ = opt; }

    /**
     * @brief Generate Savitzky-Golay coefficients.
     *
     * @param window_size Window size (must be odd)
     * @param order Polynomial order
     * @param derivative Derivative order (0 = smoothing)
     * @return Filter coefficients
     */
    static std::vector<double> computeSGCoefficients(int window_size,
                                                      int order,
                                                      int derivative = 0);

private:
    SmoothingOptions options_;

    /// Gaussian smoothing
    std::vector<double> gaussianSmooth(const std::vector<double>& data) const;

    /// Savitzky-Golay smoothing
    std::vector<double> savitzkyGolaySmooth(const std::vector<double>& data) const;

    /// Moving average smoothing
    std::vector<double> movingAverageSmooth(const std::vector<double>& data) const;

    /// Generate Gaussian kernel
    std::vector<double> gaussianKernel(int half_width, double sigma) const;

    /// Convolve signal with kernel
    std::vector<double> convolve(const std::vector<double>& signal,
                                  const std::vector<double>& kernel) const;
};

/**
 * @brief Convenience function for spectrum smoothing.
 */
inline Spectrum smoothSpectrum(const Spectrum& spectrum,
                               const SmoothingOptions& options = {}) {
    Smoother smoother(options);
    return smoother.smooth(spectrum);
}

/**
 * @brief Convenience function for chromatogram smoothing.
 */
inline Chromatogram smoothChromatogram(const Chromatogram& chromatogram,
                                       const SmoothingOptions& options = {}) {
    Smoother smoother(options);
    return smoother.smooth(chromatogram);
}

/**
 * @brief Apply Gaussian smoothing to a vector.
 */
inline std::vector<double> gaussianSmooth(const std::vector<double>& data,
                                          int window_size = 5,
                                          double sigma = 1.5) {
    SmoothingOptions opts;
    opts.method = SmoothingMethod::GAUSSIAN;
    opts.window_size = window_size;
    opts.gaussian_sigma = sigma;
    Smoother smoother(opts);
    return smoother.smooth(data);
}

/**
 * @brief Apply Savitzky-Golay smoothing to a vector.
 */
inline std::vector<double> savitzkyGolaySmooth(const std::vector<double>& data,
                                               int window_size = 5,
                                               int order = 2) {
    SmoothingOptions opts;
    opts.method = SmoothingMethod::SAVITZKY_GOLAY;
    opts.window_size = window_size;
    opts.sg_order = order;
    Smoother smoother(opts);
    return smoother.smooth(data);
}

} // namespace algorithms
} // namespace lcms
