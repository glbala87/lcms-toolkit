#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <optional>
#include <limits>
#include <cmath>

namespace lcms {

/// Mass-to-charge ratio type
using MZ = double;

/// Intensity value type
using Intensity = double;

/// Retention time in seconds
using RetentionTime = double;

/// Index type for spectra and peaks
using Index = std::size_t;

/// MS level (1 = MS1, 2 = MS/MS, etc.)
using MSLevel = std::uint8_t;

/// Charge state
using ChargeState = std::int8_t;

/// Polarity of the ion mode
enum class Polarity : std::int8_t {
    UNKNOWN = 0,
    POSITIVE = 1,
    NEGATIVE = -1
};

/// Activation method for MS/MS
enum class ActivationMethod : std::uint8_t {
    UNKNOWN = 0,
    CID,     // Collision-Induced Dissociation
    HCD,     // Higher-energy Collisional Dissociation
    ETD,     // Electron Transfer Dissociation
    ECD,     // Electron Capture Dissociation
    UVPD,    // Ultraviolet Photodissociation
    IRMPD    // Infrared Multiphoton Dissociation
};

/// Spectrum type
enum class SpectrumType : std::uint8_t {
    UNKNOWN = 0,
    PROFILE,    // Profile/continuum data
    CENTROID    // Centroided/peak-picked data
};

/// Peak shape type
enum class PeakShape : std::uint8_t {
    UNKNOWN = 0,
    GAUSSIAN,
    LORENTZIAN,
    EMG  // Exponentially Modified Gaussian
};

/// Precursor information for MS/MS spectra
struct Precursor {
    MZ mz = 0.0;
    Intensity intensity = 0.0;
    ChargeState charge = 0;
    MZ isolation_window_lower = 0.0;
    MZ isolation_window_upper = 0.0;
    ActivationMethod activation = ActivationMethod::UNKNOWN;
    double collision_energy = 0.0;

    bool hasCharge() const { return charge != 0; }
    MZ isolationWindowWidth() const {
        return isolation_window_upper - isolation_window_lower;
    }
};

/// Range template for min/max values
template<typename T>
struct Range {
    T min_value = std::numeric_limits<T>::max();
    T max_value = std::numeric_limits<T>::lowest();

    Range() = default;
    Range(T min_val, T max_val) : min_value(min_val), max_value(max_val) {}

    bool isEmpty() const { return min_value > max_value; }
    T span() const { return max_value - min_value; }
    T center() const { return (min_value + max_value) / 2; }

    bool contains(T value) const {
        return value >= min_value && value <= max_value;
    }

    void extend(T value) {
        if (value < min_value) min_value = value;
        if (value > max_value) max_value = value;
    }

    void extend(const Range<T>& other) {
        extend(other.min_value);
        extend(other.max_value);
    }

    bool overlaps(const Range<T>& other) const {
        return min_value <= other.max_value && max_value >= other.min_value;
    }
};

using MZRange = Range<MZ>;
using RTRange = Range<RetentionTime>;
using IntensityRange = Range<Intensity>;

/// Key-value metadata container
using MetaData = std::map<std::string, std::string>;

/// Tolerance for m/z matching
struct MZTolerance {
    double value = 0.0;
    bool is_ppm = false;

    MZTolerance() = default;
    MZTolerance(double val, bool ppm = false) : value(val), is_ppm(ppm) {}

    static MZTolerance Da(double daltons) { return MZTolerance(daltons, false); }
    static MZTolerance PPM(double ppm) { return MZTolerance(ppm, true); }

    double absoluteAt(MZ mz) const {
        return is_ppm ? mz * value * 1e-6 : value;
    }

    bool matches(MZ mz1, MZ mz2) const {
        return std::abs(mz1 - mz2) <= absoluteAt((mz1 + mz2) / 2);
    }
};

/// Status codes for operations
enum class Status : std::uint8_t {
    OK = 0,
    ERROR,
    FILE_NOT_FOUND,
    PARSE_ERROR,
    INVALID_FORMAT,
    OUT_OF_MEMORY,
    INVALID_PARAMETER
};

/// Convert polarity to string
inline std::string toString(Polarity p) {
    switch (p) {
        case Polarity::POSITIVE: return "positive";
        case Polarity::NEGATIVE: return "negative";
        default: return "unknown";
    }
}

/// Convert activation method to string
inline std::string toString(ActivationMethod a) {
    switch (a) {
        case ActivationMethod::CID: return "CID";
        case ActivationMethod::HCD: return "HCD";
        case ActivationMethod::ETD: return "ETD";
        case ActivationMethod::ECD: return "ECD";
        case ActivationMethod::UVPD: return "UVPD";
        case ActivationMethod::IRMPD: return "IRMPD";
        default: return "unknown";
    }
}

/// Convert spectrum type to string
inline std::string toString(SpectrumType t) {
    switch (t) {
        case SpectrumType::PROFILE: return "profile";
        case SpectrumType::CENTROID: return "centroid";
        default: return "unknown";
    }
}

} // namespace lcms
