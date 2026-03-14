#include "lcms/algorithms/spectral_matching.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>

namespace lcms {
namespace algorithms {

namespace {

struct MatchedPair {
    double int1;
    double int2;
};

std::vector<MatchedPair> matchPeaks(
    const std::vector<MZ>& mz1, const std::vector<Intensity>& int1,
    const std::vector<MZ>& mz2, const std::vector<Intensity>& int2,
    double tolerance) {

    if (mz1.empty() || mz2.empty()) return {};

    std::vector<MatchedPair> matched;
    std::vector<bool> used2(mz2.size(), false);

    for (size_t i = 0; i < mz1.size(); ++i) {
        int best_j = -1;
        double best_diff = tolerance + 1.0;

        for (size_t j = 0; j < mz2.size(); ++j) {
            if (used2[j]) continue;
            double diff = std::abs(mz1[i] - mz2[j]);
            if (diff <= tolerance && diff < best_diff) {
                best_diff = diff;
                best_j = static_cast<int>(j);
            }
        }

        if (best_j >= 0) {
            matched.push_back({int1[i], int2[best_j]});
            used2[best_j] = true;
        }
    }

    return matched;
}

double vectorNorm(const std::vector<Intensity>& v) {
    double sum = 0.0;
    for (double x : v) sum += x * x;
    return std::sqrt(sum);
}

} // anonymous namespace

std::pair<double, int> cosineSimilarity(
    const std::vector<MZ>& mz1, const std::vector<Intensity>& int1,
    const std::vector<MZ>& mz2, const std::vector<Intensity>& int2,
    double tolerance) {

    auto matched = matchPeaks(mz1, int1, mz2, int2, tolerance);
    if (matched.empty()) return {0.0, 0};

    double dot = 0.0;
    for (const auto& m : matched) {
        dot += m.int1 * m.int2;
    }

    double norm1 = vectorNorm(int1);
    double norm2 = vectorNorm(int2);

    if (norm1 == 0.0 || norm2 == 0.0) return {0.0, 0};

    double score = dot / (norm1 * norm2);
    return {std::min(score, 1.0), static_cast<int>(matched.size())};
}

std::pair<double, int> modifiedCosineSimilarity(
    const std::vector<MZ>& mz1, const std::vector<Intensity>& int1, MZ precursor1,
    const std::vector<MZ>& mz2, const std::vector<Intensity>& int2, MZ precursor2,
    double tolerance) {

    double mass_diff = precursor1 - precursor2;

    auto direct = matchPeaks(mz1, int1, mz2, int2, tolerance);

    // Shift mz2 by mass difference
    std::vector<MZ> shifted_mz2(mz2.size());
    for (size_t i = 0; i < mz2.size(); ++i) {
        shifted_mz2[i] = mz2[i] + mass_diff;
    }
    auto shifted = matchPeaks(mz1, int1, shifted_mz2, int2, tolerance);

    // Combine
    std::vector<MatchedPair> all_matched;
    all_matched.insert(all_matched.end(), direct.begin(), direct.end());
    all_matched.insert(all_matched.end(), shifted.begin(), shifted.end());

    if (all_matched.empty()) return {0.0, 0};

    double dot = 0.0;
    for (const auto& m : all_matched) {
        dot += m.int1 * m.int2;
    }

    double norm1 = vectorNorm(int1);
    double norm2 = vectorNorm(int2);

    if (norm1 == 0.0 || norm2 == 0.0) return {0.0, 0};

    double score = dot / (norm1 * norm2);
    return {std::min(score, 1.0), static_cast<int>(all_matched.size())};
}

std::pair<double, int> spectralEntropySimilarity(
    const std::vector<MZ>& mz1, const std::vector<Intensity>& int1,
    const std::vector<MZ>& mz2, const std::vector<Intensity>& int2,
    double tolerance) {

    // Normalize to probabilities
    double sum1 = 0.0, sum2 = 0.0;
    for (double v : int1) sum1 += v;
    for (double v : int2) sum2 += v;

    if (sum1 == 0.0 || sum2 == 0.0) return {0.0, 0};

    std::vector<Intensity> p1(int1.size()), p2(int2.size());
    for (size_t i = 0; i < int1.size(); ++i) p1[i] = int1[i] / sum1;
    for (size_t i = 0; i < int2.size(); ++i) p2[i] = int2[i] / sum2;

    auto matched = matchPeaks(mz1, p1, mz2, p2, tolerance);
    if (matched.empty()) return {0.0, 0};

    // Shannon entropy
    auto entropy = [](const std::vector<Intensity>& p) {
        double h = 0.0;
        for (double v : p) {
            if (v > 0) h -= v * std::log(v);
        }
        return h;
    };

    double h1 = entropy(p1);
    double h2 = entropy(p2);

    // Merged spectrum entropy
    std::vector<double> merged;
    double merged_sum = 0.0;
    for (const auto& m : matched) {
        double avg = (m.int1 + m.int2) / 2.0;
        merged.push_back(avg);
        merged_sum += avg;
    }
    if (merged_sum > 0) {
        for (auto& v : merged) v /= merged_sum;
    }

    std::vector<Intensity> merged_i(merged.begin(), merged.end());
    double h_merged = entropy(merged_i);

    double avg_entropy = (h1 + h2) / 2.0;
    if (avg_entropy == 0.0) return {0.0, static_cast<int>(matched.size())};

    double jsd = h_merged - avg_entropy;
    double max_jsd = std::log(2.0);
    double similarity = 1.0 - std::min(jsd / max_jsd, 1.0);

    return {std::max(0.0, similarity), static_cast<int>(matched.size())};
}

size_t SpectralLibrary::loadMGF(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) return 0;

    size_t count = 0;
    bool in_ions = false;
    LibrarySpectrum current;

    std::string line;
    while (std::getline(file, line)) {
        // Trim whitespace
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' '))
            line.pop_back();

        if (line == "BEGIN IONS") {
            in_ions = true;
            current = LibrarySpectrum();
        } else if (line == "END IONS") {
            if (!current.mz.empty()) {
                entries_.push_back(current);
                count++;
            }
            in_ions = false;
        } else if (in_ions) {
            auto eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value = line.substr(eq_pos + 1);
                if (key == "TITLE") {
                    current.name = value;
                } else if (key == "PEPMASS") {
                    std::istringstream iss(value);
                    iss >> current.precursor_mz;
                } else {
                    current.metadata[key] = value;
                }
            } else {
                std::istringstream iss(line);
                double mz, intensity;
                if (iss >> mz >> intensity) {
                    current.mz.push_back(mz);
                    current.intensity.push_back(intensity);
                }
            }
        }
    }

    return count;
}

size_t SpectralLibrary::loadMSP(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) return 0;

    size_t count = 0;
    LibrarySpectrum current;
    bool reading_peaks = false;

    std::string line;
    while (std::getline(file, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();

        if (line.empty()) {
            if (!current.mz.empty()) {
                entries_.push_back(current);
                count++;
            }
            current = LibrarySpectrum();
            reading_peaks = false;
            continue;
        }

        if (!reading_peaks) {
            auto colon_pos = line.find(':');
            if (colon_pos != std::string::npos) {
                std::string key = line.substr(0, colon_pos);
                std::string value = line.substr(colon_pos + 1);
                // Trim leading spaces from value
                while (!value.empty() && value[0] == ' ') value.erase(0, 1);

                if (key == "Name" || key == "NAME") {
                    current.name = value;
                } else if (key == "PrecursorMZ" || key == "PRECURSORMZ") {
                    current.precursor_mz = std::stod(value);
                } else if (key == "Num Peaks" || key == "NUM PEAKS") {
                    reading_peaks = true;
                } else {
                    current.metadata[key] = value;
                }
            }
        } else {
            std::istringstream iss(line);
            double mz, intensity;
            if (iss >> mz >> intensity) {
                current.mz.push_back(mz);
                current.intensity.push_back(intensity);
            }
        }
    }

    // Handle last entry
    if (!current.mz.empty()) {
        entries_.push_back(current);
        count++;
    }

    return count;
}

void SpectralLibrary::saveMGF(const std::string& filepath) const {
    std::ofstream file(filepath);
    for (const auto& entry : entries_) {
        file << "BEGIN IONS\n";
        file << "TITLE=" << entry.name << "\n";
        file << "PEPMASS=" << entry.precursor_mz << "\n";
        for (const auto& [key, value] : entry.metadata) {
            file << key << "=" << value << "\n";
        }
        for (size_t i = 0; i < entry.mz.size(); ++i) {
            file << entry.mz[i] << " " << entry.intensity[i] << "\n";
        }
        file << "END IONS\n\n";
    }
}

std::vector<SpectralMatch> SpectralLibrary::search(
    const std::vector<MZ>& query_mz,
    const std::vector<Intensity>& query_intensity,
    MZ query_precursor,
    const std::string& method,
    double tolerance,
    double min_score,
    int top_n,
    double precursor_tolerance) const {

    std::vector<SpectralMatch> results;

    for (const auto& entry : entries_) {
        if (precursor_tolerance > 0 && query_precursor > 0) {
            if (std::abs(entry.precursor_mz - query_precursor) > precursor_tolerance) {
                continue;
            }
        }

        double score = 0.0;
        int n_matched = 0;

        if (method == "cosine") {
            auto [s, n] = cosineSimilarity(query_mz, query_intensity,
                                            entry.mz, entry.intensity, tolerance);
            score = s;
            n_matched = n;
        } else if (method == "modified_cosine") {
            auto [s, n] = modifiedCosineSimilarity(
                query_mz, query_intensity, query_precursor,
                entry.mz, entry.intensity, entry.precursor_mz, tolerance);
            score = s;
            n_matched = n;
        } else if (method == "entropy") {
            auto [s, n] = spectralEntropySimilarity(
                query_mz, query_intensity,
                entry.mz, entry.intensity, tolerance);
            score = s;
            n_matched = n;
        }

        if (score >= min_score) {
            SpectralMatch match;
            match.name = entry.name;
            match.score = score;
            match.matched_peaks = n_matched;
            match.precursor_mz = entry.precursor_mz;
            match.metadata = entry.metadata;
            results.push_back(match);
        }
    }

    std::sort(results.begin(), results.end(),
              [](const SpectralMatch& a, const SpectralMatch& b) {
                  return a.score > b.score;
              });

    if (static_cast<int>(results.size()) > top_n) {
        results.resize(top_n);
    }

    return results;
}

} // namespace algorithms
} // namespace lcms
