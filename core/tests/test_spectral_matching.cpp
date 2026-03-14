/**
 * Tests for spectral matching and similarity scoring (C++).
 * Uses Catch2 testing framework.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

#include "lcms/algorithms/spectral_matching.hpp"

using namespace lcms;
using Catch::Matchers::WithinAbs;

TEST_CASE("Cosine similarity", "[matching]") {
    SECTION("Identical spectra score 1.0") {
        std::vector<double> mz = {100, 200, 300};
        std::vector<double> ints = {1000, 2000, 500};
        double score = cosineSimilarity(mz, ints, mz, ints, 0.01);
        REQUIRE_THAT(score, WithinAbs(1.0, 0.01));
    }

    SECTION("Non-overlapping spectra score ~0") {
        std::vector<double> mz1 = {100, 200};
        std::vector<double> ints1 = {1000, 2000};
        std::vector<double> mz2 = {300, 400};
        std::vector<double> ints2 = {1000, 2000};
        double score = cosineSimilarity(mz1, ints1, mz2, ints2, 0.01);
        REQUIRE(score < 0.1);
    }

    SECTION("Similar spectra score high") {
        std::vector<double> mz1 = {100, 200, 300};
        std::vector<double> ints1 = {1000, 2000, 500};
        std::vector<double> mz2 = {100.001, 200.001, 300.001};
        std::vector<double> ints2 = {950, 2050, 480};
        double score = cosineSimilarity(mz1, ints1, mz2, ints2, 0.01);
        REQUIRE(score > 0.9);
    }
}

TEST_CASE("Spectral library search", "[matching]") {
    SpectralLibrary lib;

    std::vector<double> refMz = {100, 200, 300};
    std::vector<double> refInts = {1000, 2000, 500};
    lib.addSpectrum("compound1", 400.0, refMz, refInts);

    REQUIRE(lib.size() == 1);

    auto matches = lib.search(refMz, refInts, 400.0, 5, 0.01);
    REQUIRE(matches.size() >= 1);
    REQUIRE(matches[0].score > 0.8);
}
