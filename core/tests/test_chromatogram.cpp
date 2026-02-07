#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "lcms/chromatogram.hpp"

using namespace lcms;
using Catch::Approx;

TEST_CASE("Chromatogram construction", "[chromatogram]") {
    SECTION("Default construction") {
        Chromatogram chrom;
        REQUIRE(chrom.size() == 0);
        REQUIRE(chrom.empty());
    }

    SECTION("Construction with data") {
        std::vector<RetentionTime> rt = {0.0, 60.0, 120.0};
        std::vector<Intensity> intensity = {100.0, 500.0, 200.0};

        Chromatogram chrom(rt, intensity);
        REQUIRE(chrom.size() == 3);
        REQUIRE_FALSE(chrom.empty());
    }
}

TEST_CASE("Chromatogram statistics", "[chromatogram]") {
    std::vector<RetentionTime> rt = {0.0, 60.0, 120.0, 180.0};
    std::vector<Intensity> intensity = {100.0, 500.0, 300.0, 150.0};
    Chromatogram chrom(rt, intensity);

    SECTION("Maximum intensity") {
        REQUIRE(chrom.maxIntensity() == Approx(500.0));
    }

    SECTION("Apex RT") {
        REQUIRE(chrom.rtAtMaxIntensity() == Approx(60.0));
    }

    SECTION("RT range") {
        REQUIRE(chrom.rtRange().min_value == Approx(0.0));
        REQUIRE(chrom.rtRange().max_value == Approx(180.0));
    }
}

TEST_CASE("Chromatogram area calculation", "[chromatogram]") {
    std::vector<RetentionTime> rt = {0.0, 60.0, 120.0};
    std::vector<Intensity> intensity = {100.0, 100.0, 100.0};
    Chromatogram chrom(rt, intensity);

    SECTION("Trapezoidal integration") {
        // Rectangle with height 100 and width 120
        double area = chrom.computeArea();
        REQUIRE(area == Approx(12000.0));
    }
}

TEST_CASE("Chromatogram interpolation", "[chromatogram]") {
    std::vector<RetentionTime> rt = {0.0, 100.0};
    std::vector<Intensity> intensity = {0.0, 1000.0};
    Chromatogram chrom(rt, intensity);

    SECTION("Interpolate at midpoint") {
        REQUIRE(chrom.interpolateAt(50.0) == Approx(500.0));
    }

    SECTION("Interpolate at boundary") {
        REQUIRE(chrom.interpolateAt(0.0) == Approx(0.0));
        REQUIRE(chrom.interpolateAt(100.0) == Approx(1000.0));
    }
}

TEST_CASE("Chromatogram extract range", "[chromatogram]") {
    std::vector<RetentionTime> rt = {0.0, 60.0, 120.0, 180.0, 240.0};
    std::vector<Intensity> intensity = {100.0, 200.0, 300.0, 200.0, 100.0};
    Chromatogram chrom(rt, intensity);

    SECTION("Extract middle") {
        auto extracted = chrom.extractRange(60.0, 180.0);
        REQUIRE(extracted.size() == 3);
        REQUIRE(extracted.rtAt(0) == Approx(60.0));
        REQUIRE(extracted.rtAt(2) == Approx(180.0));
    }
}
