#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "lcms/peak.hpp"

using namespace lcms;
using Catch::Approx;

TEST_CASE("Peak construction", "[peak]") {
    SECTION("Default construction") {
        Peak peak;
        REQUIRE(peak.mz() == Approx(0.0));
        REQUIRE(peak.intensity() == Approx(0.0));
        REQUIRE(peak.charge() == 0);
    }

    SECTION("Construction with parameters") {
        Peak peak(500.0, 10000.0, 120.0);
        REQUIRE(peak.mz() == Approx(500.0));
        REQUIRE(peak.intensity() == Approx(10000.0));
        REQUIRE(peak.rt() == Approx(120.0));
    }
}

TEST_CASE("Peak neutral mass", "[peak]") {
    Peak peak;
    peak.setMz(500.2635);  // [M+2H]2+ of a ~998 Da peptide

    SECTION("Unknown charge returns 0") {
        REQUIRE(peak.neutralMass() == Approx(0.0));
    }

    SECTION("Calculate neutral mass with charge") {
        peak.setCharge(2);
        double neutral = peak.neutralMass();
        REQUIRE(neutral == Approx(998.5124).margin(0.01));
    }
}

TEST_CASE("Peak contains", "[peak]") {
    Peak peak;
    peak.setMz(500.0);
    peak.setRt(120.0);
    peak.setMzLeft(499.5);
    peak.setMzRight(500.5);
    peak.setRtLeft(115.0);
    peak.setRtRight(125.0);

    SECTION("Contains m/z") {
        REQUIRE(peak.containsMz(500.0));
        REQUIRE(peak.containsMz(499.5));
        REQUIRE(peak.containsMz(500.5));
        REQUIRE_FALSE(peak.containsMz(499.0));
        REQUIRE_FALSE(peak.containsMz(501.0));
    }

    SECTION("Contains RT") {
        REQUIRE(peak.containsRt(120.0));
        REQUIRE(peak.containsRt(115.0));
        REQUIRE_FALSE(peak.containsRt(110.0));
    }

    SECTION("Contains 2D") {
        REQUIRE(peak.contains(500.0, 120.0));
        REQUIRE_FALSE(peak.contains(500.0, 110.0));
        REQUIRE_FALSE(peak.contains(501.0, 120.0));
    }
}

TEST_CASE("PeakList operations", "[peaklist]") {
    PeakList peaks;

    SECTION("Add peaks") {
        peaks.add(Peak(100.0, 1000.0));
        peaks.add(Peak(200.0, 5000.0));
        peaks.add(Peak(300.0, 2000.0));

        REQUIRE(peaks.size() == 3);
    }

    SECTION("Sort by m/z") {
        peaks.add(Peak(300.0, 2000.0));
        peaks.add(Peak(100.0, 1000.0));
        peaks.add(Peak(200.0, 5000.0));

        peaks.sortByMz();
        REQUIRE(peaks[0].mz() == Approx(100.0));
        REQUIRE(peaks[1].mz() == Approx(200.0));
        REQUIRE(peaks[2].mz() == Approx(300.0));
    }

    SECTION("Sort by intensity") {
        peaks.add(Peak(100.0, 1000.0));
        peaks.add(Peak(200.0, 5000.0));
        peaks.add(Peak(300.0, 2000.0));

        peaks.sortByIntensity();
        REQUIRE(peaks[0].intensity() == Approx(5000.0));
        REQUIRE(peaks[1].intensity() == Approx(2000.0));
        REQUIRE(peaks[2].intensity() == Approx(1000.0));
    }

    SECTION("Find in m/z range") {
        peaks.add(Peak(100.0, 1000.0));
        peaks.add(Peak(200.0, 5000.0));
        peaks.add(Peak(300.0, 2000.0));

        auto filtered = peaks.findInMzRange(150.0, 250.0);
        REQUIRE(filtered.size() == 1);
        REQUIRE(filtered[0].mz() == Approx(200.0));
    }

    SECTION("Find nearest m/z") {
        peaks.add(Peak(100.0, 1000.0));
        peaks.add(Peak(200.0, 5000.0));
        peaks.add(Peak(300.0, 2000.0));

        const Peak* nearest = peaks.findNearestMz(210.0);
        REQUIRE(nearest != nullptr);
        REQUIRE(nearest->mz() == Approx(200.0));
    }
}
