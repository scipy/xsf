#include "../testing_utils.h"

#include <cmath>
#include <xsf/sphd_wave.h>

// Regression tests for gh-241: all _cv spheroidal wave functions must reject
// inputs where (n - m) > 198 by returning NaN, rather than crashing.

TEST_CASE("prolate_aswfa rejects n - m > 198", "[sphd_wave][xsf_tests]") {
    double s1f, s1d;

    // n - m = 220 > 198: must return NaN, not crash
    xsf::prolate_aswfa(0.0, 220.0, 220.0, 1.0, 0.5, s1f, s1d);
    REQUIRE(std::isnan(s1f));
    REQUIRE(std::isnan(s1d));

    // n - m = 199 > 198: boundary, must also return NaN
    xsf::prolate_aswfa(0.0, 199.0, 10.0, 1.0, 0.5, s1f, s1d);
    REQUIRE(std::isnan(s1f));
    REQUIRE(std::isnan(s1d));
}

TEST_CASE("prolate_aswfa accepts n - m == 198", "[sphd_wave][xsf_tests]") {
    double s1f, s1d;

    xsf::prolate_aswfa(0.0, 198.0, 10.0, 1.0, 0.5, s1f, s1d);
    REQUIRE(std::isfinite(s1f));
    REQUIRE(std::isfinite(s1d));
}

TEST_CASE("oblate_aswfa rejects n - m > 198", "[sphd_wave][xsf_tests]") {
    double s1f, s1d;

    xsf::oblate_aswfa(0.0, 220.0, 220.0, 1.0, 0.5, s1f, s1d);
    REQUIRE(std::isnan(s1f));
    REQUIRE(std::isnan(s1d));

    xsf::oblate_aswfa(0.0, 199.0, 10.0, 1.0, 0.5, s1f, s1d);
    REQUIRE(std::isnan(s1f));
    REQUIRE(std::isnan(s1d));
}

TEST_CASE("oblate_aswfa accepts n - m == 198", "[sphd_wave][xsf_tests]") {
    double s1f, s1d;

    xsf::oblate_aswfa(0.0, 198.0, 10.0, 1.0, 0.5, s1f, s1d);
    REQUIRE(std::isfinite(s1f));
    REQUIRE(std::isfinite(s1d));
}

TEST_CASE("prolate_radial1 rejects n - m > 198", "[sphd_wave][xsf_tests]") {
    double r1f, r1d;

    // Exact reproducer from gh-241: pro_rad1_cv(0, 220, 220.0, 1.0, 1.1)
    xsf::prolate_radial1(0.0, 220.0, 220.0, 1.0, 1.1, r1f, r1d);
    REQUIRE(std::isnan(r1f));
    REQUIRE(std::isnan(r1d));

    xsf::prolate_radial1(0.0, 199.0, 10.0, 1.0, 1.1, r1f, r1d);
    REQUIRE(std::isnan(r1f));
    REQUIRE(std::isnan(r1d));
}

TEST_CASE("prolate_radial1 accepts n - m == 198", "[sphd_wave][xsf_tests]") {
    double r1f, r1d;

    xsf::prolate_radial1(0.0, 198.0, 10.0, 1.0, 1.1, r1f, r1d);
    REQUIRE(std::isfinite(r1f));
    REQUIRE(std::isfinite(r1d));
}

TEST_CASE("prolate_radial2 rejects n - m > 198", "[sphd_wave][xsf_tests]") {
    double r2f, r2d;

    xsf::prolate_radial2(0.0, 220.0, 220.0, 1.0, 1.1, r2f, r2d);
    REQUIRE(std::isnan(r2f));
    REQUIRE(std::isnan(r2d));

    xsf::prolate_radial2(0.0, 199.0, 10.0, 1.0, 1.1, r2f, r2d);
    REQUIRE(std::isnan(r2f));
    REQUIRE(std::isnan(r2d));
}

TEST_CASE("prolate_radial2 accepts n - m == 198", "[sphd_wave][xsf_tests]") {
    double r2f, r2d;

    xsf::prolate_radial2(0.0, 198.0, 10.0, 1.0, 1.1, r2f, r2d);
    REQUIRE(std::isfinite(r2f));
    REQUIRE(std::isfinite(r2d));
}

TEST_CASE("oblate_radial1 rejects n - m > 198", "[sphd_wave][xsf_tests]") {
    double r1f, r1d;

    xsf::oblate_radial1(0.0, 220.0, 220.0, 1.0, 1.0, r1f, r1d);
    REQUIRE(std::isnan(r1f));
    REQUIRE(std::isnan(r1d));

    xsf::oblate_radial1(0.0, 199.0, 10.0, 1.0, 1.0, r1f, r1d);
    REQUIRE(std::isnan(r1f));
    REQUIRE(std::isnan(r1d));
}

TEST_CASE("oblate_radial1 accepts n - m == 198", "[sphd_wave][xsf_tests]") {
    double r1f, r1d;

    xsf::oblate_radial1(0.0, 198.0, 10.0, 1.0, 1.0, r1f, r1d);
    REQUIRE(std::isfinite(r1f));
    REQUIRE(std::isfinite(r1d));
}

TEST_CASE("oblate_radial2 rejects n - m > 198", "[sphd_wave][xsf_tests]") {
    double r2f, r2d;

    xsf::oblate_radial2(0.0, 220.0, 220.0, 1.0, 1.0, r2f, r2d);
    REQUIRE(std::isnan(r2f));
    REQUIRE(std::isnan(r2d));

    xsf::oblate_radial2(0.0, 199.0, 10.0, 1.0, 1.0, r2f, r2d);
    REQUIRE(std::isnan(r2f));
    REQUIRE(std::isnan(r2d));
}

TEST_CASE("oblate_radial2 accepts n - m == 198", "[sphd_wave][xsf_tests]") {
    double r2f, r2d;

    xsf::oblate_radial2(0.0, 198.0, 10.0, 1.0, 1.0, r2f, r2d);
    REQUIRE(std::isfinite(r2f));
    REQUIRE(std::isfinite(r2d));
}

TEST_CASE("cv functions reject n - m > 198 with nonzero m", "[sphd_wave][xsf_tests]") {
    double out0, out1;
