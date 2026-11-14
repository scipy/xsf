#include "../testing_utils.h"
#include <tuple>
#include <xsf/cpu/sphd_wave.h>

TEST_CASE("prolate radial Wronskian", "[sphd_wave][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double>;
    // R1 R2' - R1' R2 = 1 / (c (x^2 - 1))
    auto [m, n, c, x] =
        GENERATE(test_case{1.0, 4.0, 1.0, 2.0}, test_case{2.0, 4.0, 15.0, 1.1}, test_case{1.0, 3.0, 8.0, 5.0});
    constexpr double rtol = 1e-10;
    double r1f, r1d, r2f, r2d;
    xsf::prolate_radial1_nocv(m, n, c, x, r1f, r1d);
    xsf::prolate_radial2_nocv(m, n, c, x, r2f, r2d);
    const auto wronskian = r1f * r2d - r1d * r2f;
    const auto expected = 1.0 / (c * (x * x - 1.0));
    const auto rel_error = xsf::extended_relative_error(wronskian, expected);
    CAPTURE(m, n, c, x, r1f, r1d, r2f, r2d, wronskian, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("oblate radial Wronskian", "[sphd_wave][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double>;
    // R1 R2' - R1' R2 = 1 / (c (x^2 + 1))
    auto [m, n, c, x] = GENERATE(test_case{1.0, 4.0, 0.5, 2.0}, test_case{2.0, 3.0, 3.0, 5.0});
    constexpr double rtol = 1e-10;
    double r1f, r1d, r2f, r2d;
    xsf::oblate_radial1_nocv(m, n, c, x, r1f, r1d);
    xsf::oblate_radial2_nocv(m, n, c, x, r2f, r2d);
    const auto wronskian = r1f * r2d - r1d * r2f;
    const auto expected = 1.0 / (c * (x * x + 1.0));
    const auto rel_error = xsf::extended_relative_error(wronskian, expected);
    CAPTURE(m, n, c, x, r1f, r1d, r2f, r2d, wronskian, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}
