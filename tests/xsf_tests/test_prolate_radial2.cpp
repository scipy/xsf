#include "../testing_utils.h"

#include <xsf/cpu/sphd_wave.h>

TEST_CASE("prolate_radial2 odd m", "[prolate_radial2][xsf_tests]") {
    // Regression test for gh-551: lpmns applied the Ferrers phase (-1)^m for x > 1,
    // making prolate_radial2 wrong for odd m. Check the Wronskian
    // R1 * R2' - R1' * R2 = 1 / (c * (x**2 - 1)).
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [m, n, c, x] = GENERATE(test_case{1, 1, 1, 1.1}, test_case{1, 1, 10, 1.1}, test_case{1, 2, 5, 1.5});

    double r1f, r1d, r2f, r2d;
    xsf::prolate_radial1_nocv(m, n, c, x, r1f, r1d);
    xsf::prolate_radial2_nocv(m, n, c, x, r2f, r2d);
    const double out = (r1f * r2d - r1d * r2f) * c * (x * x - 1);
    const double error = xsf::extended_relative_error(out, 1.0);
    const double rtol = 1e-9;
    CAPTURE(m, n, c, x, r1f, r1d, r2f, r2d, out, error, rtol);
    REQUIRE(error <= rtol);
}
