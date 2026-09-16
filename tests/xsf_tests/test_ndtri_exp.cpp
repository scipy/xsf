/* C++ translation of tests from scipy/special/tests/test_ndtri_exp.py
 * https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_ndtri_exp.py
 */

#include "../testing_utils.h"

#include <xsf/ndtri_exp.h>
#include <xsf/stats.h>

namespace {

void check_roundtrip(double y, double rtol) {
    // Check that the roundtrip from y to ndtri_exp(y) and then to log_ndtr is accurate.
    const double out = xsf::log_ndtr(xsf::ndtri_exp(y));
    const double error = xsf::extended_relative_error(out, y);
    CAPTURE(y, out, rtol, error);
    REQUIRE(error <= rtol);
}

} // namespace

TEST_CASE("ndtri_exp very small arguments", "[ndtri_exp][xsf_tests]") {
    SET_FP_FORMAT()
    const auto scale = GENERATE(-1e1, -1e2, -1e10, -1e20, -std::numeric_limits<double>::max());

    for (const auto point : linspace(1e-10, 1.0 - 1e-10, 1000)) {
        check_roundtrip(scale * (0.5 * point + 0.5), 1e-14);
    }
}

TEST_CASE("ndtri_exp intervals", "[ndtri_exp][xsf_tests]") {
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double>;
    auto [left, right, rtol] = GENERATE(
        test_case{-10.0, -2.0, 1e-14}, test_case{-2.0, -0.14542, 1e-12}, test_case{-0.14542, -1e-6, 1e-10},
        test_case{-1e-6, 0.0, 1e-6}
    );

    for (const auto point : linspace(1e-10, 1.0 - 1e-10, 1000)) {
        check_roundtrip((right - left) * point + left, rtol);
    }
}

TEST_CASE("ndtri_exp extreme", "[ndtri_exp][xsf_tests]") {
    SET_FP_FORMAT()
    double bigneg = std::numeric_limits<double>::lowest();

    for (int i = 0; i < 4; ++i) {
        bigneg = std::nextafter(bigneg, 0.0);
    }
    check_roundtrip(-std::numeric_limits<double>::min(), 1e-12);
    check_roundtrip(bigneg, 1e-12);
}

TEST_CASE("ndtri_exp asymptotes", "[ndtri_exp][xsf_tests]") {
    REQUIRE(xsf::ndtri_exp(-std::numeric_limits<double>::infinity()) == -std::numeric_limits<double>::infinity());
    REQUIRE(xsf::ndtri_exp(0.0) == std::numeric_limits<double>::infinity());
}

TEST_CASE("ndtri_exp outside domain", "[ndtri_exp][xsf_tests]") { REQUIRE(std::isnan(xsf::ndtri_exp(1.0))); }
