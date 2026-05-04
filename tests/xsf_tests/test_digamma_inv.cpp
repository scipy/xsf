#include "../testing_utils.h"
#include <xsf/cephes/psi.h>
#include <xsf/digamma_inv.h>

#include <cmath>
#include <limits>

TEST_CASE("digamma_inv round-trip", "[digamma_inv][xsf_tests]") {
    const double rtol = 50 * std::numeric_limits<double>::epsilon();
    const std::vector<double> xs = logspace<double>(-20, 20, 500);
    for (double x : xs) {
        const double y = xsf::cephes::psi(x);
        const auto result = xsf::digamma_inv(y);
        const auto relative_error = xsf::extended_relative_error(result, x);
        CAPTURE(x, y, result, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("digamma_inv edge cases", "[digamma_inv][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double log_max_double = std::log(std::numeric_limits<double>::max());

    REQUIRE(std::isnan(xsf::digamma_inv(nan)));
    REQUIRE(xsf::digamma_inv(inf) == inf);
    REQUIRE(xsf::digamma_inv(-inf) == 0.0);
    REQUIRE(xsf::digamma_inv(std::nextafter(log_max_double, inf)) == inf);
}

TEST_CASE("digamma_inv precision", "[digamma_inv][xsf_tests]") {
    // Some test cases taken from https://en.wikipedia.org/wiki/Digamma_function#Special_values
    using test_case = std::tuple<double, double>;
    // Reference values for digamma were computed with the Python library mpmath.
    auto [x, ref_digamma] = GENERATE(
        test_case{1.0, -0.5772156649015329}, test_case{1. / 2., -1.9635100260214235},
        test_case{1. / 3., -3.1320337800208065}, test_case{1. / 4., -4.2274535333762655},
        test_case{1.461632144968362341, 0.0}, // positive root of digamma
        test_case{10, 2.251752589066721}, test_case{1e5, 11.512920464961896}, test_case{1e50, 115.12925464970229},
        test_case{1e200, 460.51701859880916}, test_case{1e-3, -1000.5755719318103},
        test_case{1e-1, -10.423754940411076}, test_case{1e-10, -10000000000.577215}, test_case{1e-20, -1e+20},
        test_case{1e-50, -1e+50}, test_case{1e-100, -1e+100}
    );
    double x_result = xsf::digamma_inv(ref_digamma);
    const auto rel_error = xsf::extended_relative_error(x_result, x);

    CAPTURE(x, ref_digamma, x_result, rel_error);
    REQUIRE(rel_error <= 1e-13);
}
