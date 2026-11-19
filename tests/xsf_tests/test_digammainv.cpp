#include "../testing_utils.h"
#include "xsf/config.h"
#include <random>
#include <xsf/cephes/psi.h>
#include <xsf/digammainv.h>

TEST_CASE("digammainv round-trip", "[digammainv][xsf_tests]") {
    const double rtol = 50 * std::numeric_limits<double>::epsilon();
    const std::vector<double> xs = logspace<double>(-20, 20, 500);
    for (double x : xs) {
        const double y = xsf::cephes::psi(x);
        const auto result = xsf::digammainv(y);
        const auto relative_error = xsf::extended_relative_error(result, x);
        CAPTURE(x, y, result, rtol, relative_error);
        REQUIRE(relative_error <= rtol);
    }
}

TEST_CASE("digammainv edge cases", "[digammainv][xsf_tests]") {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double log_max_double = std::log(std::numeric_limits<double>::max());

    REQUIRE(std::isnan(xsf::digammainv(nan)));
    REQUIRE(xsf::digammainv(inf) == inf);
    REQUIRE(xsf::digammainv(-inf) == 0.0);
    REQUIRE(xsf::digammainv(std::nextafter(log_max_double, inf)) == inf);
}

TEST_CASE("digammainv precision", "[digammainv][xsf_tests]") {
    // Some test cases taken from https://en.wikipedia.org/wiki/Digamma_function#Special_values
    using test_case = std::tuple<double, double>;
    // Reference values for digamma were computed with the Python library mpmath.
    auto [x, ref_digamma] = GENERATE(
        test_case{1e-100, -1e+100}, test_case{1e-50, -1e+50}, test_case{1e-20, -1e+20},
        test_case{1e-10, -10000000000.577215}, test_case{1e-3, -1000.5755719318103},
        test_case{1e-1, -10.423754940411076}, test_case{1. / 4., -4.2274535333762655},
        test_case{1. / 3., -3.1320337800208065}, test_case{1. / 2., -1.9635100260214235},
        test_case{1.0, -0.5772156649015329}, test_case{1.461632144968362341, 0.0}, // positive root of digamma
        test_case{10, 2.251752589066721}, test_case{1e5, 11.512920464961896}, test_case{1e50, 115.12925464970229},
        test_case{1e200, 460.51701859880916}
    );
    double x_result = xsf::digammainv(ref_digamma);
    const auto rel_error = xsf::extended_relative_error(x_result, x);

    CAPTURE(x, ref_digamma, x_result, rel_error);
    REQUIRE(rel_error <= 1e-13);
}

TEST_CASE("random_inputs", "[digammainv][xsf_tests]") {
    // round trip test with random values of y in [0, 40]
    const double y_min = 0.0;
    const double y_max = 40.0;
    std::mt19937 gen(42088764423867);
    std::uniform_real_distribution<double> dist(y_min, y_max);
    for (int i = 0; i < 10000; i++) {
        double y = dist(gen);
        double x = xsf::digammainv(y);
        double y_check = xsf::cephes::psi(x);
        CAPTURE(y, x, y_check);
        REQUIRE(xsf::extended_relative_error(y_check, y) <= 1e-13);
    }
}

TEST_CASE("gh-361", "[digammainv][xsf_tests]") {
    // gh-361 reported that digammainv yielded nan for certain inputs
    const double y = 9.220027678716601;
    const double x = xsf::digammainv(y);
    const double y_check = xsf::cephes::psi(x);
    CAPTURE(y, x, y_check);
    REQUIRE(xsf::extended_relative_error(y_check, y) <= 1e-14);
}
