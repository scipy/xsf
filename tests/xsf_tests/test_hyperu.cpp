/* C++ translation of tests from scipy/special/tests/test_hypergeometric.py
 * https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_hypergeometric.py
 */

#include "../testing_utils.h"

#include <cmath>
#include <limits>
#include <tuple>
#include <xsf/hyperu.h>

TEST_CASE("hyperu negative x", "[hyperu][xsf_tests]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (const double a : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
        for (const double b : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
            for (const double x : linspace(-100.0, -1.0, 10)) {
                CAPTURE(a, b, x);
                REQUIRE(std::isnan(xsf::hyperu(a, b, x)));
            }
        }
    }

    REQUIRE(std::isnan(xsf::hyperu(nan, 1.0, -1.0)));
}

TEST_CASE("hyperu special cases", "[hyperu][xsf_tests]") { REQUIRE(xsf::hyperu(0.0, 1.0, 1.0) == 1.0); }

TEST_CASE("hyperu nan inputs", "[hyperu][xsf_tests]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (const double a : {0.5, 1.0, nan}) {
        for (const double b : {1.0, 2.0, nan}) {
            for (const double x : {0.25, 3.0, nan}) {
                const bool expected_nan = std::isnan(a) || std::isnan(b) || std::isnan(x);
                CAPTURE(a, b, x);
                REQUIRE(std::isnan(xsf::hyperu(a, b, x)) == expected_nan);
            }
        }
    }
}

TEST_CASE("hyperu gh-15650 mpmath reference values", "[hyperu][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double>;

    auto [a, b, x, expected] = GENERATE(
        test_case{0.21581740448533887, 1.0, 1e-05, 3.6030558839391325},
        test_case{0.21581740448533887, 1.0, 0.00021544346900318823, 2.8783254988948976},
        test_case{0.21581740448533887, 1.0, 0.004641588833612777, 2.154928216691109},
        test_case{0.21581740448533887, 1.0, 0.1, 1.446546638718792},
        test_case{0.0030949064301273865, 1.0, 1e-05, 1.0356696454116199},
        test_case{0.0030949064301273865, 1.0, 0.00021544346900318823, 1.0261510362481985},
        test_case{0.0030949064301273865, 1.0, 0.004641588833612777, 1.0166326903402296},
        test_case{0.0030949064301273865, 1.0, 0.1, 1.0071174207698674},
        test_case{0.1509924314279033, 1.0, 1e-05, 2.806173846998948},
        test_case{0.1509924314279033, 1.0, 0.00021544346900318823, 2.3092158526816124},
        test_case{0.1509924314279033, 1.0, 0.004641588833612777, 1.812905980588048},
        test_case{0.1509924314279033, 1.0, 0.1, 1.3239738117634872},
        test_case{-0.010678995342969011, 1.0, 1e-05, 0.8775194903781114},
        test_case{-0.010678995342969011, 1.0, 0.00021544346900318823, 0.9101008998540128},
        test_case{-0.010678995342969011, 1.0, 0.004641588833612777, 0.9426854294058609},
        test_case{-0.010678995342969011, 1.0, 0.1, 0.9753065150174902},
        test_case{-0.06556622211831487, 1.0, 1e-05, 0.26435429752668904},
        test_case{-0.06556622211831487, 1.0, 0.00021544346900318823, 0.4574756033875781},
        test_case{-0.06556622211831487, 1.0, 0.004641588833612777, 0.6507121093358457},
        test_case{-0.06556622211831487, 1.0, 0.1, 0.8453129788602187},
        test_case{-0.21628242470175185, 1.0, 1e-05, -1.2318314201114489},
        test_case{-0.21628242470175185, 1.0, 0.00021544346900318823, -0.6704694233529538},
        test_case{-0.21628242470175185, 1.0, 0.004641588833612777, -0.10795098653682857},
        test_case{-0.21628242470175185, 1.0, 0.1, 0.4687227684115524}
    );

    const double result = xsf::hyperu(a, b, x);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(a, b, x, result, expected, error);
    REQUIRE(error <= 1e-13);
}

TEST_CASE("hyperu gh-15650 sanity", "[hyperu][xsf_tests]") {
    for (const double a : linspace(-0.5, 0.5, 500)) {
        for (const double x : linspace(1e-6, 1e-1, 500)) {
            const double result = xsf::hyperu(a, 1.0, x);
            CAPTURE(a, x, result);
            REQUIRE(std::abs(result) < 1e3);
        }
    }
}

/* gh-2287: chguit integrates exp(-x t) t^(a-1) (1+t)^(b-a-1) by composite
 * Gauss-Legendre. For non-integer a the factor t^(a-1) is an algebraic endpoint
 * singularity at t = 0, against which Gauss-Legendre converges only
 * algebraically, so refining the panel count does not recover the lost digits.
 * Reference values from mpmath at 50 decimal digits.
 *
 * Before the t = s^(1/a) substitution these inputs lost up to 7 digits; the
 * first entry below returned 0.05412514120003843 against a true
 * 0.054125064626921962.
 */
TEST_CASE("hyperu gh-2287 endpoint singularity mpmath reference values", "[hyperu][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double>;

    auto [a, b, x, expected] = GENERATE(
        test_case{1.094823, 0.503515, 12.863076, 0.054125064626921962},
        test_case{0.77, 0.7, 21.3, 0.091484136843193521},
        test_case{1.0057, 1.1009, 17.8632, 0.052523025175485847},
        test_case{0.35, 0.15, 15.0, 0.37766112753745428},
        test_case{1.2, 0.7, 9.5, 0.057080950328595002},
        test_case{0.6, 2.35, 25.0, 0.14754510472505289},
        test_case{1.3, 1.8, 30.0, 0.011768924460569544},
        test_case{0.9, 0.25, 6.0, 0.16353060681990762}
    );

    const double result = xsf::hyperu(a, b, x);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(a, b, x, result, expected, error);
    REQUIRE(error <= 1e-8);
}
