/* C++ translation of tests from scipy/special/tests/test_hypergeometric.py
 * https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_hypergeometric.py
 */

#include "../testing_utils.h"

#include <cmath>
#include <limits>
#include <tuple>
#include <xsf/cpu/hyperu.h>

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

TEST_CASE("hyperu tiny a", "[hyperu][xsf_tests]") {
    // Regression test for gh-350: the first (tiny a) or second (a near -1)
    // term of the large-x series rounds to 0, and log10(0) was converted to
    // int, so this returned NaN.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    for a, b, x in [(1e-17, 0.5, 100.0), (-1e-17, 0.5, 100.0), (-1 + 2**-53, 0.5, 200.0)]:
        print(f"hyperu({a!r}, {b!r}, {x!r}) = {float(mpmath.hyperu(a, b, x))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, b, x, desired] = GENERATE(
        test_case{1e-17, 0.5, 100.0, 1.0}, test_case{-1e-17, 0.5, 100.0, 1.0},
        test_case{-0.9999999999999999, 0.5, 200.0, 199.4999999999999}
    );

    const double out = xsf::hyperu(a, b, x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-16;
    CAPTURE(a, b, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("hyperu small a small-x series", "[hyperu][xsf_tests]") {
    // Regression test for gh-436: the small-x series rounded `a + 1.0 - 1.0`,
    // which is inaccurate for small `a`.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for a, b, x in [(1e-15, 0.5, 50.0), (1e-10, 0.5, 30.0), (-1e-12, 2.5, 25.0)]:
        print(f"hyperu({a!r}, {b!r}, {x!r}) = {float(mpmath.hyperu(a, b, x))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, b, x, desired] = GENERATE(
        test_case{1e-15, 0.5, 50.0, 0.9999999999999961}, test_case{1e-10, 0.5, 30.0, 0.9999999996582531},
        test_case{-1e-12, 2.5, 25.0, 1.0000000000031584}
    );

    const double out = xsf::hyperu(a, b, x);
    const double error = xsf::extended_relative_error(out, desired);
    // the two series cancel, so allow more than a few eps
    const double rtol = 1e-11;
    CAPTURE(a, b, x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
