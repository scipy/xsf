#include "../testing_utils.h"
#include <xsf/toms708.h>

// xsf::toms708::bratio only supports b = 1/2 (what the Student's t CDF needs).
TEST_CASE("toms708 bratio", "[bratio][xsf_tests]") {
    SECTION("regularized incomplete beta and its complement, b = 1/2") {
        // a, b, x, y = 1 - x, I_x(a, b), 1 - I_x(a, b)
        // Reference values computed with mpmath (dps=50), see mpmath.betainc(a, 0.5, 0, x, regularized=True).
        // The rows with large a are the Student's t CDF at t = 3 for nu = 2a.
        using test_case = std::tuple<double, double, double, double, double, double>;
        auto [a, b, x, y, expected, expected_c] = GENERATE(
            test_case{1e-20, 0.5, 0.3, 0.7, 1.0, 2.4198702426718918e-20},
            test_case{1e-20, 0.5, 0.8, 0.19999999999999996, 1.0, 9.624236501192067e-21},
            test_case{0.01, 0.5, 0.05, 0.95, 0.9575227020613565, 0.04247729793864357},
            test_case{0.01, 0.5, 0.9, 0.09999999999999998, 0.9935423956747499, 0.006457604325250123},
            test_case{0.3, 0.5, 0.1, 0.9, 0.3712349819272945, 0.6287650180727056},
            test_case{0.3, 0.5, 0.95, 0.050000000000000044, 0.9006314139336828, 0.0993685860663172},
            test_case{0.5, 0.5, 0.5, 0.5, 0.5, 0.5},
            test_case{1.0, 0.5, 0.2, 0.8, 0.10557280900008413, 0.8944271909999159},
            test_case{3.0, 0.5, 0.05, 0.95, 3.981768725535919e-05, 0.9999601823127446},
            test_case{3.0, 0.5, 0.9, 0.09999999999999998, 0.4454155553479705, 0.5545844446520295},
            test_case{10, 0.5, 0.9, 0.09999999999999998, 0.15164090963470997, 0.84835909036529},
            test_case{20.0, 0.5, 0.02, 0.98, 1.3273091787181266e-35, 1.0},
            test_case{20.0, 0.5, 0.95, 0.050000000000000044, 0.1545907814334381, 0.8454092185665619},
            test_case{
                50000.0, 0.5, 0.9999100080992711, 8.999190072893439e-05, 0.002700460884064719, 0.9972995391159353
            },
            test_case{1000000.0, 0.5, 0.999999, 1.0000000000287557e-06, 0.15729915515583529, 0.8427008448441647},
            test_case{500000000.0, 0.5, 0.999999991, 8.999999919e-09, 0.0026997961297379157, 0.9973002038702621}
        );
        const auto result = xsf::toms708::bratio(a, b, x, y);
        const auto rel_error = xsf::extended_relative_error(result.d1, expected);
        const auto rel_error_c = xsf::extended_relative_error(result.d2, expected_c);
        CAPTURE(a, b, x, y, result.d1, result.d2, expected, expected_c, rel_error, rel_error_c);
        REQUIRE(result.i1 == 0);
        REQUIRE(rel_error <= 1e-13);
        REQUIRE(rel_error_c <= 1e-13);
    }

    SECTION("unsupported or invalid arguments set the error flag") {
        REQUIRE(xsf::toms708::bratio(2.0, 3.0, 0.5, 0.5).i1 == 8); // b != 1/2
        REQUIRE(xsf::toms708::bratio(0.0, 0.5, 0.0, 1.0).i1 == 6); // x = a = 0
    }

    SECTION("endpoints") {
        auto lo = xsf::toms708::bratio(2.0, 0.5, 0.0, 1.0);
        REQUIRE(lo.d1 == 0.0);
        REQUIRE(lo.d2 == 1.0);
        auto hi = xsf::toms708::bratio(2.0, 0.5, 1.0, 0.0);
        REQUIRE(hi.d1 == 1.0);
        REQUIRE(hi.d2 == 0.0);
        auto zero_a = xsf::toms708::bratio(0.0, 0.5, 0.5, 0.5);
        REQUIRE(zero_a.d1 == 1.0);
        REQUIRE(zero_a.d2 == 0.0);
    }
}
