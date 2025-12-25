#include "../testing_utils.h"
#include <tuple>
#include <xsf/legendre.h>

TEST_CASE("assoc_legendre_p scipy/gh-23101", "[assoc_legendre_p][xsf_tests]") {
    using test_case = std::tuple<int, int, double, double, double>;
    // Reference values were computed with the Python library mpmath.
    auto [n, m, z, ref, rtol] = GENERATE(
        test_case{1, 0, -1.0, -1.0, 1e-15}, test_case{1, 0, std::nextafter(-1.0, -2.0), -1.0000000000000002, 1e-15},
        test_case{1, 0, std::nextafter(-1.0, 0.0), -0.9999999999999999, 1e-15}, test_case{2, 0, -1.0, 1.0, 1e-15},
        test_case{2, 0, std::nextafter(-1.0, -2.0), 1.0000000000000007, 1e-15},
        test_case{2, 0, std::nextafter(-1.0, 0.0), 0.9999999999999997, 1e-15}, test_case{3, 0, -1.0, -1.0, 1e-15},
        test_case{3, 0, std::nextafter(-1.0, -2.0), -1.0000000000000013, 1e-15},
        test_case{3, 0, std::nextafter(-1.0, 0.0), -0.9999999999999993, 1e-15}, test_case{4, 0, -1.0, 1.0, 1e-15},
        test_case{4, 0, std::nextafter(-1.0, -2.0), 1.0000000000000022, 4e-15},
        test_case{4, 0, std::nextafter(-1.0, 0.0), 0.9999999999999989, 1e-15}
    );
    const double w = xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, 2);
    const auto rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(n, m, z, w, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("assoc_legendre_p norm m0 gh-78", "[assoc_legendre_p][xsf_tests]") {
    const int n_max = 10;
    const int m = 0;
    const int num_points = 1000;
    const double left = -1.0;
    const double right = 1.0;

    std::vector<double> z(num_points);
    for (int i = 0; i < num_points; ++i) {
        z[i] = left + (right - left) * i / (num_points - 1);
    }

    for (int n = 0; n <= n_max; ++n) {
        for (const auto z_val : z) {
            // Compute unnormalized and normalized versions
            const double leg_p = xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z_val, 1);
            const double leg_p_norm = xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z_val, 1);

            // Expected relationship: p_norm = sqrt((2*n + 1) / 2) * p
            const double factor = std::sqrt((2.0 * n + 1.0) / 2.0);
            const double expected = factor * leg_p;

            const double rel_error = xsf::extended_relative_error(leg_p_norm, expected);
            CAPTURE(n, m, z_val, leg_p, leg_p_norm, expected, factor, rel_error);
            REQUIRE(rel_error <= 1e-8);
        }
    }
}
