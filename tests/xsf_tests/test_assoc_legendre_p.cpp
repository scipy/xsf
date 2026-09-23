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

TEST_CASE("assoc_legendre_p second derivative sign at z = +-1 scipy/gh-26229", "[assoc_legendre_p][xsf_tests]") {
    // The second derivative of P_n^m at x = -1 is (-1)^n times its value at x = 1, so for
    // odd n the z = -1 endpoint must flip sign. References were computed by differentiating
    // P_n^m(x) = (-1)^m (1 - x^2)^(m/2) * d^m P_n(x) / dx^m symbolically and evaluating at the endpoints.
    using test_case = std::tuple<int, int, double, double, double>;
    auto [n, m, z, ref, rtol] = GENERATE(
        test_case{3, 0, -1.0, -15.0, 1e-13}, test_case{3, 0, 1.0, 15.0, 1e-13},
        test_case{4, 0, -1.0, 45.0, 1e-13}, test_case{4, 0, 1.0, 45.0, 1e-13},
        test_case{5, 0, -1.0, -105.0, 1e-13},
        test_case{3, 2, -1.0, 90.0, 1e-12}, test_case{3, 2, 1.0, -90.0, 1e-12},
        test_case{4, 2, -1.0, -510.0, 1e-12}, test_case{5, 2, -1.0, 1890.0, 1e-12},
        test_case{7, 0, -1.0, -378.0, 1e-12}, test_case{7, 2, -1.0, 13356.0, 1e-11},
        test_case{5, 4, -1.0, -7560.0, 1e-11}
    );
    xsf::dual<double, 2> z_dual(z);
    const auto w = xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z_dual, 2);
    const auto rel_error = xsf::extended_relative_error(w[2], ref);

    CAPTURE(n, m, z, (double)w[2], ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("assoc_legendre_p norm m0 gh-78", "[assoc_legendre_p][xsf_tests]") {
    const int n_max = 10;
    const int m = 0;
    const int num_points = 1000;
    const double left = -1.0;
    const double right = 1.0;

    const std::vector<double> z = linspace(left, right, num_points);

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
