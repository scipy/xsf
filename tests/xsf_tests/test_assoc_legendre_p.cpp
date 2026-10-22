#include "../testing_utils.h"
#include <tuple>
#include <xsf/cpu/legendre.h>

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

TEST_CASE("assoc_legendre_p negative order second derivative at z = +-1", "[assoc_legendre_p][xsf_tests]") {
    // Reference values were computed with SymPy by differentiating
    // P_n^m(x) = (n - |m|)! / (n + |m|)! * (1 - x^2)^(|m|/2) d^|m| P_n(x) / dx^|m| (branch cut 2)
    // twice and evaluating at the endpoints. For branch cut 3, (1 - x^2) is replaced by (x^2 - 1).
    using test_case = std::tuple<int, int, double, double>;
    auto [n, m, z, ref] = GENERATE(
        test_case{2, -2, 1.0, -1.0 / 4}, test_case{2, -2, -1.0, -1.0 / 4}, test_case{3, -2, 1.0, -3.0 / 4},
        test_case{3, -2, -1.0, 3.0 / 4}, test_case{4, -2, 1.0, -17.0 / 12}, test_case{4, -2, -1.0, -17.0 / 12},
        test_case{7, -2, 1.0, -53.0 / 12}, test_case{7, -2, -1.0, 53.0 / 12}, test_case{4, -4, 1.0, 1.0 / 48},
        test_case{4, -4, -1.0, 1.0 / 48}, test_case{5, -4, 1.0, 1.0 / 48}, test_case{5, -4, -1.0, -1.0 / 48},
        test_case{8, -4, 1.0, 1.0 / 48}, test_case{8, -4, -1.0, 1.0 / 48}
    );
    const double rtol = 1e-14;

    const auto w = xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, xsf::dual<double, 2>(z), 2);
    const auto rel_error = xsf::extended_relative_error(w[2], ref);
    CAPTURE(n, m, z, w[2], ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);

    // Branch cut 3 flips the sign for m = -2, since (x^2 - 1) = -(1 - x^2).
    const std::complex<double> ref_cut3 = (m == -2) ? -ref : ref;
    const auto w_cut3 = xsf::assoc_legendre_p(
        xsf::assoc_legendre_unnorm, n, m, xsf::dual<std::complex<double>, 2>(std::complex<double>(z, 0)), 3
    );
    const auto rel_error_cut3 = xsf::extended_relative_error(w_cut3[2], ref_cut3);
    CAPTURE(w_cut3[2], ref_cut3, rel_error_cut3);
    REQUIRE(rel_error_cut3 <= rtol);
}

TEST_CASE("assoc_legendre_p norm derivatives at z = +-1", "[assoc_legendre_p][xsf_tests]") {
    // At the endpoints, the derivatives of the normalized polynomial must be
    // sqrt((2 n + 1) / 2 * (n - m)! / (n + m)!) times those of the unnormalized one.
    for (int n = 0; n <= 10; ++n) {
        for (int m = -std::min(n, 4); m <= std::min(n, 4); ++m) {
            const double factor =
                std::sqrt((2.0 * n + 1.0) / 2.0 * std::exp(std::lgamma(n - m + 1.0) - std::lgamma(n + m + 1.0)));
            for (const double z : {-1.0, 1.0}) {
                const auto leg_p = xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, xsf::dual<double, 2>(z), 2);
                const auto leg_p_norm =
                    xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, xsf::dual<double, 2>(z), 2);
                for (int k = 1; k <= 2; ++k) {
                    const double expected = factor * leg_p[k];
                    const double rel_error = xsf::extended_relative_error(leg_p_norm[k], expected);
                    CAPTURE(n, m, z, k, leg_p[k], leg_p_norm[k], expected, rel_error);
                    REQUIRE(rel_error <= 1e-13);
                }
            }
        }
    }
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
