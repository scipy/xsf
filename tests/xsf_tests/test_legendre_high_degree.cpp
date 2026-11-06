#include <array>
#include <cmath>
#include <complex>
#include <tuple>
#include <type_traits>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <xsf/cpu/legendre.h>
#include <xsf/cpu/sph_harm.h>

namespace {
template <typename T>
void require_close(T actual, T expected, double rtol) {
    CAPTURE(actual, expected);
    REQUIRE(std::abs(actual - expected) <= rtol * std::abs(expected));
}
} // namespace

TEMPLATE_TEST_CASE(
    "normalized Legendre high degree scipy/gh-26574", "[assoc_legendre_p][sph_legendre_p][xsf_tests]", float, double,
    std::complex<double>
) {
    using T = TestType;
    const double rtol = std::is_same_v<T, float> ? 5e-4 : 1e-10;
    for (int n : {645, 646, 1000, 1235}) {
        for (double angle : {0.2, 0.75, 1.0}) {
            CAPTURE(n, angle);
            T theta(angle);
            if constexpr (!std::is_floating_point_v<T>) {
                theta += T(0, 0.1);
            }
            const T z = std::cos(theta);
            const T p = xsf::legendre_p(n, z);
            const T norm = std::sqrt(T((2.0 * n + 1.0) / 2.0));
            require_close(xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, 0, z, 2), norm * p, rtol);
            if constexpr (std::is_floating_point_v<T>) {
                require_close(xsf::sph_legendre_p(n, 0, theta), norm * p / T(std::sqrt(2.0 * M_PI)), rtol);
            }
        }
    }
}

TEST_CASE("normalized Legendre high order scipy/gh-26574", "[assoc_legendre_p][sph_harm_y][xsf_tests]") {
    using test_case = std::tuple<int, int, double>;
    // mpmath, 60 decimal digits: sqrt((2*n+1)/2 * (n-m)!/(n+m)!) * legenp(n, m, 0.25).
    const std::array<test_case, 4> cases = {
        test_case{802, 800, 1.723635946013349242034372e-9},
        test_case{802, -800, 1.723635946013349242034372e-9},
        test_case{1235, 1200, 0.609699707843673049492721},
        test_case{1235, -1200, 0.609699707843673049492721},
    };
    for (auto [n, m, ref] : cases) {
        CAPTURE(n, m);
        require_close(xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, 0.25, 2), ref, 1e-10);
        const auto harmonic = xsf::sph_harm_y(n, m, std::acos(0.25), 0.0);
        require_close(harmonic.real(), ref / std::sqrt(2.0 * M_PI), 1e-10);
        REQUIRE(harmonic.imag() == 0.0);
    }
}

TEMPLATE_TEST_CASE(
    "normalized Legendre high degree derivatives scipy/gh-26574", "[assoc_legendre_p][sph_legendre_p][xsf_tests]",
    double, std::complex<double>
) {
    using T = TestType;
    const xsf::dual<T, 2> z{T(0.75), T(1)};
    const xsf::dual<T, 2> theta{T(0.75), T(1)};
    for (int n : {645, 646, 1000}) {
        const T norm = std::sqrt(T((2.0 * n + 1.0) / 2.0));
        const auto assoc_ref = xsf::legendre_p(n, z) * norm;
        const auto assoc = xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, 0, z, 2);
        for (unsigned k = 0; k <= 2; ++k) {
            CAPTURE(n, k);
            require_close(assoc[k], assoc_ref[k], 1e-10);
        }
        if constexpr (std::is_floating_point_v<T>) {
            const auto spherical_ref = xsf::legendre_p(n, xsf::cos(theta)) * (norm / T(std::sqrt(2.0 * M_PI)));
            const auto spherical = xsf::sph_legendre_p(n, 0, theta);
            for (unsigned k = 0; k <= 2; ++k) {
                CAPTURE(n, k);
                require_close(spherical[k], spherical_ref[k], 1e-10);
            }
        }
    }
}
