#include <array>
#include <cmath>
#include <complex>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <xsf/cpu/legendre.h>

#define MDSPAN_USE_PAREN_OPERATOR 1
#include <xsf/third_party/kokkos/mdspan.hpp>

TEST_CASE("lqn complex recurrence near the unit circle", "[lqn][complex]") {
    // scipy/scipy#26452: the recurrence depth could be negative, leaving
    // its normalization factor zero and returning NaNs and infinities.
    const std::array<std::complex<double>, 7> points = {std::complex<double>(0.0, 1.01),
                                                        {0.02, 1.02},
                                                        {-0.02, 1.02},
                                                        {-1.01, 0.0001},
                                                        {1.01, 0.0001},
                                                        {0.0, 0.9},
                                                        {0.0, 1.5}};
    for (auto point : points) {
        for (auto z : {point, std::conj(point)}) {
            for (int n : {1, 2, 10, 25}) {
                CAPTURE(z, n);
                std::vector<std::complex<double>> q(n + 1), d(n + 1);
                xsf::lqn(z, std::mdspan(q.data(), n + 1), std::mdspan(d.data(), n + 1));
                for (int k = 0; k <= n; ++k) {
                    CAPTURE(k);
                    REQUIRE(std::isfinite(q[k].real()));
                    REQUIRE(std::isfinite(q[k].imag()));
                    REQUIRE(std::isfinite(d[k].real()));
                    REQUIRE(std::isfinite(d[k].imag()));
                }
                // Preserve the existing branch convention on either side of
                // the unit circle, and check against the closed forms.
                double sign = std::abs(z) > 1 ? -1 : 1;
                auto q0 = std::log(sign * (1.0 + z) / (1.0 - z)) / 2.0;
                auto d0 = 1.0 / (1.0 - z * z);
                REQUIRE(std::abs(q[0] - q0) <= 1e-13 * std::abs(q0));
                REQUIRE(std::abs(q[1] - (z * q0 - 1.0)) <= 1e-13 * std::abs(z * q0 - 1.0));
                REQUIRE(std::abs(d[0] - d0) <= 1e-13 * std::abs(d0));
                REQUIRE(std::abs(d[1] - (q0 + z * d0)) <= 1e-13 * std::abs(q0 + z * d0));
            }
        }
    }
}

TEST_CASE("lqn higher complex degrees", "[lqn][complex]") {
    // mpmath.legenq(10, 0, z, type=3), with 60 decimal digits; derivatives
    // were evaluated with mpmath.diff at the same precision.
    const std::array<std::array<std::complex<double>, 3>, 3> cases = {
        {{{{0.0, 1.01}, {0.0, 2.857768011837148e-5}, {-2.18375040422366e-4, 0.0}}},
         {{{0.02, 1.02},
           {-4.008855535997653e-6, 2.615522369066548e-5},
           {-1.9859355603925883e-4, -3.2411989440671124e-5}}},
         {{{-1.01, 0.0001}, {-0.21784568075233668, -0.0021022847306103425}, {-21.020256704935, -0.3304735734924109}}}}
    };
    for (const auto &test : cases) {
        auto z = test[0];
        CAPTURE(z);
        std::array<std::complex<double>, 11> q, d;
        xsf::lqn(z, std::mdspan(q.data(), 11), std::mdspan(d.data(), 11));
        REQUIRE(std::abs(q[10] - test[1]) <= 1e-12 * std::abs(test[1]));
        REQUIRE(std::abs(d[10] - test[2]) <= 1e-12 * std::abs(test[2]));
    }
}
