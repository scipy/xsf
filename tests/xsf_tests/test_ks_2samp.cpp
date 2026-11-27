#include "../testing_utils.h"
#include <algorithm>
#include <numeric>
#include <xsf/stats.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

namespace {

// The path counting functions take a caller-provided work array of length
// at least min(m, n) + 1.
std::vector<double> make_work(int64_t m, int64_t n) {
    return std::vector<double>(static_cast<std::size_t>(std::min(m, n) + 1));
}

} // namespace

/*
# Reference values computed with scipy.stats._stats_py and scipy.stats._stats_pythran
import math
import numpy as np
from scipy.stats._stats_py import (_attempt_exact_2kssamp, _compute_prob_outside_square,
                                   _count_paths_outside_method)
from scipy.stats._stats_pythran import _compute_outer_prob_inside_method

cases = [  # n1, n2, d, alternative
    (5, 5, 0.0, 'two-sided'), (5, 5, 0.6, 'two-sided'), (20, 20, 0.15, 'two-sided'),
    (100, 100, 0.25, 'two-sided'), (7, 3, 0.5, 'two-sided'), (40, 25, 0.3, 'two-sided'),
    (200, 150, 0.1, 'two-sided'), (1000, 999, 0.05, 'two-sided'),
    (5, 5, 0.6, 'less'), (20, 20, 0.15, 'greater'), (100, 100, 0.25, 'less'),
    (7, 3, 0.5, 'greater'), (40, 25, 0.3, 'less'), (200, 150, 0.1, 'greater'),
    (1000, 1, 1.0, 'less'), (1100, 1099, 1/(1100*1099), 'less'), (2000, 1000, 1/2000, 'greater'),
]
for n1, n2, d, alt in cases:
    g = math.gcd(n1, n2)
    with np.errstate(all='ignore'):
        s, dd, p = _attempt_exact_2kssamp(n1, n2, g, d, alt)
    print(f"{{{n1}, {n2}, {d!r}, {str(alt == 'two-sided').lower()}, "
          f"{str(bool(s)).lower()}, {float(dd)!r}, {float(p)!r}}},")

for n, h in [(1, 1), (5, 3), (20, 3), (100, 25), (500, 7)]:
    print(f"{{{n}, {h}, {_compute_prob_outside_square(n, h)!r}}},")

for m, n, h in [(1, 1, 1), (7, 3, 11), (40, 25, 60), (200, 150, 60), (150, 200, 60)]:
    g = math.gcd(m, n)
    print(f"{{{m}, {n}, {g}, {h}, {float(_compute_outer_prob_inside_method(m, n, g, h))!r}}},")

for m, n, h in [(7, 3, 11), (40, 25, 60), (200, 150, 60), (1000, 1, 1001), (3, 7, 11)]:
    g = math.gcd(m, n)
    print(f"{{{m}, {n}, {g}, {h}, {float(_count_paths_outside_method(m, n, g, h))!r}}},")

# The (40, 25) and (200, 150) one-sided cases use exact references instead, computed
# with integer arithmetic (SciPy's values inherit the roundoff error of binom).
from fractions import Fraction

def count_outside(m, n, h):
    g = math.gcd(m, n)
    A = [[0] * (n + 1) for _ in range(m + 1)]
    for x in range(m + 1):
        for y in range(n + 1):
            if n * x - m * y >= h * g:
                continue
            A[x][y] = 1 if x == y == 0 else (A[x - 1][y] if x else 0) + (A[x][y - 1] if y else 0)
    return math.comb(m + n, n) - A[m][n]

for m, n, h in [(40, 25, 60), (200, 150, 60)]:
    c = count_outside(m, n, h)
    print(f"{{{m}, {n}, {math.gcd(m, n)}, {h}, {float(c)!r}}},",
          float(Fraction(c, math.comb(m + n, n))))
*/

namespace {

struct AttemptExactCase {
    int64_t n1;
    int64_t n2;
    double d;
    bool two_sided;
    bool expected_success;
    double expected_d;
    double expected_prob;
    double rtol = 1e-14;
};

struct OutsideSquareCase {
    int64_t n;
    int64_t h;
    double expected;
};

struct LatticeCase {
    int64_t m;
    int64_t n;
    int64_t g;
    int64_t h;
    double expected;
    double rtol = 1e-14;
};

constexpr double qnan = std::numeric_limits<double>::quiet_NaN();

} // namespace

TEST_CASE("attempt_exact_2kssamp test", "[attempt_exact_2kssamp][xsf_tests]") {
    const std::vector<AttemptExactCase> cases = {
        {5, 5, 0.0, true, true, 0.0, 1.0},
        {5, 5, 0.6, true, true, 0.6, 0.35714285714285715},
        {20, 20, 0.15, true, true, 0.15, 0.9831368772656193},
        {100, 100, 0.25, true, true, 0.25, 0.0037294923618311367},
        {7, 3, 0.5, true, true, 0.47619047619047616, 0.5333333333333334},
        {40, 25, 0.3, true, true, 0.3, 0.1013339946846},
        {200, 150, 0.1, true, true, 0.1, 0.33868681219914437},
        {1000, 999, 0.05, true, true, 0.05, 0.15592065000426333},
        {5, 5, 0.6, false, true, 0.6, 0.17857142857142855},
        {20, 20, 0.15, false, true, 0.15, 0.6437041219649916},
        {100, 100, 0.25, false, true, 0.25, 0.001864746185927891},
        {7, 3, 0.5, false, true, 0.47619047619047616, 0.26666666666666666},
        {40, 25, 0.3, false, true, 0.3, 0.05067042227907703},
        // Exact reference. Divides by binom(350, 200), which is only accurate to ~1e-13
        // and the error depends on whether the compiler contracts to FMA.
        {200, 150, 0.1, false, true, 0.1, 0.17015140184061314, 1e-12},
        {1000, 1, 1.0, false, true, 1.0, 0.000999000999000999},
        // The path count overflows, so the exact calculation is unsuccessful.
        {1100, 1099, 8.271982794275788e-07, false, false, 8.271982794275788e-07, qnan},
        {2000, 1000, 0.0005, false, false, 0.0005, qnan},
    };

    for (const auto &c : cases) {
        const int64_t g = std::gcd(c.n1, c.n2);
        auto work = make_work(c.n1, c.n2);
        const auto [success, d, prob] =
            xsf::attempt_exact_2kssamp(c.n1, c.n2, g, c.d, c.two_sided, std::mdspan(work.data(), work.size()));
        const auto rel_error = xsf::extended_relative_error(prob, c.expected_prob);
        CAPTURE(c.n1, c.n2, c.d, c.two_sided, success, d, prob, c.expected_prob, c.rtol, rel_error);
        REQUIRE(success == c.expected_success);
        REQUIRE(d == c.expected_d);
        REQUIRE(rel_error <= c.rtol);
    }
}

TEST_CASE("compute_prob_outside_square test", "[compute_prob_outside_square][xsf_tests]") {
    const double rtol = 1e-14;
    const std::vector<OutsideSquareCase> cases = {
        {1, 1, 1.0},
        {5, 3, 0.35714285714285715},
        {20, 3, 0.9831368772656193},
        {100, 25, 0.0037294923618311367},
        {500, 7, 0.9999999998936827},
    };

    for (const auto &c : cases) {
        const double result = xsf::compute_prob_outside_square(c.n, c.h);
        const auto rel_error = xsf::extended_relative_error(result, c.expected);
        CAPTURE(c.n, c.h, result, c.expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("compute_outer_prob_inside_method test", "[compute_outer_prob_inside_method][xsf_tests]") {
    const double rtol = 1e-14;
    const std::vector<LatticeCase> cases = {
        {1, 1, 1, 1, 1.0},
        {7, 3, 1, 11, 0.5333333333333334},
        {40, 25, 5, 60, 0.1013339946846},
        {200, 150, 50, 60, 0.33868681219914437},
        {150, 200, 50, 60, 0.33868681219914437},
    };

    for (const auto &c : cases) {
        auto work = make_work(c.m, c.n);
        const double result =
            xsf::compute_outer_prob_inside_method(c.m, c.n, c.g, c.h, std::mdspan(work.data(), work.size()));
        const auto rel_error = xsf::extended_relative_error(result, c.expected);
        CAPTURE(c.m, c.n, c.g, c.h, result, c.expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("count_paths_outside_method test", "[count_paths_outside_method][xsf_tests]") {
    const std::vector<LatticeCase> cases = {
        {7, 3, 1, 11, 32.0},
        {40, 25, 5, 60, 3.3021289646854336e+16},
        // Exact reference. The recursion subtracts large binomial coefficients.
        {200, 150, 50, 60, 4.666891420789811e+101, 1e-12},
        {1000, 1, 1, 1001, 1001.0},
        {3, 7, 1, 11, 32.0},
    };

    for (const auto &c : cases) {
        auto work = make_work(c.m, c.n);
        const double result =
            xsf::count_paths_outside_method(c.m, c.n, c.g, c.h, std::mdspan(work.data(), work.size()));
        const auto rel_error = xsf::extended_relative_error(result, c.expected);
        CAPTURE(c.m, c.n, c.g, c.h, result, c.expected, c.rtol, rel_error);
        REQUIRE(rel_error <= c.rtol);
    }

    // The path count overflows; the caller is expected to check for a finite answer.
    auto work = make_work(1100, 1099);
    REQUIRE_FALSE(
        std::isfinite(xsf::count_paths_outside_method(1100, 1099, 1, 1, std::mdspan(work.data(), work.size())))
    );
    work = make_work(2000, 1000);
    REQUIRE_FALSE(
        std::isfinite(xsf::count_paths_outside_method(2000, 1000, 1, 1, std::mdspan(work.data(), work.size())))
    );
}

/*
Cases ported from TestKSTwoSamples in scipy/stats/tests/test_stats.py, which test
stats.ks_2samp on data sets whose sample sizes and statistic select the exact method.
Here the sample sizes and statistic are passed to attempt_exact_2kssamp directly and the
reference p-values are those of the SciPy tests. "less" and "greater" both map to
two_sided = false. A NaN reference means the exact calculation is unsuccessful and
ks_2samp falls back to the asymptotic method.

The tolerances match the SciPy tests: xp_assert_close with its default
rtol = 4 * sqrt(eps) and atol = 1e-14. Note that the references for testLarge
(two-sided) and test_gh11184 are not accurate: the exact p-values are 6.26e-26 and
3.85e-17. They agree only because of the absolute tolerance.
*/
TEST_CASE("attempt_exact_2kssamp scipy ks_2samp cases", "[attempt_exact_2kssamp][xsf_tests]") {
    struct ScipyCase {
        const char *test_name;
        int64_t n1;
        int64_t n2;
        double d;
        bool two_sided;
        double expected_prob;
    };
    const double rtol = 4 * std::sqrt(std::numeric_limits<double>::epsilon());
    const double atol = 1e-14;
    const std::vector<ScipyCase> cases = {
        {"testSmall", 1, 1, 1.0, true, 1.0},
        {"testSmall", 1, 1, 1.0, false, 0.5},
        {"testSmall", 1, 1, 0.0, false, 1.0},
        {"testTwoVsThree", 2, 3, 1.0 / 3, true, 1.0},
        {"testTwoVsThree", 2, 3, 1.0 / 3, false, 0.7},
        {"testTwoVsThree", 2, 3, 2.0 / 3, true, 0.6},
        {"testTwoVsThree", 2, 3, 2.0 / 3, false, 0.3},
        {"testTwoVsThree", 2, 3, 0.0, false, 1.0},
        {"testTwoVsFour", 2, 4, 2.0 / 4, true, 14.0 / 15},
        {"testTwoVsFour", 2, 4, 2.0 / 4, false, 8.0 / 15},
        {"testTwoVsFour", 2, 4, 1.0 / 4, false, 12.0 / 15},
        {"testTwoVsFour", 2, 4, 3.0 / 4, true, 6.0 / 15},
        {"testTwoVsFour", 2, 4, 3.0 / 4, false, 3.0 / 15},
        {"testTwoVsFour", 2, 4, 0.0, false, 1.0},
        {"test100_100", 100, 100, 3.0 / 100, true, 0.999999999996206},
        {"test100_100", 100, 100, 3.0 / 100, false, 0.9143290114276248},
        {"test100_100", 100, 100, 0.0, false, 1.0},
        {"test100_100", 100, 100, 2.0 / 100, true, 1.0},
        {"test100_100", 100, 100, 2.0 / 100, false, 0.960978450786184},
        {"test100_110", 100, 110, 232.0 / 1100, true, 0.015739183865607353},
        {"test100_110", 100, 110, 232.0 / 1100, false, 0.007869594319053203},
        {"test100_110", 100, 110, 0.0, false, 1.0},
        {"test100_110", 100, 110, 229.0 / 1100, true, 0.017803803861026313},
        {"test100_110", 100, 110, 229.0 / 1100, false, 0.008901905958245056},
        {"testRepeatedValues", 16, 16, 5.0 / 16, true, 0.4262934613454952},
        {"testRepeatedValues", 16, 16, 5.0 / 16, false, 0.21465428276573786},
        {"testRepeatedValues", 16, 16, 0.0, false, 1.0},
        {"testRepeatedValues", 21, 26, 190.0 / 21 / 26, true, 0.0919245790168125},
        {"testRepeatedValues", 21, 26, 190.0 / 21 / 26, false, 0.0459633806858544},
        {"testRepeatedValues", 21, 26, 70.0 / 21 / 26, false, 0.6121593130022775},
        {"testEqualSizes", 3, 3, 1.0 / 3, true, 1.0},
        {"testEqualSizes", 3, 3, 1.0 / 3, false, 0.75},
        {"testEqualSizes", 3, 3, 0.0, false, 1.0},
        {"testMiddlingBoth", 500, 600, 2000.0 / 500 / 600, true, 1.0},
        {"testMiddlingBoth", 500, 600, 2000.0 / 500 / 600, false, qnan},
        {"testMiddlingBoth", 500, 600, 500.0 / 500 / 600, false, qnan},
        {"testMediumBoth", 1000, 1100, 6600.0 / 1000 / 1100, true, 1.0},
        {"testMediumBoth", 1000, 1100, 6600.0 / 1000 / 1100, false, qnan},
        {"testMediumBoth", 1000, 1100, 1000.0 / 1000 / 1100, false, qnan},
        {"testLarge", 10000, 110, 55275.0 / 110000, true, 4.2188474935755949e-15},
        {"testLarge", 10000, 110, 561.0 / 110000, false, 0.99115454582047591},
        {"testLarge", 10000, 110, 55275.0 / 110000, false, 3.1317328311518713e-26},
        {"test_gh11184", 3000, 3001, 0.11292880151060758, true, 2.7755575615628914e-15},
        {"testLargeBoth", 10000, 11000, 563.0 / 110000, true, 0.9990456491488628},
        {"testLargeBoth", 10000, 11000, 563.0 / 110000, false, qnan},
        {"testLargeBoth", 10000, 11000, 10.0 / 110000, false, qnan},
    };

    for (const auto &c : cases) {
        const int64_t g = std::gcd(c.n1, c.n2);
        auto work = make_work(c.n1, c.n2);
        const auto [success, d, prob] =
            xsf::attempt_exact_2kssamp(c.n1, c.n2, g, c.d, c.two_sided, std::mdspan(work.data(), work.size()));
        const bool expected_success = !std::isnan(c.expected_prob);
        CAPTURE(c.test_name, c.n1, c.n2, c.d, c.two_sided, success, d, prob, c.expected_prob, rtol, atol);
        REQUIRE(success == expected_success);
        REQUIRE(std::abs(d - c.d) <= rtol * std::abs(c.d));
        if (expected_success) {
            REQUIRE(std::abs(prob - c.expected_prob) <= atol + rtol * std::abs(c.expected_prob));
        } else {
            REQUIRE(std::isnan(prob));
        }
    }
}
