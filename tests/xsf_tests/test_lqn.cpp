#include "../testing_utils.h"

#include <xsf/cpu/legendre.h>
#include <xsf/mdspan.h>

TEST_CASE("lqn x < -1", "[lqn][real][xsf_tests]") {
    // Regression test for lqn at x < -1: the forward recurrence is unstable
    // there, e.g. Q_30(-951.41) was 9.6e79 instead of -6.97e-103.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    # Q_n(-x) = (-1)^(n+1) Q_n(x) and Q_n'(-x) = (-1)^n Q_n'(x) for x > 1
    for x, n in [(951.41, 30), (1.5, 30), (10.0, 10)]:
        X = mpmath.mpf(x)
        q = mpmath.legenq(n, 0, X, type=3).real
        d = mpmath.diff(lambda t: mpmath.legenq(n, 0, t, type=3).real, X)
        s = (-1) ** (n + 1)
        print(f"x={-x!r}, n={n}: Q={float(s*q)!r} dQ={float(-s*d)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, int, double, double>;
    auto [x, n, desired_q, desired_d] = GENERATE(
        test_case{-951.41, 30, -6.971400693327489e-103, -2.2715080014285592e-104},
        test_case{-1.5, 30, -3.811280512308448e-14, -1.0624790205129686e-12},
        test_case{-10.0, 10, -2.7163939560017536e-15, -3.00374423063823e-15}
    );

    std::vector<double> q(n + 1), d(n + 1);
    xsf::lqn(x, xsf::cxx::mdspan(q.data(), q.size()), xsf::cxx::mdspan(d.data(), d.size()));

    const double error_q = xsf::extended_relative_error(q[n], desired_q);
    const double error_d = xsf::extended_relative_error(d[n], desired_d);
    const double rtol = 1e-13;
    CAPTURE(x, n, q[n], desired_q, error_q, d[n], desired_d, error_d, rtol);
    REQUIRE(error_q <= rtol);
    REQUIRE(error_d <= rtol);
}
