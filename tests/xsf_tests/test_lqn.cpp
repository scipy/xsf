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

TEST_CASE("lqn lqmn Q_0 for |x| << 1 and |x| >> 1", "[lqn][lqmn][real][xsf_tests]") {
    // Regression test: log(|(1 + x) / (1 - x)|) / 2 lost digits for |x| << 1 and |x| >> 1,
    // e.g. relative error 8.3e-8 for lqn(0, 1e-10).

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for x, kind in [(1e-10, 2), (-1e-6, 2), (1e6, 3)]:
        q = mpmath.legenq(0, 0, mpmath.mpf(x), type=kind).real
        print(f"x={x!r}: Q_0={float(q)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, bool>;
    auto [x, desired, use_lqn] = GENERATE(
        test_case{1e-10, 1e-10, true}, test_case{-1e-6, -1.0000000000003333e-06, true},
        test_case{1e6, 1.0000000000003333e-06, false}
    );

    // lqmn needs space for at least 2x2 values
    std::vector<double> q(4), d(4);
    if (use_lqn) {
        xsf::lqn(x, xsf::cxx::mdspan(q.data(), 2), xsf::cxx::mdspan(d.data(), 2));
    } else {
        xsf::lqmn(x, xsf::cxx::mdspan(q.data(), 2, 2), xsf::cxx::mdspan(d.data(), 2, 2));
    }

    const double error = xsf::extended_relative_error(q[0], desired);
    const double rtol = 5e-16;
    CAPTURE(x, q[0], desired, error, rtol);
    REQUIRE(error <= rtol);
}
