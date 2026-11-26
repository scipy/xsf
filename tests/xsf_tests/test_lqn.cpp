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

TEST_CASE("lqn lqmn Q_0 for complex z with |z| << 1 or tiny Re z", "[lqn][lqmn][complex][xsf_tests]") {
    // Regression test: log((1 + z) / (1 - z)) / 2 lost the digits of the real part for |z| << 1,
    // e.g. relative error 8.6e-8 for lqn(0, 1e-9+2e-9j), and for tiny Re z with 1/2 <= |z| <= 2,
    // e.g. relative error 1.9e-6 for lqn(0, 1e-10+1.5j). Also |z| >> 1.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 700

    for z, kind in [(1e-9 + 2e-9j, 2), (1e-300 + 1.5j, 3), (1e6 + 3e5j, 3)]:
        q = mpmath.legenq(0, 0, mpmath.mpc(z), type=kind)
        print(f"z={z!r}: Q_0=({float(q.real)!r}, {float(q.imag)!r})")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<std::complex<double>, std::complex<double>, bool>;
    auto [z, desired, use_lqn] = GENERATE(
        test_case{{1e-9, 2e-9}, {1e-09, 2e-09}, true},
        test_case{{1e-10, 1.5}, {3.076923076923077e-11, -0.5880026035475675}, true},
        test_case{{1e-300, 1.5}, {3.076923076923077e-301, -0.5880026035475675}, true},
        test_case{{1e-300, 1.5}, {3.076923076923077e-301, -0.5880026035475675}, false},
        test_case{{1e6, 3e5}, {9.174311926607384e-07, -2.7522935779838986e-07}, true}
    );

    // lqmn needs space for at least 2x2 values
    std::vector<std::complex<double>> q(4), d(4);
    if (use_lqn) {
        xsf::lqn(z, xsf::cxx::mdspan(q.data(), 2), xsf::cxx::mdspan(d.data(), 2));
    } else {
        xsf::lqmn(z, xsf::cxx::mdspan(q.data(), 2, 2), xsf::cxx::mdspan(d.data(), 2, 2));
    }

    const double error_re = xsf::extended_relative_error(q[0].real(), desired.real());
    const double error_im = xsf::extended_relative_error(q[0].imag(), desired.imag());
    const double rtol = 5e-16;
    CAPTURE(z, q[0], desired, error_re, error_im, rtol);
    REQUIRE(error_re <= rtol);
    REQUIRE(error_im <= rtol);
}

TEST_CASE("lqn |x| just above 1.021", "[lqn][real][xsf_tests]") {
    // Regression test for scipy/xsf#552 (comment 7052780497): the series for the start values
    // of the backward recurrence was cut off after 500 terms, e.g. relative error 6.4e-5 for
    // Q_0(-1.02101) with n = 100.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    # Q_n(-x) = (-1)^(n+1) Q_n(x) for x > 1
    for x, k in [(-1.02101, 0), (1.02101, 0), (-1.02101, 100)]:
        s = (-1) ** (k + 1) if x < 0 else 1
        q = s * mpmath.legenq(k, 0, abs(mpmath.mpf(x)), type=3).real
        print(f"x={x!r}, k={k}: Q_k={float(q)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, int, int, double>;
    auto [x, n, k, desired] = GENERATE(
        test_case{-1.02101, 100, 0, -2.2831770752899967}, test_case{1.02101, 200, 0, 2.2831770752899967},
        test_case{-1.02101, 100, 100, -3.20577818739551e-10}
    );

    std::vector<double> q(n + 1), d(n + 1);
    xsf::lqn(x, xsf::cxx::mdspan(q.data(), q.size()), xsf::cxx::mdspan(d.data(), d.size()));

    const double error = xsf::extended_relative_error(q[k], desired);
    const double rtol = 1e-13;
    CAPTURE(x, n, k, q[k], desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("lqn complex z with Re z = 1", "[lqn][complex][xsf_tests]") {
    // Regression test: lqn returned the 1e300 pole sentinel on the whole line Re z = 1,
    // and without it, the backward recurrence started at km = 0 for z = 1 + 0.4j.
    // Q_0'(z) = 1 / (1 - z^2) = 0.2403846153846154 + 1.2019230769230769j
    SET_FP_FORMAT()
    const std::complex<double> z(1, 0.4);
    std::vector<std::complex<double>> q(2), d(2);
    xsf::lqn(z, xsf::cxx::mdspan(q.data(), q.size()), xsf::cxx::mdspan(d.data(), d.size()));

    const double rtol = 1e-15;
    const double error_re = xsf::extended_relative_error(d[0].real(), 0.2403846153846154);
    const double error_im = xsf::extended_relative_error(d[0].imag(), 1.2019230769230769);
    CAPTURE(d[0], error_re, error_im, rtol);
    REQUIRE(error_re <= rtol);
    REQUIRE(error_im <= rtol);
}
