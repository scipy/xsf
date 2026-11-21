#include "../testing_utils.h"

#include <algorithm>
#include <complex>
#include <vector>

#include <xsf/cpu/legendre.h>
#include <xsf/mdspan.h>

TEST_CASE("lqmn large |x|", "[lqmn][xsf_tests]") {
    // Regression test: the backward recurrences overflowed, giving NaN.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for x, m, n in [(1e7, 0, 0), (1e8, 0, 0), (907.57, 26, 30)]:
        print(x, m, n, repr(float(mpmath.legenq(n, m, mpmath.mpf(x), type=3).real)))
    z = mpmath.mpc(1e8, 3e7)
    q = mpmath.legenq(3, 2, z, type=3)
    print(repr(float(q.real)), repr(float(q.imag)))
    */
    SET_FP_FORMAT()
    const double rtol = 1e-13;

    SECTION("real") {
        using test_case = std::tuple<double, int, int, double>;
        auto [x, m, n, desired] = GENERATE(
            test_case{1e7, 0, 0, 1.0000000000000034e-07}, test_case{1e8, 0, 0, 1e-08},
            test_case{907.57, 26, 30, 8.065789716549577e-60}
        );

        std::vector<double> q((std::max(m, 1) + 1) * (n + 1)), d((std::max(m, 1) + 1) * (n + 1));
        xsf::lqmn(
            x, xsf::cxx::mdspan(q.data(), std::max(m, 1) + 1, n + 1),
            xsf::cxx::mdspan(d.data(), std::max(m, 1) + 1, n + 1)
        );

        const double error = xsf::extended_relative_error(q[m * (n + 1) + n], desired);
        CAPTURE(x, m, n, q[m * (n + 1) + n], desired, error, rtol);
        REQUIRE(error <= rtol);
    }

    SECTION("complex") {
        // looser, since Q_0 loses digits to cancellation in the logarithm for large |z|
        const double rtol = 1e-7;
        const int m = 2, n = 3;
        const std::complex<double> z(1e8, 3e7), desired(3.78987247199564e-33, -8.841146634093655e-33);

        std::vector<std::complex<double>> q((m + 1) * (n + 1)), d((m + 1) * (n + 1));
        xsf::lqmn(z, xsf::cxx::mdspan(q.data(), m + 1, n + 1), xsf::cxx::mdspan(d.data(), m + 1, n + 1));

        const double error = xsf::extended_relative_error(q[m * (n + 1) + n], desired);
        CAPTURE(q[m * (n + 1) + n], desired, error, rtol);
        REQUIRE(error <= rtol);
    }
}

TEST_CASE("lqmn tiny Q_n^0 for large |x|", "[lqmn][xsf_tests]") {
    // Regression test: Q_n^0 and Q_n^1, the start values of the recurrence in m, underflowed.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for x, m, n in [(1e8, 40, 40), (1e8, 6, 37), (2, 175, 565)]:
        print(x, m, n, repr(float(mpmath.legenq(n, m, mpmath.mpf(x), type=3).real)))
    q = mpmath.legenq(40, 40, mpmath.mpc(1e8, 1e7), type=3)
    print(repr(float(q.real)), repr(float(q.imag)))
    */
    SET_FP_FORMAT()
    const double rtol = 1e-13;

    SECTION("real") {
        using test_case = std::tuple<double, int, int, double>;
        auto [x, m, n, desired] = GENERATE(
            test_case{1e8, 40, 40, 1.107541161989153e-270}, test_case{1e8, 6, 37, 4.606612264703314e-307},
            test_case{2.0, 175, 565, -2.2702708407624725e+169}
        );

        std::vector<double> q((m + 1) * (n + 1)), d((m + 1) * (n + 1));
        xsf::lqmn(x, xsf::cxx::mdspan(q.data(), m + 1, n + 1), xsf::cxx::mdspan(d.data(), m + 1, n + 1));

        const double error = xsf::extended_relative_error(q[m * (n + 1) + n], desired);
        CAPTURE(x, m, n, q[m * (n + 1) + n], desired, error, rtol);
        REQUIRE(error <= rtol);
    }

    SECTION("complex") {
        const int m = 40, n = 40;
        const std::complex<double> z(1e8, 1e7), desired(-5.291583069659328e-271, 7.319261139487894e-271);

        std::vector<std::complex<double>> q((m + 1) * (n + 1)), d((m + 1) * (n + 1));
        xsf::lqmn(z, xsf::cxx::mdspan(q.data(), m + 1, n + 1), xsf::cxx::mdspan(d.data(), m + 1, n + 1));

        const double error = xsf::extended_relative_error(q[m * (n + 1) + n], desired);
        CAPTURE(q[m * (n + 1) + n], desired, error, rtol);
        REQUIRE(error <= rtol);
    }
}
