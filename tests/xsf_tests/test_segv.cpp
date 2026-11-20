#include "../testing_utils.h"

#include <xsf/cpu/sphd_wave.h>

TEST_CASE("prolate_segv, oblate_segv negative c", "[prolate_segv][oblate_segv][xsf_tests]") {
    // Regression test: the characteristic values depend only on c^2 (DLMF 30.3),
    // but c < 0 returned n * (n + 1), e.g. prolate_segv(1, 1, -10) was 2.
    // A |c| too large for the arrays in segv gives NaN.

    /*
    Reference values with mpmath: eigenvalues of l (l + 1) + s c^2 x^2 in the
    normalized associated Legendre basis (s = 1 prolate, s = -1 oblate), where
    x Pbar_l^m = a_l Pbar_{l+1}^m + a_{l-1} Pbar_{l-1}^m:

    import mpmath as mp
    mp.mp.dps = 40

    def cv(m, n, c, s, N=60):
        a = lambda l: mp.sqrt(mp.mpf((l + 1) ** 2 - m * m) / ((2 * l + 1) * (2 * l + 3)))
        M = mp.zeros(N, N)
        for i in range(N):
            l = m + i
            M[i, i] = l * (l + 1) + s * c**2 * (a(l) ** 2 + (a(l - 1) ** 2 if i else 0))
            if i + 2 < N:
                M[i, i + 2] = M[i + 2, i] = s * c**2 * a(l) * a(l + 1)
        return float(sorted(mp.eigsy(M, eigvals_only=True))[n - m])

    print(repr(cv(1, 1, 10, 1)), repr(cv(0, 2, 5, 1)), repr(cv(0, 2, 5, -1)))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<bool, double, double, double, double>;
    auto [prolate, m, n, c, desired] = GENERATE(
        test_case{true, 1, 1, -10, 10.287768767391468}, test_case{true, 0, 2, -5, 20.17691472053327},
        test_case{false, 0, 2, -5, -2.448598903320186},
        test_case{true, 0, 0, -1000, std::numeric_limits<double>::quiet_NaN()}
    );

    const double out = prolate ? xsf::prolate_segv(m, n, c) : xsf::oblate_segv(m, n, c);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 5e-13;
    CAPTURE(prolate, m, n, c, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("prolate_aswfa_nocv large c", "[prolate_aswfa_nocv][xsf_tests]") {
    // Regression test: a |c| too large for segv gives NaN, not the c = 0 value.
    double s1f, s1d;
    xsf::prolate_aswfa_nocv(0.0, 0.0, -1000.0, 0.5, s1f, s1d);
    CAPTURE(s1f, s1d);
    REQUIRE(std::isnan(s1f));
    REQUIRE(std::isnan(s1d));
}
