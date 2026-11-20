#include "../testing_utils.h"

#include <xsf/cpu/exp.h>

TEST_CASE("expm1 complex where exp(Re(z)) overflows", "[expm1][xsf_tests]") {
    // Regression test for complex expm1 where e^Re(z) overflows (gh-1730 in
    // jorenham/xsf-rust): (expm1(x) + 1) * sin(y) was inf*0 = nan or inf*tiny = inf.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    zs = [(709.9, 1e-300), (800, 1e-200), (709.9, 1.5707963267948966), (1420, -5e-324)]
    for zr, zi in zs:
        r = mpmath.expm1(mpmath.mpc(zr, zi))
        print(f"expm1({zr!r}, {zi!r}) = ({float(r.real)!r}, {float(r.imag)!r})")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [zr, zi, desired_re, desired_im] = GENERATE(
        test_case{709.9, 1e-300, std::numeric_limits<double>::infinity(), 202140205.6119564},
        test_case{800.0, 1e-200, std::numeric_limits<double>::infinity(), 2.7263745721125666e+147},
        test_case{709.9, 1.5707963267948966, 1.2377517789083512e+292, std::numeric_limits<double>::infinity()},
        test_case{1420.0, -5e-324, std::numeric_limits<double>::infinity(), -2.4657495327684952e+293}
    );

    const std::complex<double> out = xsf::expm1(std::complex<double>{zr, zi});
    const double rtol = 5e-15;
    CAPTURE(zr, zi, out, desired_re, desired_im, rtol);
    for (auto [o, d] : {std::pair{out.real(), desired_re}, std::pair{out.imag(), desired_im}}) {
        if (std::isinf(d)) {
            REQUIRE(o == d);
        } else {
            REQUIRE(xsf::extended_relative_error(o, d) <= rtol);
        }
    }
}
