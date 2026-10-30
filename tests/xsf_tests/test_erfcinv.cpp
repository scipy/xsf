#include "../testing_utils.h"

#include <xsf/cephes/erfinv.h>

TEST_CASE("cephes::erfcinv subnormal arguments", "[cephes::erfcinv][xsf_tests]") {
    // Regression test for gh-309 (https://github.com/scipy/xsf/issues/309):
    // 0.5 * y lost precision (or underflowed to 0) for subnormal y, e.g.
    // erfcinv(5e-324) was inf and erfcinv(1.5e-323) was 27.187827569944027.
    // DBL_MIN = 2**-1022 and 2 * DBL_MIN = 2**-1021 check either side of the
    // y < 2 * DBL_MIN branch.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for y in [5e-324, 1e-323, 1.5e-323, 2**-1022, 2**-1021]:
        # erfc(x) underflows in double precision here, so solve log(erfc(x)) = log(y).
        x = mpmath.findroot(lambda x: mpmath.log(mpmath.erfc(x)) - mpmath.log(y), 27)
        print(f"erfcinv({y!r}) = {float(x)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [y, desired] = GENERATE(
        test_case{5e-324, 27.21329321081295}, test_case{1e-323, 27.200563366536258},
        test_case{1.5e-323, 27.19311412620397}, test_case{std::numeric_limits<double>::min(), 26.54325845425098},
        test_case{2 * std::numeric_limits<double>::min(), 26.530207562276377}
    );

    const double out = xsf::cephes::erfcinv(y);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-16;
    CAPTURE(y, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
