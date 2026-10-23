#include "../testing_utils.h"

#include <xsf/cephes/erfinv.h>

TEST_CASE("cephes::erfcinv subnormal arguments", "[cephes::erfcinv][xsf_tests]") {
    // 0.5 * y used to lose precision (or underflow to 0) for subnormal y,
    // e.g. erfcinv(5e-324) was inf and erfcinv(1.5e-323) was 27.187827569944027.
    // The expected values were calculated with mpmath:
    //   mpmath.mp.dps = 60
    //   mpmath.findroot(lambda x: mpmath.log(mpmath.erfc(x)) - mpmath.log(y), 27)
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [y, desired] = GENERATE(test_case{5e-324, 27.213293210812949}, test_case{1.5e-323, 27.193114126203970});

    const double out = xsf::cephes::erfcinv(y);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-16;
    CAPTURE(y, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
