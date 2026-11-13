#include "../testing_utils.h"
#include <tuple>
#include <xsf/cpu/struve.h>

TEST_CASE("itstruve0 large x", "[itstruve0][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;
    // Reference values were computed with the Python library mpmath.
    auto [x, ref, rtol] = GENERATE(
        test_case{31.0, 2.946198073862365, 1e-13}, test_case{40.0, 3.148417627447186, 1e-13},
        test_case{100.0, 3.720914252854685, 1e-13}, test_case{300.0, 4.47323362427753, 1e-13}
    );
    const auto out = xsf::itstruve0(x);
    const auto rel_error = xsf::extended_relative_error(out, ref);
    CAPTURE(x, out, ref, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("itstruve0 moderate x", "[itstruve0][xsf_tests]") {
    // Regression test for gh-399 (https://github.com/scipy/xsf/issues/399):
    // the power series cancelled for moderate x.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 30

    for x in [20.0, 29.75]:
        y = mpmath.quad(lambda t: mpmath.struveh(0, t), mpmath.linspace(0, x, 31))
        print(f"itstruve0({x!r}) = {float(y)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(test_case{20.0, 2.5484516922939573}, test_case{29.75, 3.0836546700490843});

    const double out = xsf::itstruve0(x);
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 2e-15;
    CAPTURE(x, out, desired, error, rtol);
    REQUIRE(error <= rtol);
}
