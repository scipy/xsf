#include "../testing_utils.h"

#include <xsf/gamma.h>
#include <xsf/loggamma.h>

TEST_CASE("loggamma, gammaln near the zeros at 1 and 2", "[loggamma][gammaln][xsf_tests]") {
    // Regression test for gh-529: the relative error was about eps / |x - x0|.
    /*
    import mpmath

    mpmath.mp.dps = 60

    for x in [2 - 1e-10, 1 + 1e-9, 1 - 1e-9]:
        print(repr(x), repr(float(mpmath.loggamma(mpmath.mpf(x)))))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double>;
    auto [x, desired] = GENERATE(
        test_case{1.9999999999, -4.2278437004755317e-11}, test_case{1.000000001, -5.772157118381039e-10},
        test_case{0.999999999, 5.77215649399226e-10}
    );

    const double rtol = 1e-15;
    const double out_loggamma = xsf::loggamma(x);
    const double out_gammaln = xsf::gammaln(x);
    const double error_loggamma = xsf::extended_relative_error(out_loggamma, desired);
    const double error_gammaln = xsf::extended_relative_error(out_gammaln, desired);
    CAPTURE(x, out_loggamma, out_gammaln, desired, error_loggamma, error_gammaln, rtol);
    REQUIRE(error_loggamma <= rtol);
    REQUIRE(error_gammaln <= rtol);
}
