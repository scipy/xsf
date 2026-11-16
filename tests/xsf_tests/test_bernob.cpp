#include "../testing_utils.h"

#include <xsf/cpu/specfun/specfun.h>

TEST_CASE("specfun::bernob small even n", "[bernob][xsf_tests]") {
    // The zeta(n) sum for B_n was cut off too early, giving relative errors of
    // 1.7e-12 for B_4, 6e-14 for B_6 and 1e-14 for B_8.

    /*
    import mpmath

    mpmath.mp.dps = 50

    for n in [4, 6, 8, 10, 30]:
        print(f"B_{n} = {float(mpmath.bernoulli(n))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, double>;
    auto [n, desired] = GENERATE(
        test_case{4, -0.03333333333333333}, test_case{6, 0.023809523809523808}, test_case{8, -0.03333333333333333},
        test_case{10, 0.07575757575757576}, test_case{30, 601580873.9006424}
    );

    double bn[31];
    xsf::specfun::bernob(30, bn);
    const double error = xsf::extended_relative_error(bn[n], desired);
    const double rtol = 2e-15;
    CAPTURE(n, bn[n], desired, error, rtol);
    REQUIRE(error <= rtol);
}
