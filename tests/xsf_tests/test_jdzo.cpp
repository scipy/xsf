#include "../testing_utils.h"

#include <vector>

#include <xsf/cpu/specfun/specfun.h>

TEST_CASE("specfun::jdzo missed zeros of Jn'", "[jdzo][xsf_tests]") {
    // Regression test: jdzo did not search high enough orders (nt = 210, 831), and missed
    // zeros of Jn' whose initial guess was above xm (nt = 1121), so the nt-th zero was wrong.

    /*
    import mpmath as mp

    mp.mp.dps = 40
    for n, m in [(26, 1), (54, 1), (56, 2)]:
        print(n, m, repr(float(mp.besseljzero(n, m, derivative=1))))
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<int, int, int, double>;
    auto [nt, desired_n, desired_m, desired] = GENERATE(
        test_case{210, 26, 1, 28.418074834560365}, test_case{831, 54, 1, 57.07462600766977},
        test_case{1121, 56, 2, 66.37255694220194}
    );

    // jdzo can compute up to nt + 24 zeros, and zo is 1-based
    std::vector<double> zo(nt + 32);
    std::vector<int> n(nt + 32), m(nt + 32), p(nt + 32);
    xsf::specfun::jdzo(nt, zo.data(), n.data(), m.data(), p.data());

    const double out = zo[nt];
    const double error = xsf::extended_relative_error(out, desired);
    const double rtol = 1e-14;
    CAPTURE(nt, n[nt - 1], m[nt - 1], p[nt - 1], out, desired, error, rtol);
    REQUIRE(n[nt - 1] == desired_n);
    REQUIRE(m[nt - 1] == desired_m);
    REQUIRE(p[nt - 1] == 1);
    REQUIRE(error <= rtol);
}
