#include "../testing_utils.h"
#include <tuple>
#include <xsf/cpu/airy.h>

TEST_CASE("itairy", "[itairy][xsf_tests]") {
    using test_case = std::tuple<double, double, double, double, double, double>;
    // Reference values were computed with the Python library mpmath.
    auto [x, ref_apt, ref_bpt, ref_ant, ref_bnt, rtol] = GENERATE(
        test_case{1.0, 0.2363173419171098, 0.8727691167380082, 0.4656739834670686, 0.37300500963429495, 1e-13},
        test_case{4.0, 0.3328926453861221, 48.21947526380343, 0.47800749642926166, -0.05076005849528754, 1e-12},
        test_case{6.0, 0.33332945170523853, 2838.369692679631, 0.7326752562674856, -0.1303810644524539, 1e-11},
        test_case{8.5, 0.33333332966571266, 1760815.3976228302, 0.6673920742577031, -0.11303861602723536, 1e-7},
        test_case{9.0, 0.33333333253066366, 7377921.170522001, 0.5588197489447188, -0.010330373551680197, 1e-7}
    );
    double apt, bpt, ant, bnt;
    xsf::itairy(x, apt, bpt, ant, bnt);

    const auto rel_error_apt = xsf::extended_relative_error(apt, ref_apt);
    const auto rel_error_bpt = xsf::extended_relative_error(bpt, ref_bpt);
    const auto rel_error_ant = xsf::extended_relative_error(ant, ref_ant);
    const auto rel_error_bnt = xsf::extended_relative_error(bnt, ref_bnt);

    CAPTURE(x, apt, ref_apt, rel_error_apt, bpt, ref_bpt, rel_error_bpt, rtol);
    CAPTURE(ant, ref_ant, rel_error_ant, bnt, ref_bnt, rel_error_bnt);
    REQUIRE(rel_error_apt <= rtol);
    REQUIRE(rel_error_bpt <= rtol);
    REQUIRE(rel_error_ant <= rtol);
    REQUIRE(rel_error_bnt <= rtol);
}

TEST_CASE("itairy overflow", "[itairy][xsf_tests]") {
    // Regression test for gh-448 (https://github.com/scipy/xsf/issues/448):
    // bpt overflowed for 104.2 < x <= 104.66, and x > 5.9e204 gave NaNs instead of the limits.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for x in [104.5]:
        x = mpmath.mpf(x)
        apt = mpmath.airyai(x, derivative=-1)
        bpt = mpmath.airybi(x, derivative=-1)
        ant = -mpmath.airyai(-x, derivative=-1)
        bnt = -mpmath.airybi(-x, derivative=-1)
        print(f"{float(x)!r}, {float(apt)!r}, {float(bpt)!r}, {float(ant)!r}, {float(bnt)!r}")
    */
    SET_FP_FORMAT()
    const double inf = std::numeric_limits<double>::infinity();
    using test_case = std::tuple<double, double, double, double, double>;
    auto [x, ref_apt, ref_bpt, ref_ant, ref_bnt] = GENERATE_COPY(
        test_case{104.5, 0.3333333333333333, 3.3773949845225984e+307, 0.683625159831577, 0.0032219271426297096},
        // the limits, DLMF 9.10.11, 9.10.12
        test_case{1e205, 1.0 / 3.0, inf, 2.0 / 3.0, 0.0}, test_case{inf, 1.0 / 3.0, inf, 2.0 / 3.0, 0.0}
    );
    double apt, bpt, ant, bnt;
    xsf::itairy(x, apt, bpt, ant, bnt);

    const double rtol = 1e-12;
    CAPTURE(x, apt, ref_apt, bpt, ref_bpt, ant, ref_ant, bnt, ref_bnt, rtol);
    REQUIRE(xsf::extended_relative_error(apt, ref_apt) <= rtol);
    REQUIRE(xsf::extended_relative_error(bpt, ref_bpt) <= rtol);
    REQUIRE(xsf::extended_relative_error(ant, ref_ant) <= rtol);
    REQUIRE(xsf::extended_relative_error(bnt, ref_bnt) <= rtol);
}
