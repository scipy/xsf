#include "../testing_utils.h"
#include <xsf/cosine.h>

namespace {

inline constexpr double M_PI64 = 0x1.921fb54442d18p+1;

}

TEST_CASE("cosine_cdf exact cases", "[cosine_cdf][xsf_tests]") {
    // Mirrors scipy/special/tests/test_cosine_distr.py::test_cosine_cdf_exact
    using test_case = std::tuple<double, double>;
    auto [x, expected] =
        GENERATE(test_case{-4.0, 0.0}, test_case{0.0, 0.5}, test_case{M_PI64, 1.0}, test_case{4.0, 1.0});

    double output = xsf::cosine_cdf(x);
    CAPTURE(x, output, expected);
    REQUIRE(output == expected);
}

TEST_CASE("cosine_cdf close cases", "[cosine_cdf][xsf_tests]") {
    // Mirrors scipy/special/tests/test_cosine_distr.py::test_cosine_cdf
    constexpr double rtol = 5e-15;
    using test_case = std::tuple<double, double>;
    auto [x, expected] = GENERATE(
        test_case{3.1409, 0.999999999991185}, test_case{2.25, 0.9819328173287907},
        test_case{-1.599, 0.08641959838382553}, test_case{-1.601, 0.086110582992713},
        test_case{-2.0, 0.0369709335961611}, test_case{-3.0, 7.522387241801384e-05},
        test_case{-3.1415, 2.109869685443648e-14}, test_case{-3.14159, 4.956444476505336e-19},
        test_case{-M_PI64, 4.871934450264861e-50}
    );

    double output = xsf::cosine_cdf(x);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(x, output, expected, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("cosine_invcdf exact cases", "[cosine_invcdf][xsf_tests]") {
    // Mirrors scipy/special/tests/test_cosine_distr.py::test_cosine_invcdf_exact
    using test_case = std::tuple<double, double>;
    auto [p, expected] = GENERATE(test_case{0.0, -M_PI64}, test_case{0.5, 0.0}, test_case{1.0, M_PI64});

    double output = xsf::cosine_invcdf(p);
    CAPTURE(p, output, expected);
    REQUIRE(output == expected);
}

TEST_CASE("cosine_invcdf invalid p", "[cosine_invcdf][xsf_tests]") {
    // Mirrors scipy/special/tests/test_cosine_distr.py::test_cosine_invcdf_invalid_p
    REQUIRE(std::isnan(xsf::cosine_invcdf(-0.1)));
    REQUIRE(std::isnan(xsf::cosine_invcdf(1.1)));
}

TEST_CASE("cosine_invcdf close cases", "[cosine_invcdf][xsf_tests]") {
    // Mirrors scipy/special/tests/test_cosine_distr.py::test_cosine_invcdf
    constexpr double rtol = 1e-14;
    using test_case = std::tuple<double, double>;
    auto [p, expected] = GENERATE(
        test_case{1e-50, -M_PI64}, test_case{1e-14, -3.1415204137058454}, test_case{1e-08, -3.1343686589124524},
        test_case{0.0018001, -2.732563923138336}, test_case{0.010, -2.41276589008678},
        test_case{0.060, -1.7881244975330157}, test_case{0.125, -1.3752523669869274},
        test_case{0.250, -0.831711193579736}, test_case{0.400, -0.3167954512395289},
        test_case{0.419, -0.25586025626919906}, test_case{0.421, -0.24947570750445663},
        test_case{0.750, 0.831711193579736}, test_case{0.940, 1.7881244975330153},
        test_case{0.9999999996, 3.1391220839917167}
    );

    double output = xsf::cosine_invcdf(p);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(p, output, expected, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}
