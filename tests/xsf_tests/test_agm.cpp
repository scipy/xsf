#include "../testing_utils.h"

#include <xsf/agm.h>

TEST_CASE("agm simple", "[agm][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4444-L4455
    const double rtol = 1e-13;

    const auto gauss_constant = 1 / xsf::agm(1, std::sqrt(2));
    const auto gauss_constant_expected = 0.834626841674073186;
    const auto gauss_constant_rel_error = xsf::extended_relative_error(gauss_constant, gauss_constant_expected);
    CAPTURE(gauss_constant, gauss_constant_expected, rtol, gauss_constant_rel_error);
    REQUIRE(gauss_constant_rel_error <= rtol);

    using test_case = std::tuple<double, double, double>;
    auto [a_wolfram, b_wolfram, expected_wolfram] = GENERATE(
        test_case{1, 1, 1}, test_case{1, 3, 1.863616783244897}, test_case{1, 5, 2.604008190530940},
        test_case{3, 1, 1.863616783244897}, test_case{3, 3, 3}, test_case{3, 5, 3.936235503649555}
    );

    const auto output_wolfram = xsf::agm(a_wolfram, b_wolfram);
    const auto rel_error_wolfram = xsf::extended_relative_error(output_wolfram, expected_wolfram);
    CAPTURE(a_wolfram, b_wolfram, output_wolfram, expected_wolfram, rtol, rel_error_wolfram);
    REQUIRE(rel_error_wolfram <= rtol);
}

TEST_CASE("agm finite cases", "[agm][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4457-L4480
    const double rtol = 1e-13;
    using test_case = std::tuple<double, double, double>;
    auto [a, b, expected] = GENERATE(
        test_case{1, 2, 1.4567910310469068}, test_case{2, 1, 1.4567910310469068},
        test_case{-1, -2, -1.4567910310469068},
        test_case{24, 6, 13.458171481725614}, test_case{13, 123456789.5, 11111458.498599306},
        test_case{1e30, 1, 2.229223055945383e+28}, test_case{1e-22, 1, 0.030182566420169886},
        test_case{1e150, 1e180, 2.229223055945383e+178},
        test_case{1e180, 1e-150, 2.0634722510162677e+177},
        test_case{1e-150, 1e-170, 3.3112619670463756e-152},
        test_case{std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), 1.9892072050015473e+305},
        test_case{
            0.75 * std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 1.564904312298045e+308
        },
        test_case{
            std::numeric_limits<double>::min(), 3 * std::numeric_limits<double>::min(), 4.1466849866735005e-308
        }
    );

    const auto output = xsf::agm(a, b);
    const auto rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(a, b, output, expected, rtol, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("agm zero, nan, and inf cases", "[agm][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.0rc1/scipy/special/tests/test_basic.py#L4482-L4499
    REQUIRE(xsf::agm(0.0, 0.0) == 0);
    REQUIRE(xsf::agm(99.0, 0.0) == 0);

    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    REQUIRE(std::isnan(xsf::agm(-1.0, 10.0)));
    REQUIRE(std::isnan(xsf::agm(0, inf)));
    REQUIRE(std::isnan(xsf::agm(inf, 0)));
    REQUIRE(std::isnan(xsf::agm(0, -inf)));
    REQUIRE(std::isnan(xsf::agm(-inf, 0)));
    REQUIRE(std::isnan(xsf::agm(inf, -inf)));
    REQUIRE(std::isnan(xsf::agm(-inf, inf)));
    REQUIRE(std::isnan(xsf::agm(1, nan)));
    REQUIRE(std::isnan(xsf::agm(nan, -1)));

    REQUIRE(xsf::agm(1, inf) == inf);
    REQUIRE(xsf::agm(inf, 1) == inf);
    REQUIRE(xsf::agm(-1, -inf) == -inf);
    REQUIRE(xsf::agm(-inf, -1) == -inf);
}
