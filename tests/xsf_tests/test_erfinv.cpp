#include "../testing_utils.h"

#include <xsf/cephes/erfinv.h>

TEST_CASE("cephes::erfinv central and tail accuracy", "[cephes::erfinv][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    // Reference values computed with mpmath at 120 digits.
    auto [y, expected] = GENERATE(
        test_case{9.999999e-8, 8.862268368300678e-8}, test_case{1e-7, 8.862269254527602e-8},
        test_case{1.0000001e-7, 8.86227014075453e-8}, test_case{1e-8, 8.86226925452758e-9},
        test_case{0.49999999, 0.47693626507862175}, test_case{0.5, 0.4769362762044699},
        test_case{0.75, 0.8134198475976185}, test_case{0.9999999999999999, 5.8635847487551676}
    );
    const double sign = GENERATE(-1.0, 1.0);
    const double actual = xsf::cephes::erfinv(sign * y);
    CAPTURE(y, sign, actual, expected);
    CHECK(xsf::extended_relative_error(actual, sign * expected) <= 1e-15);
}

TEST_CASE("cephes::erfinv domain and origin", "[cephes::erfinv][xsf_tests]") {
    const double sign = GENERATE(-1.0, 1.0);
    CHECK(xsf::cephes::erfinv(sign) == sign * std::numeric_limits<double>::infinity());
    CHECK(std::isnan(xsf::cephes::erfinv(sign * 1.01)));
    CHECK(std::isnan(xsf::cephes::erfinv(std::numeric_limits<double>::quiet_NaN())));
    const double zero = xsf::cephes::erfinv(sign * 0.0);
    CHECK(zero == 0.0);
    CHECK(std::signbit(zero) == (sign < 0));
    const double tiny = sign * std::numeric_limits<double>::denorm_min();
    CHECK(xsf::cephes::erfinv(tiny) == tiny);
}
