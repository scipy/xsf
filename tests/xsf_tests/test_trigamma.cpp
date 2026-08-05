#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/trigamma.h>


TEST_CASE("trigamma real", "[trigamma][xsf_tests]") {
    using test_case = std::tuple<double, double, double>;

    // Reference values from https://github.com/JuliaMath/SpecialFunctions.jl/blob/c9ff36085c9008b9bdaa457e745f459a1827e52e/test/gamma.jl#L30-L44.
    auto [z, ref_trigamma, rtol] = GENERATE(
        test_case{0.1, 101.433299150792758817, 1e-14},
        test_case{0.5, M_PI*M_PI/2, 1e-14},
        test_case{1.0, M_PI*M_PI/6, 1e-14},
        test_case{2.0, M_PI*M_PI/6 - 1, 1e-14},
        test_case{3.0, M_PI*M_PI/6 - 5.0/4, 1e-14},
        test_case{4.0, M_PI*M_PI/6 - 49.0/36, 1e-14},
        test_case{5.0, M_PI*M_PI/6 - 205.0/144, 1e-14},
        test_case{10.0, 0.10516633568168565, 1e-14}
    );

    double result_trigamma = xsf::trigamma(z);
    double rel_err_trigamma = xsf::extended_relative_error(result_trigamma, ref_trigamma);

    CAPTURE(z, result_trigamma, ref_trigamma, rel_err_trigamma, rtol);
    REQUIRE(rel_err_trigamma <= rtol);
}

TEST_CASE("trigamma real poles", "[trigamma][xsf_tests]") {
    auto z = GENERATE(-1.0, 0.0);
    double result_trigamma = xsf::trigamma(z);
    CAPTURE(z, result_trigamma);
    REQUIRE(std::isinf(result_trigamma));
}