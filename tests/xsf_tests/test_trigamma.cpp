#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/gamma.h>
#include <xsf/trigamma.h>
#include <xsf/zeta.h>

TEST_CASE("trigamma real", "[trigamma][xsf_tests]") {
    double rtol = 1e-14;
    using test_case = std::tuple<double, double>;

    auto [z, ref_trigamma] = GENERATE(
        // Reference values from
        // https://github.com/JuliaMath/SpecialFunctions.jl/blob/c9ff36085c9008b9bdaa457e745f459a1827e52e/test/gamma.jl#L30-L44.
        test_case{0.1, 101.433299150792758817}, test_case{0.5, M_PI * M_PI / 2}, test_case{1.0, M_PI * M_PI / 6},
        test_case{2.0, M_PI * M_PI / 6 - 1}, test_case{3.0, M_PI * M_PI / 6 - 5.0 / 4},
        test_case{4.0, M_PI * M_PI / 6 - 49.0 / 36}, test_case{5.0, M_PI * M_PI / 6 - 205.0 / 144},
        test_case{10.0, 0.10516633568168565},
        // Reference values from mpmath
        test_case{-7.5, 9.744766282170433}, test_case{-3.2, 28.29818640219408}
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

TEST_CASE("trigamma gamma-zeta consistency", "[trigamma][xsf_tests]") {
    const double rtol = 1e-14;
    auto x = GENERATE(-7.5, -3.2, -0.5, 1e-6, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0, 1e6);
    double result_trigamma = xsf::trigamma(x);
    double ref_trigamma = xsf::gamma(2.0) * xsf::zeta(2.0, x);
    double rel_err_trigamma = xsf::extended_relative_error(result_trigamma, ref_trigamma);

    CAPTURE(x, result_trigamma, ref_trigamma, rel_err_trigamma, rtol);
    REQUIRE(rel_err_trigamma <= rtol);
}

TEST_CASE("trigamma complex", "[trigamma][xsf_tests]") {
    double rtol = 1e-14;
    using test_case = std::tuple<std::complex<double>, std::complex<double>>;

    auto [z, ref_trigamma] = GENERATE(
        // Reference values from
        // https://github.com/JuliaMath/SpecialFunctions.jl/blob/c9ff36085c9008b9bdaa457e745f459a1827e52e/test/gamma.jl#L261-L294
        test_case{{8.0, 0.0}, {0.133137014694031425134546685920401606452509991909746283540546, 0.0}},
        test_case{
            {0.0, 8.0},
            {-0.0078125000000000000029194973110119898029284994355721719150,
             -0.12467345030312762782439017882063360876391046513966063947}
        },
        test_case{
            {-3.2, 0.1},
            {15.2073506449733631753218003030676132587307964766963426965699,
             15.7081038855113567966903832015076316497656334265029416039199}
        }
    );

    std::complex<double> result_trigamma = xsf::trigamma(z);
    double rel_err_trigamma = xsf::extended_relative_error(result_trigamma, ref_trigamma);

    CAPTURE(z, result_trigamma, ref_trigamma, rel_err_trigamma, rtol);
    REQUIRE(rel_err_trigamma <= rtol);
}

TEST_CASE("trigamma float overloads", "[trigamma][xsf_tests]") {
    const double rtol = 1e-6;

    SECTION("real") {
        using test_case = std::tuple<float, double>;
        // Reference values from the double-precision tests above.
        auto [z, ref_trigamma] = GENERATE(
            test_case{0.5f, M_PI * M_PI / 2}, test_case{10.0f, 0.10516633568168565}, test_case{-7.5f, 9.744766282170433}
        );
        float result_trigamma = xsf::trigamma(z);
        double rel_err_trigamma = xsf::extended_relative_error(static_cast<double>(result_trigamma), ref_trigamma);
        CAPTURE(z, result_trigamma, ref_trigamma, rel_err_trigamma, rtol);
        REQUIRE(rel_err_trigamma <= rtol);
    }

    SECTION("complex") {
        using test_case = std::tuple<std::complex<float>, std::complex<double>>;
        // Reference values from the double-precision tests above.
        auto [z, ref_trigamma] = GENERATE(
            test_case{{0.0f, 8.0f}, {-0.0078125, -0.12467345030312763}},
            test_case{{-3.2f, 0.1f}, {15.207350644973363, 15.708103885511357}}
        );
        std::complex<float> result_trigamma = xsf::trigamma(z);
        double rel_err_trigamma =
            xsf::extended_relative_error(static_cast<std::complex<double>>(result_trigamma), ref_trigamma);
        CAPTURE(z, result_trigamma, ref_trigamma, rel_err_trigamma, rtol);
        REQUIRE(rel_err_trigamma <= rtol);
    }
}
