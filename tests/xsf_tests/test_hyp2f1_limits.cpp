#include "../testing_utils.h"

#include <xsf/hyp2f1.h>

TEST_CASE("hyp2f1 bounds fixed parameter sums", "[hyp2f1][xsf_tests]") {
    using C = std::complex<double>;
    using test_case = std::tuple<double, double, double, C>;
    auto [a, b, c, z] = GENERATE(
        test_case{-1e15, 1.0, 2.0, C{0.5, 0.0}}, test_case{1.0, -1e15, 2.0, C{0.5, 0.0}},
        test_case{1e15 + 0.5, 1.25, 0.5, C{0.5, 0.0}}, test_case{1.25, 1.5, 1e17, C{0.8, 0.7}},
        test_case{1e15 + 0.25, 1.25, 0.5, C{0.8, 0.5}}, test_case{1e15 + 0.25, 0.25, 0.75, C{2.0, 0.5}},
        test_case{-0x1p64, 1.0, 2.0, C{0.5, 0.0}}
    );
    const auto actual = xsf::hyp2f1(a, b, c, z);
    CAPTURE(a, b, c, z, actual);
    CHECK(std::isnan(actual.real()));
    CHECK(std::isnan(actual.imag()));
}

TEST_CASE("hyp2f1 large c in the unit disk", "[hyp2f1][xsf_tests]") {
    using C = std::complex<double>;
    using test_case = std::tuple<double, double, C>;
    // Gauss-series reference values computed at 120 digits.
    auto [b, c, expected] = GENERATE(
        test_case{2.14, 1e6, C{1.0000060125429742, 6.106556144237689e-6}},
        test_case{2.14, 1e12, C{1.0000000000060125, 6.106490000066144e-12}}, test_case{2.14, 1e17, C{1.0, 6.10649e-17}},
        test_case{-2.14, 1e6, C{0.9999939874556276, -6.1064659859879145e-6}},
        test_case{-2.14, 1e12, C{0.9999999999939875, -6.106489999975986e-12}},
        test_case{-2.14, 1e17, C{0.9999999999999999, -6.106489999999999e-17}}
    );
    const auto actual = xsf::hyp2f1(4.39, b, c, C{0.64, 0.65});
    CAPTURE(b, c, actual, expected);
    CHECK(xsf::extended_relative_error(actual, expected) <= 1e-14);
    CHECK(xsf::extended_relative_error(actual.imag(), expected.imag()) <= 1e-14);
}

TEST_CASE("hyp2f1 keeps the shorter terminating polynomial", "[hyp2f1][xsf_tests]") {
    using C = std::complex<double>;
    using test_case = std::tuple<double, double, double, C, C>;
    auto [a, b, c, z, expected] = GENERATE(
        test_case{-1e15, -2.0, 1.5, C{0.5, 0.0}, C{6.6666666666667265e28, 0.0}},
        test_case{1.5, 1e15 + 0.5, 0.5, C{0.0, 1e-15}, C{-1.142639663747656, 1.9220755965441734}}
    );
    const auto actual = xsf::hyp2f1(a, b, c, z);
    CAPTURE(a, b, c, z, actual, expected);
    CHECK(xsf::extended_relative_error(actual, expected) <= 1e-14);
}
