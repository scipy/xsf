#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <tuple>
#include <xsf/gamma.h>

TEST_CASE("gammainc large-shape lower tail (scipy gh-26535)", "[gammainc][xsf_tests]") {
    using case_type = std::tuple<double, double, double>;
    // References computed with mpmath.gammainc(a, 0, x, regularized=True) at 100 digits.
    auto [a, x, expected] = GENERATE(
        case_type{0x1.1bb9d5a718bb3p+12, 0x1.49d468f0c82d9p+11, 1.0535831861145413e-246},
        case_type{0x1.74d80466f9134p+12, 0x1.b86ace35d5d52p+11, 2.8344710632721627e-306},
        case_type{0x1.1b4e937e46dcbp+12, 0x1.5305961513509p+11, 6.239827729196435e-223},
        case_type{0x1.1f805570f193bp+12, 0x1.4f838a55aa9c9p+11, 1.001757841288208e-246}
    );
    double actual = xsf::gammainc(a, x);
    CAPTURE(a, x, actual, expected);
    REQUIRE(std::abs((actual - expected) / expected) < 1e-12);
}

TEST_CASE("incomplete gamma prefactor switching boundaries", "[gammainc][gammaincc][xsf_tests]") {
    using case_type = std::tuple<double, double, double, double>;
    // References computed with mpmath.gammainc at 100 digits for both tails.
    auto [a, x, expected_p, expected_q] = GENERATE(
        case_type{199.999999, 100, 9.3431566174655e-19, 1.0}, case_type{200, 100, 9.34315007298839e-19, 1.0},
        case_type{200.000001, 100, 9.343143528515819e-19, 1.0}, case_type{300, 90, 7.152710770691828e-68, 1.0},
        case_type{1000, 490, 1.1988275832107112e-90, 1.0}, case_type{1000, 500, 3.2982727970671e-86, 1.0},
        case_type{1000, 550, 1.7446109896636843e-66, 1.0}, case_type{1000, 599.999999, 2.324058502834769e-50, 1.0},
        case_type{1000, 600, 2.324060057960917e-50, 1.0}, case_type{1000, 600.000001, 2.324061613088099e-50, 1.0},
        case_type{200, 280, 0.9999997945244374, 2.0547556257375757e-07},
        case_type{1000, 1400, 1.0, 8.04132540818483e-30}, case_type{1000, 1400.000001, 1.0, 8.041323090850252e-30},
        case_type{1000, 1500, 1.0, 2.204698611388996e-43}, case_type{1000, 1800, 1.0, 1.0801622506932697e-94},
        case_type{4539.614661308867, 6809.4219919633, 1.0, 4.948947801153583e-189}
    );
    double actual_p = xsf::gammainc(a, x);
    double actual_q = xsf::gammaincc(a, x);
    CAPTURE(a, x, actual_p, expected_p, actual_q, expected_q);
    REQUIRE(std::abs((actual_p - expected_p) / expected_p) < 1e-12);
    REQUIRE(std::abs((actual_q - expected_q) / expected_q) < 1e-12);
}

TEST_CASE("incomplete gamma large-shape underflow", "[gammainc][gammaincc][xsf_tests]") {
    REQUIRE(xsf::gammainc(200.0, 1e-6) == 0.0);
    REQUIRE(xsf::gammaincc(200.0, 1e-6) == 1.0);
    REQUIRE(xsf::gammainc(200.0, 20000.0) == 1.0);
    REQUIRE(xsf::gammaincc(200.0, 20000.0) == 0.0);
    REQUIRE(xsf::gammainc(1e300, 1e-100) == 0.0);
    REQUIRE(xsf::gammaincc(1e300, 1e-100) == 1.0);
    REQUIRE(xsf::gammainc(200.0, 1e308) == 1.0);
    REQUIRE(xsf::gammaincc(200.0, 1e308) == 0.0);
}
