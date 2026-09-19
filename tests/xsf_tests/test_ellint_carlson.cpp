#include "../testing_utils.h"

#include <xsf/cpu/ellint_carlson.h>
#include <xsf/ellip.h>

namespace {

using C = std::complex<double>;

} // namespace

TEST_CASE("Carlson RC", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L1868-L1886
    SECTION("real values") {
        using test_case = std::tuple<double, double, double>;
        auto [x, y, expected] = GENERATE(
            test_case{1.0, 1.0, 1.0}, test_case{0.0, 0.25, M_PI}, test_case{2.25, 2.0, std::log(2.0)},
            test_case{0.25, -2.0, std::log(2.0) / 3.0}
        );
        const double result = xsf::cpu::elliprc(x, y);
        const double error = xsf::extended_relative_error(result, expected);
        CAPTURE(x, y, result, expected, error);
        REQUIRE(error <= 1e-13);
    }

    REQUIRE(xsf::cpu::elliprc(1.0, std::numeric_limits<double>::infinity()) == 0.0);
    REQUIRE(std::isnan(xsf::cpu::elliprc(1.0, 0.0)));
    REQUIRE(xsf::cpu::elliprc(C{1.0}, C{1.0, std::numeric_limits<double>::infinity()}) == C{0.0});

    using test_case = std::tuple<C, C, C>;
    auto [x, y, expected] = GENERATE(
        test_case{0.0, 0.25, M_PI}, test_case{2.25, 2.0, std::log(2.0)},
        test_case{0.0, C{0.0, 1.0}, C{1.1107207345396, -1.1107207345396}},
        test_case{C{0.0, -1.0}, C{0.0, 1.0}, C{1.2260849569072, -0.34471136988768}},
        test_case{0.25, -2.0, std::log(2.0) / 3.0}, test_case{C{0.0, 1.0}, -1.0, C{0.77778596920447, 0.19832484993429}}
    );
    const C result = xsf::cpu::elliprc(x, y);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(x, y, result, expected, error);
    REQUIRE(error <= 1e-13);
}

TEST_CASE("Carlson RD", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L1888-L1910
    SECTION("real values") {
        using test_case = std::tuple<double, double, double, double>;
        auto [x, y, z, expected] = GENERATE(
            test_case{1.0, 1.0, 1.0, 1.0}, test_case{0.0, 2.0, 1.0, 3.0 * 0.59907011736779610371},
            test_case{2.0, 3.0, 4.0, 0.16510527294261}
        );
        const double result = xsf::cpu::elliprd(x, y, z);
        const double error = xsf::extended_relative_error(result, expected);
        CAPTURE(x, y, z, result, expected, error);
        REQUIRE(error <= 1e-13);
    }

    REQUIRE(xsf::cpu::elliprd(1.0, 1.0, std::numeric_limits<double>::infinity()) == 0.0);
    REQUIRE(std::isinf(xsf::cpu::elliprd(1.0, 1.0, 0.0)));
    REQUIRE(std::isnan(xsf::cpu::elliprd(1.0, 1.0, -1.0)));
    REQUIRE(std::isinf(xsf::cpu::elliprd(C{1.0}, C{1.0}, C{0.0}).real()));
    REQUIRE(std::isinf(xsf::cpu::elliprd(C{0.0}, C{1.0}, C{0.0}).real()));
    REQUIRE(std::isnan(xsf::cpu::elliprd(1.0, 1.0, -std::numeric_limits<double>::min() / 2.0)));
    REQUIRE(std::isnan(xsf::cpu::elliprd(C{1.0}, C{1.0}, C{-1.0}).real()));

    using test_case = std::tuple<C, C, C, C>;
    auto [x, y, z, expected] = GENERATE(
        test_case{0.0, 2.0, 1.0, 1.7972103521034}, test_case{2.0, 3.0, 4.0, 0.16510527294261},
        test_case{C{0.0, 1.0}, C{0.0, -1.0}, 2.0, 0.65933854154220},
        test_case{0.0, C{0.0, 1.0}, C{0.0, -1.0}, C{1.2708196271910, 2.7811120159521}},
        test_case{0.0, C{-1.0, 1.0}, C{0.0, 1.0}, C{-1.8577235439239, -0.96193450888839}},
        test_case{C{-2.0, -1.0}, C{0.0, -1.0}, C{-1.0, 1.0}, C{1.8249027393704, -1.2218475784827}}
    );
    const C result = xsf::cpu::elliprd(x, y, z);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(x, y, z, result, expected, error);
    REQUIRE(error <= 1e-13);
}

TEST_CASE("Carlson RF", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L1912-L1935
    SECTION("real values") {
        using test_case = std::tuple<double, double, double, double>;
        auto [x, y, z, expected] =
            GENERATE(test_case{1.0, 1.0, 1.0, 1.0}, test_case{0.0, 1.0, 2.0, 1.31102877714605990523});
        const double result = xsf::cpu::elliprf(x, y, z);
        const double error = xsf::extended_relative_error(result, expected);
        CAPTURE(x, y, z, result, expected, error);
        REQUIRE(error <= 1e-13);
    }

    REQUIRE(xsf::cpu::elliprf(1.0, std::numeric_limits<double>::infinity(), 1.0) == 0.0);
    REQUIRE(std::isinf(xsf::cpu::elliprf(0.0, 1.0, 0.0)));
    REQUIRE(std::isnan(xsf::cpu::elliprf(1.0, 1.0, -1.0)));
    REQUIRE(xsf::cpu::elliprf(C{std::numeric_limits<double>::infinity()}, C{0.0}, C{1.0}) == C{0.0});
    REQUIRE(std::isnan(xsf::cpu::elliprf(C{1.0}, C{1.0}, C{-std::numeric_limits<double>::infinity(), 1.0}).real()));

    using test_case = std::tuple<C, C, C, C>;
    auto [x, y, z, expected] = GENERATE(
        test_case{1.0, 2.0, 0.0, 1.3110287771461}, test_case{C{0.0, 1.0}, C{0.0, -1.0}, 0.0, 1.8540746773014},
        test_case{0.5, 1.0, 0.0, 1.8540746773014},
        test_case{C{-1.0, 1.0}, C{0.0, 1.0}, 0.0, C{0.79612586584234, -1.2138566698365}},
        test_case{2.0, 3.0, 4.0, 0.58408284167715}, test_case{C{0.0, 1.0}, C{0.0, -1.0}, 2.0, 1.0441445654064},
        test_case{C{-1.0, 1.0}, C{0.0, 1.0}, C{1.0, -1.0}, C{0.93912050218619, -0.53296252018635}}
    );
    const C result = xsf::cpu::elliprf(x, y, z);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(x, y, z, result, expected, error);
    REQUIRE(error <= 1e-13);
}

TEST_CASE("Carlson RG", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L1937-L1956
    SECTION("real values") {
        using test_case = std::tuple<double, double, double, double>;
        auto [x, y, z, expected] =
            GENERATE(test_case{1.0, 1.0, 1.0, 1.0}, test_case{0.0, 0.0, 1.0, 0.5}, test_case{0.0, 16.0, 16.0, M_PI});
        const double result = xsf::cpu::elliprg(x, y, z);
        const double error = xsf::extended_relative_error(result, expected);
        CAPTURE(x, y, z, result, expected, error);
        REQUIRE(error <= 1e-13);
    }

    REQUIRE(xsf::cpu::elliprg(0.0, 0.0, 0.0) == 0.0);
    REQUIRE(std::isinf(xsf::cpu::elliprg(1.0, std::numeric_limits<double>::infinity(), 1.0)));
    REQUIRE(std::isinf(xsf::cpu::elliprg(C{std::numeric_limits<double>::infinity()}, C{1.0}, C{1.0}).real()));

    using test_case = std::tuple<C, C, C, C>;
    auto [x, y, z, expected] = GENERATE(
        test_case{0.0, 16.0, 16.0, M_PI}, test_case{2.0, 3.0, 4.0, 1.7255030280692},
        test_case{0.0, C{0.0, 1.0}, C{0.0, -1.0}, 0.42360654239699},
        test_case{C{-1.0, 1.0}, C{0.0, 1.0}, 0.0, C{0.44660591677018, 0.70768352357515}},
        test_case{C{0.0, -1.0}, C{-1.0, 1.0}, C{0.0, 1.0}, C{0.36023392184473, 0.40348623401722}},
        test_case{0.0, 0.0796, 4.0, 1.0284758090288}
    );
    const C result = xsf::cpu::elliprg(x, y, z);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(x, y, z, result, expected, error);
    REQUIRE(error <= 1e-13);
}

TEST_CASE("Carlson RJ", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L1958-L1983
    SECTION("real values") {
        using test_case = std::tuple<double, double, double, double, double>;
        auto [x, y, z, p, expected] = GENERATE(
            test_case{1.0, 1.0, 1.0, 1.0, 1.0}, test_case{0.0, 1.0, 2.0, 3.0, 0.77688623778582},
            test_case{2.0, 3.0, 4.0, -0.5, 0.24723819703052}, test_case{2.0, 3.0, 4.0, -5.0, -0.12711230042964}
        );
        const double result = xsf::cpu::elliprj(x, y, z, p);
        const double error = xsf::extended_relative_error(result, expected);
        CAPTURE(x, y, z, p, result, expected, error);
        REQUIRE(error <= 1e-13);
    }

    REQUIRE(xsf::cpu::elliprj(1.0, 1.0, 1.0, std::numeric_limits<double>::infinity()) == 0.0);
    REQUIRE(std::isnan(xsf::cpu::elliprj(-1.0, 1.0, 1.0, 1.0)));
    REQUIRE(xsf::cpu::elliprj(1.0, 1.0, std::numeric_limits<double>::infinity(), 1.0) == 0.0);
    REQUIRE(std::isnan(xsf::cpu::elliprj(1.0, 0.0, 0.0, 0.0)));

    using test_case = std::tuple<C, C, C, C, C>;
    auto [x, y, z, p, expected] = GENERATE(
        test_case{0.0, 1.0, 2.0, 3.0, 0.77688623778582}, test_case{2.0, 3.0, 4.0, 5.0, 0.14297579667157},
        test_case{2.0, 3.0, 4.0, C{-1.0, 1.0}, C{0.13613945827771, -0.38207561624427}},
        test_case{C{0.0, 1.0}, C{0.0, -1.0}, 0.0, 2.0, 1.6490011662711},
        test_case{C{-1.0, 1.0}, C{-1.0, -1.0}, 1.0, 2.0, 0.94148358841220},
        test_case{C{0.0, 1.0}, C{0.0, -1.0}, 0.0, C{1.0, -1.0}, C{1.8260115229009, 1.2290661908643}},
        test_case{C{-1.0, 1.0}, C{-1.0, -1.0}, 1.0, C{-3.0, 1.0}, C{-0.61127970812028, -1.0684038390007}},
        // Cauchy principal values.
        test_case{2.0, 3.0, 4.0, -0.5, 0.24723819703052}, test_case{2.0, 3.0, 4.0, -5.0, -0.12711230042964}
    );
    const C result = xsf::cpu::elliprj(x, y, z, p);
    const double error = xsf::extended_relative_error(result, expected);
    CAPTURE(x, y, z, p, result, expected, error);
    REQUIRE(error <= 1e-13);
}

TEST_CASE("Carlson RJ difficult arguments", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L1985-L1998
    // SciPy marks this test xfail for insufficient accuracy on 32-bit.
    using test_case = std::tuple<double, double, double, double, double>;
    auto [x, y, z, p, expected] = GENERATE(
        test_case{
            6.483625725195452e-08, 1.1649136528196886e-27, 3.6767340167168e13, 0.493704617023468,
            8.63426920644241857617477551054e-6
        },
        test_case{
            14.375105857849121, 9.993988969725365e-11, 1.72844262269944e-26, 5.898871222598245e-06,
            829774.1424801627252574054378691828
        }
    );
    const double result = xsf::cpu::elliprj(x, y, z, p);
    CAPTURE(x, y, z, p, result, expected);
    REQUIRE(std::abs(result - expected) <= 1e-20 + 5e-15 * std::abs(expected));
}

TEST_CASE("Legendre-Carlson K and E identities", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L2001-L2024
    // and
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L2037-L2042
    auto check_identities = [](double m) {
        const double k = xsf::ellipk(m);
        const double rf = xsf::cpu::elliprf(0.0, 1.0 - m, 1.0);
        const double e = xsf::ellipe(m);
        const double rg = 2.0 * xsf::cpu::elliprg(0.0, 1.0 - m, 1.0);
        CAPTURE(m, k, rf, e, rg);
        REQUIRE(xsf::extended_relative_error(k, rf) <= 1e-7);
        REQUIRE(xsf::extended_relative_error(e, rg) <= 1e-7);
    };
    check_identities(std::numeric_limits<double>::lowest());
    for (int exponent = 1023; exponent > 0; --exponent) {
        check_identities(-std::ldexp(1.0, exponent));
    }
    for (int i = 0; i < 200; ++i) {
        check_identities(-1.0 + 0.01 * i);
    }
}

TEST_CASE("Legendre-Carlson Km1 identity", "[ellint_carlson][xsf_tests]") {
    // Mirrors
    // https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/tests/test_basic.py#L2026-L2035
    for (int exponent = -1022; exponent < 0; ++exponent) {
        const double m1 = std::ldexp(1.0, exponent);
        const double result = xsf::ellipkm1(m1);
        const double expected = xsf::cpu::elliprf(0.0, m1, 1.0);
        const double error = xsf::extended_relative_error(result, expected);
        CAPTURE(m1, result, expected, error);
        REQUIRE(error <= 1e-7);
    }
}
