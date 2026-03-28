#include "../testing_utils.h"
#include <tuple>
#include <xsf/lambertw.h>

TEST_CASE("lambertw k=0 branch", "[lambertw][xsf_tests]") {

    using test_case = std::tuple<double, double>;
    // Reference values were computed with the Python library mpmath.
    // import mpmath
    // x0 = [-1/mpmath.e+1e-3, -0.3, -0.1, -1e-6, 0.0, 1e-6, 1, 10, 1000, 1000000]
    // expected_k0 = [float(mpmath.lambertw(x, k=0)) for x in x0]

    auto [x, ref_w] = GENERATE(
        test_case{-0.366879441171442, -0.9280201500545675},
        test_case{-0.3, -0.4894022271802149},
        test_case{-0.1, -0.11183255915896297},
        test_case{-1e-06, -1.0000010000014999e-06},
        test_case{0.0, 0.0},
        test_case{1e-06, 9.999990000015e-07},
        test_case{1, 0.5671432904097838},
        test_case{10, 1.7455280027406994},
        test_case{1000, 5.249602852401596},
        test_case{1000000, 11.383358086140053}
    );

    double w_d = xsf::lambertw(x, 0, 1e-14);
    const auto abs_error_d = xsf::extended_absolute_error(w_d, ref_w);
    CAPTURE(x, w_d, ref_w, abs_error_d);
    REQUIRE(abs_error_d <= 1e-11);

    float w_f = xsf::lambertw(static_cast<float>(x), 0, 1e-14);
    const auto abs_error_f = xsf::extended_absolute_error(static_cast<double>(w_f), ref_w);
    CAPTURE(x, w_f, ref_w, abs_error_f);
    REQUIRE(abs_error_f <= 1e-6);
}

TEST_CASE("lambertw k=-1 branch", "[lambertw][xsf_tests]") {

    using test_case = std::tuple<double, double>;
    // Reference values were computed with the Python library mpmath.
    // import mpmath
    // xn1 = [-1/mpmath.e+1e-3, -0.3, -0.1, -1e-3]
    // expected_kn1 = [float(mpmath.lambertw(x, k=-1)) for x in xn1]

    auto [x, ref_w] = GENERATE(
        test_case{-0.366879441171442, -1.0756089411866245},
        test_case{-0.3, -1.7813370234216277},
        test_case{-0.1, -3.577152063957297},
        test_case{-0.001, -9.11800647040274}
    );

    double w_d = xsf::lambertw(x, -1, 1e-14);
    const auto abs_error_d = xsf::extended_absolute_error(w_d, ref_w);
    CAPTURE(x, w_d, ref_w, abs_error_d);
    REQUIRE(abs_error_d <= 1e-11);

    float w_f = xsf::lambertw(static_cast<float>(x), -1, 1e-14);
    const auto abs_error_f = xsf::extended_absolute_error(static_cast<double>(w_f), ref_w);
    CAPTURE(x, w_f, ref_w, abs_error_f);
    REQUIRE(abs_error_f <= 1e-6);
}

TEST_CASE("lambertw nan inf", "[lambertw][xsf_tests]") {

    double nan = std::numeric_limits<double>::quiet_NaN();
    double inf = std::numeric_limits<double>::infinity();

    double w = xsf::lambertw(nan, 0, 1e-14);
    REQUIRE(std::isnan(w));

    w = xsf::lambertw(inf, 0, 1e-14);
    REQUIRE(w == std::numeric_limits<double>::infinity());

    w = xsf::lambertw(-inf, 0, 1e-14);
    REQUIRE(w == std::numeric_limits<double>::infinity());
}