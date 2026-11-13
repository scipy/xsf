#include "../testing_utils.h"

#include <xsf/cpu/par_cyl.h>

TEST_CASE("pbdv and pbvv derivative when |v| + 1 rounds up", "[pbdv][pbvv][xsf_tests]") {
    // Regression test for gh-401: the work buffer was one element too short.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import math

    import mpmath

    mpmath.mp.dps = 50


    def pbdv_d(v, x):
        return mpmath.diff(lambda t: mpmath.pcfd(v, t), x)


    def pbvv_d(v, x):
        return mpmath.diff(lambda t: mpmath.pcfv(-mpmath.mpf(v) - 0.5, t), x)


    for fn, v, x in [
        (pbdv_d, math.nextafter(2, 0), 2.0),
        (pbdv_d, -math.nextafter(2, 0), 2.0),
        (pbvv_d, math.nextafter(1, 0), 2.0),
        (pbvv_d, -math.nextafter(1, 0), 2.0),
    ]:
        print(f"{fn.__name__}({v!r}, {x!r}) = {float(fn(v, x))!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<bool, double, double, double>;
    auto [is_pbdv, v, x, desired] = GENERATE(
        test_case{true, 1.9999999999999998, 2.0, 0.367879441171442},
        test_case{true, -1.9999999999999998, 2.0, -0.09715978862055018},
        test_case{false, 0.9999999999999999, 2.0, 0.780820901800019},
        test_case{false, -0.9999999999999999, 2.0, 2.168875102838455}
    );

    double f, d;
    if (is_pbdv) {
        xsf::pbdv(v, x, f, d);
    } else {
        xsf::pbvv(v, x, f, d);
    }
    const double error = xsf::extended_relative_error(d, desired);
    const double rtol = 1e-14;
    CAPTURE(is_pbdv, v, x, d, desired, error, rtol);
    REQUIRE(error <= rtol);
}

TEST_CASE("pbdv and pbvv at v = -0.0", "[pbdv][pbvv][xsf_tests]") {
    // Regression test for gh-418: v = -0.0 took the negative-order path.
    SET_FP_FORMAT()
    double f_pos, d_pos, f_neg, d_neg;

    xsf::pbdv(0.0, 50.0, f_pos, d_pos);
    xsf::pbdv(-0.0, 50.0, f_neg, d_neg);
    CHECK(f_neg == f_pos);
    CHECK(d_neg == d_pos);

    xsf::pbvv(0.0, 5.0, f_pos, d_pos);
    xsf::pbvv(-0.0, 5.0, f_neg, d_neg);
    CHECK(f_neg == f_pos);
    CHECK(d_neg == d_pos);
}

TEST_CASE("pbwa accuracy", "[pbwa][xsf_tests]") {
    // Regression test for gh-450: cancellation where W(a, x) is small, and inaccurate a = 0 case.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50

    for a, x in [(5.0, 4.5), (0.0, -5.0), (0.0, 0.25)]:
        w = mpmath.pcfw(a, x)
        wd = mpmath.diff(lambda t: mpmath.pcfw(a, t), x)
        print(f"{a!r}, {x!r}, {float(w)!r}, {float(wd)!r}")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [a, x, desired_f, desired_d] = GENERATE(
        test_case{5.0, 4.5, 0.00025418079605278414, -0.0002554687616971141},
        test_case{0.0, -5.0, 0.6599633712356283, -1.7552224081475016},
        test_case{0.0, 0.25, 0.9004707748272274, -0.49008286316884153}
    );

    double f, d;
    xsf::pbwa(a, x, f, d);
    const double error_f = xsf::extended_relative_error(f, desired_f);
    const double error_d = xsf::extended_relative_error(d, desired_d);
    const double rtol = 2e-14;
    CAPTURE(a, x, f, desired_f, error_f, d, desired_d, error_d, rtol);
    CHECK(error_f <= rtol);
    CHECK(error_d <= rtol);
}

TEST_CASE("pbwa accuracy float", "[pbwa][xsf_tests]") {
    // The same values as above (gh-450), for the float instantiation.
    SET_FP_FORMAT()
    using test_case = std::tuple<float, float, double, double>;
    auto [a, x, desired_f, desired_d] = GENERATE(
        test_case{5.0f, 4.5f, 0.00025418079605278414, -0.0002554687616971141},
        test_case{0.0f, -5.0f, 0.6599633712356283, -1.7552224081475016},
        test_case{0.0f, 0.25f, 0.9004707748272274, -0.49008286316884153}
    );

    float f, d;
    xsf::pbwa(a, x, f, d);
    const double error_f = xsf::extended_relative_error(static_cast<double>(f), desired_f);
    const double error_d = xsf::extended_relative_error(static_cast<double>(d), desired_d);
    const double rtol = 1e-5;
    CAPTURE(a, x, f, desired_f, error_f, d, desired_d, error_d, rtol);
    CHECK(error_f <= rtol);
    CHECK(error_d <= rtol);
}

TEST_CASE("pbvv half-integer v and negative x", "[pbvv][xsf_tests]") {
    // Regression test for gh-520: V_v(x) is recessive here and lost all digits.

    /*
    If you want to check the expected values below with mpmath, this is the script:

    import mpmath

    mpmath.mp.dps = 50


    def pbvv(v, x):
        v, x = mpmath.mpf(v), mpmath.mpf(x)
        f = lambda t: mpmath.pcfv(-v - 0.5, t)
        return f(x), mpmath.diff(f, x)


    for v, x in [(-0.5, -8.0), (-0.5, -12.0), (0.5, -12.0), (0.5, -8.0), (-0.5, -4.0)]:
        f, d = pbvv(v, x)
        print(f"pbvv({v!r}, {x!r}) = ({float(f)!r}, {float(d)!r})")
    */
    SET_FP_FORMAT()
    using test_case = std::tuple<double, double, double, double>;
    auto [v, x, desired_f, desired_d] = GENERATE(
        test_case{-0.5, -8.0, 2.232022159445888e-08, 9.064507251631539e-08},
        test_case{-0.5, -12.0, 3.76805638742432e-17, 2.2763748814178562e-16},
        test_case{0.5, -12.0, -9.074417427744898e-16, -5.406969892772695e-15}
    );

    double f, d;
    xsf::pbvv(v, x, f, d);
    const double error_f = xsf::extended_relative_error(f, desired_f);
    const double error_d = xsf::extended_relative_error(d, desired_d);
    const double rtol = 1e-8;
    CAPTURE(v, x, f, d, desired_f, desired_d, error_f, error_d, rtol);
    REQUIRE(error_f <= rtol);
    REQUIRE(error_d <= rtol);
}

TEST_CASE("pbvv half-integer v and negative x, float", "[pbvv][xsf_tests]") {
    // Regression test for gh-520, float instantiation. Reference values from the script above:
    // pbvv(0.5, -8.0) = (-3.5985191778830186e-07, -1.4170874495587485e-06)
    // pbvv(-0.5, -4.0) = (0.005059046431512129, 0.01070212790594681)
    SET_FP_FORMAT()
    using test_case = std::tuple<float, float, double, double>;
    auto [v, x, desired_f, desired_d] = GENERATE(
        test_case{0.5f, -8.0f, -3.5985191778830186e-07, -1.4170874495587485e-06},
        test_case{-0.5f, -4.0f, 0.005059046431512129, 0.01070212790594681}
    );
    float f, d;
    xsf::pbvv(v, x, f, d);
    const double rtol = 1e-5;
    CAPTURE(v, x, f, d, rtol);
    REQUIRE(xsf::extended_relative_error(static_cast<double>(f), desired_f) <= rtol);
    REQUIRE(xsf::extended_relative_error(static_cast<double>(d), desired_d) <= rtol);
}
