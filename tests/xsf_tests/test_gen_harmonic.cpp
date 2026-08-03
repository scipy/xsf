#include "../testing_utils.h"

#include <xsf/gen_harmonic.h>

static constexpr double rtol = 5e-15;
static const double INF = std::numeric_limits<double>::infinity();
static const double NaN = std::numeric_limits<double>::quiet_NaN();

TEST_CASE("gen harmonic", "[gen_harmonic][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_gen_harmonic.py#L11-L39
    using test_case = std::tuple<int, double, double>;
    auto [n, a, expected] = GENERATE(
        test_case{8, 9.0, 1.0020083884212339}, test_case{1000, 2.5, 1.3414661912046497},
        test_case{10, 1.5, 1.9953364933456017}, test_case{10000, 1.25, 4.1951168257387765},
        test_case{10000, 1.00001, 9.787182620770265}, test_case{80, 1.000002, 4.965460167788836},
        test_case{75, 1 + 1e-12, 4.901355630543771}, test_case{100, 1 + 1e-14, 5.187377517639515},
        test_case{100, 1 + 8e-16, 5.187377517639611}, test_case{100, 1.0, 5.187377517639621},
        test_case{7, 1.0, 2.592857142857143}, test_case{8000, 1.0, 9.564474984261423},
        test_case{5, 1 - 1e-12, 2.2833333333347143}, test_case{25000, 1 - 1e-12, 10.703866768669737},
        test_case{1000, 0.995, 7.6058022857089975}, test_case{1000, 0.75, 19.055178975831392},
        test_case{10000, 0.25, 1332.5700547197382}, test_case{5, 1e-8, 4.999999952125083},
        test_case{15, 1e-16, 14.999999999999996}, test_case{100, 0.0, 100.0}, test_case{4, -1.0, 10.0},
        test_case{75, -1.5, 19811.38815892374}
    );

    double output = xsf::gen_harmonic(n, a);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(n, a, output, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("gen harmonic exact cases", "[gen_harmonic][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_gen_harmonic.py#L42-L53
    using test_case = std::tuple<int, double, double>;
    auto [n, a, expected] = GENERATE(
        test_case{10, INF, 1.0}, test_case{1, NaN, 1.0}, test_case{1, -INF, 1.0}, test_case{3, NaN, NaN},
        test_case{-3, 1.0, NaN}
    );

    double output = xsf::gen_harmonic(n, a);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(n, a, output, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("gen harmonic nan inputs", "[gen_harmonic][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_gen_harmonic.py#L56-L58
    using test_case = std::tuple<double, double, double>;
    auto [n, a, expected] = GENERATE(test_case{NaN, 0.75, NaN});

    double output = xsf::gen_harmonic(n, a);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(n, a, output, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("normalized gen harmonic", "[gen_harmonic][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_gen_harmonic.py#L61-L76
    using test_case = std::tuple<int, int, int, double, double>;
    auto [j, k, n, a, expected] = GENERATE(
        test_case{400, 5000, 5000, 10.0, 4.2821759663214485e-25}, test_case{400, 5000, 5000, 3.5, 1.11086549102426e-07},
        test_case{1, 2, 3, 1.5, 0.8755176866163012}, test_case{300, 500, 500, 1 + 1e-14, 0.07559343891632035},
        test_case{1500, 2500, 3000, 1 - 1e-12, 0.05957291246371843}, test_case{10, 12, 16, 0.5, 0.13601665344521513},
        test_case{16, 16, 20, 0.125, 0.04583107002260924}, test_case{10, 12, 16, -0.5, 0.22359306724308234},
        test_case{1, 8000, 10000, -1.5, 0.5724512895513029}
    );

    double output = xsf::normalized_gen_harmonic(j, k, n, a);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(j, k, n, a, output, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("normalized gen harmonic exact cases", "[gen_harmonic][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_gen_harmonic.py#L79-L94
    using test_case = std::tuple<int, int, int, double, double>;
    auto [j, k, n, a, expected] = GENERATE(
        test_case{1, 1, 1, 0.5, 1.0}, test_case{1, 1, 1, NaN, 1.0}, test_case{1, 2, 5, NaN, NaN},
        test_case{1, 2, 1, 1.25, NaN}, test_case{1, 2, 3, INF, 1.0}, test_case{2, 3, 4, INF, 0.0},
        test_case{1, 1, 10, -INF, 0.0}, test_case{2, 3, 4, -INF, NaN}, test_case{3, 6, 8, 0.0, 0.5}
    );

    double output = xsf::normalized_gen_harmonic(j, k, n, a);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(j, k, n, a, output, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}

TEST_CASE("normalied gen harmonic nan inputs", "[gen_harmonic][xsf_tests]") {
    // https://github.com/scipy/scipy/blob/a125578782dd3213fb57fda0f4b97c70dd054b1d/scipy/special/tests/test_gen_harmonic.py#L97-L99
    using test_case = std::tuple<double, double, double, double, double>;
    auto [j, k, n, a, expected] = GENERATE(test_case{1.0, NaN, 10.0, 1.05, NaN});

    double output = xsf::normalized_gen_harmonic(j, k, n, a);
    double rel_error = xsf::extended_relative_error(output, expected);
    CAPTURE(j, k, n, a, output, expected, rel_error);
    REQUIRE(rel_error <= rtol);
}
