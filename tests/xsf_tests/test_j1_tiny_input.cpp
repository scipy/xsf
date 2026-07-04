#include "../testing_utils.h"
#include <xsf/cephes/j1.h>

TEST_CASE("j1 left tail subnormal-input", "[j1][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    // Reference values computed with mpmath with 1000 digits of precision:
    //
    // import math
    // import mpmath as mp
    // mp.mp.dps = 1000
    // xs = np.logspace(np.log10(2**(-1074)), np.log10(np.finfo(np.float64).tiny), num=10)
    // for x in xs:
    //     print(x, mp.nstr(mp.besselj(1, x), 15))

    auto [x, ref] = GENERATE(
        test_case{5e-324, 0.0}, test_case{2.7e-322, 1.4e-322}, test_case{1.487e-320, 7.436e-321},
        test_case{8.159e-319, 4.0795e-319}, test_case{4.4763326e-317, 2.2381663e-317},
        test_case{2.4558778e-315, 1.2279389e-315}, test_case{1.3473833403e-313, 6.736916702e-314},
        test_case{7.39223207106e-312, 3.696116035527e-312}, test_case{4.0556457361725e-310, 2.02782286808625e-310},
        test_case{2.2250738585071875e-308, 1.112536929253594e-308}
    );
    const double w = xsf::cephes::j1(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rel_error);
    REQUIRE(rel_error <= 5e-16);
}

TEST_CASE("j1 left tail normal-input", "[j1][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    // Reference values computed with mpmath with 1000 digits of precision:
    //
    // import math
    // import mpmath as mp
    // mp.mp.dps = 1000
    // xs = np.logspace(-307, np.log10(np.nextafter(np.sqrt(np.finfo(np.float64).eps), np.inf)), 20, endpoint=True)
    // for x in xs:
    //     print(x, mp.nstr(mp.besselj(1, x), 15))

    auto [x, ref] = GENERATE(
        test_case{1e-307, 5e-308}, test_case{5.571330898367484e-292, 2.785665449183742e-292},
        test_case{3.1039727979104237e-276, 1.5519863989552119e-276},
        test_case{1.7293259556692778e-260, 8.646629778346389e-261},
        test_case{9.634647130169125e-245, 4.817323565084562e-245},
        test_case{5.367780725117885e-229, 2.6838903625589425e-229},
        test_case{2.990568260951069e-213, 1.4952841304755344e-213},
        test_case{1.6661445355914895e-197, 8.330722677957448e-198},
        test_case{9.282642532287614e-182, 4.641321266143807e-182},
        test_case{5.171667315863417e-166, 2.5858336579317085e-166},
        test_case{2.881306991294709e-150, 1.4406534956473544e-150},
        test_case{1.6052714668283512e-134, 8.026357334141756e-135},
        test_case{8.943498523408486e-119, 4.471749261704243e-119},
        test_case{4.982718966297293e-103, 2.4913594831486466e-103},
        test_case{2.77603761348138e-87, 1.38801880674069e-87}, test_case{1.5466224131020154e-71, 7.733112065510077e-72},
        test_case{8.616745238222938e-56, 4.308372619111469e-56},
        test_case{4.800673898907548e-40, 2.400336949453774e-40},
        test_case{2.674614282596992e-24, 1.337307141298496e-24},
        test_case{1.4901161193847666e-08, 7.450580596923833e-09}
    );
    const double w = xsf::cephes::j1(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rel_error);
    REQUIRE(rel_error <= 1e-15);
}
