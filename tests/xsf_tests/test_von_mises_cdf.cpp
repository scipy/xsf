#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("von Mises CDF test", "[von_mises_cdf][xsf_tests]") {
    // Reference values from scipy.stats._stats
    //
    // import numpy as np
    // from scipy.stats._stats import von_mises_cdf
    //
    // rng = np.random.default_rng(123456789)
    //
    // ks = rng.uniform(1e-5, 100.0, size=30)
    // xs = rng.uniform(-10.0, 10.0, size=30)
    //
    // print("von_mises_cdf(k, x) = y")
    // for k, x in zip(ks, xs):
    //     y = von_mises_cdf(k, x)
    //     print({k, x, y})
    SECTION("von Mises CDF scipy reference values") {
        using test_case = std::tuple<double, double, double>;
        auto [k, x, expected] = GENERATE(
            test_case{2.771283651124301, 3.035602639562689, 0.9997388651495406},
            test_case{90.67000648140171, 2.293446349343553, 1.0},
            test_case{88.13935665603788, -2.468954013494673, 2.122421301009252e-70},
            test_case{62.48973129311811, -3.681298253520844, 0.0},
            test_case{79.07148320264592, -8.504555095600582, -1.0},
            test_case{82.59080317722926, 7.475230407392111, 2.0},
            test_case{84.17058518158716, -8.427660362483639, -1.0},
            test_case{47.17280005458051, -7.298676164429338, -0.9999999999860489},
            test_case{95.72287841265161, -4.843563745488528, 0.0},
            test_case{94.65915383214558, 9.826110526459306, 2.0},
            test_case{52.985976440270356, 5.285135683873441, 1.000000000001819},
            test_case{7.977659725009496, 8.36402231256455, 1.9999882231205492},
            test_case{26.984432994165303, -2.984375644364212, 0.0},
            test_case{5.379993175825972, -2.1288745245032525, 6.371283853173759e-05},
            test_case{56.66185333254235, -0.03067047118158639, 0.4089111625414147},
            test_case{1.891783544480924, -9.162003142812924, -0.9969533342868221},
            test_case{20.53737003734821, 8.589497604841984, 2.0},
            test_case{24.93549047448925, -6.241188917880091, -0.41737059956329736},
            test_case{89.87780270746866, -2.2883249936323824, 6.995837213590395e-67},
            test_case{10.779738243151687, 4.826092542543064, 1.0000085304585589},
            test_case{74.86781378235544, 2.4434840647672407, 1.0},
            test_case{71.4703644287635, 0.9509492993052443, 0.19632196274422142},
            test_case{74.0779107888269, 5.7564413695173, 1.000003842629081},
            test_case{94.51174842697316, 8.489953936625268, 2.0},
            test_case{45.94609255068201, 7.60320354814775, 2.0},
            test_case{23.163504610347342, 0.25318137444421307, -0.1389481133079018},
            test_case{46.49738252955772, -7.629461259574195, -1.0},
            test_case{63.071055091009974, 7.753408442909997, 2.0},
            test_case{96.55498800547863, 5.051930923931616, 1.0},
            test_case{99.70821785882904, -8.449848608452879, -1.0}
        );
        constexpr double rtol = 1e-10;
        constexpr double atol = 1e-10;
        const double output = xsf::von_mises_cdf(k, x);
        CAPTURE(k, x, output, expected, rtol, atol);
        if (std::abs(expected) < 1e-10) {
            REQUIRE(std::abs(output - expected) <= atol);
        } else {
            REQUIRE(xsf::extended_relative_error(output, expected) <= rtol);
        }
    }

    SECTION("von Mises CDF normal approximation scipy reference values") {
        // Reference values from scipy.stats._stats
        //
        // from scipy.stats._stats import von_mises_cdf
        //
        // cases = [
        //     (50.0, -0.5),
        //     (50.0, 0.5),
        //     (75.0, 1.25),
        //     (100.0, -1.25),
        //     (100.0, 8.0),
        //     (250.0, -8.0),
        // ]
        //
        // print("von_mises_cdf(k, x) = y")
        // for k, x in cases:
        //     y = von_mises_cdf(k, x)
        //     print({k, x, y})
        using test_case = std::tuple<double, double, double>;
        auto [k, x, expected] = GENERATE(
            test_case{50, -0.5, 0.00024151470946535982}, test_case{50, 0.5, 0.99975848529053468},
            test_case{75, 1.25, 1}, test_case{100, -1.25, 7.4049490815546377e-32}, test_case{100, 8, 2},
            test_case{250, -8, -1}
        );
        constexpr double rtol = 1e-10;
        constexpr double atol = 1e-10;
        const double output = xsf::von_mises_cdf(k, x);
        CAPTURE(k, x, output, expected, rtol, atol);
        if (std::abs(expected) < 1e-10) {
            REQUIRE(std::abs(output - expected) <= atol);
        } else {
            REQUIRE(xsf::extended_relative_error(output, expected) <= rtol);
        }
    }
}
