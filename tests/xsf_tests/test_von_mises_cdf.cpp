#include "../testing_utils.h"
#include <cmath>
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
            test_case{95.72287841265161, -4.843563745488528, 0.0}, test_case{94.65915383214558, 9.826110526459306, 2.0},
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
            test_case{94.51174842697316, 8.489953936625268, 2.0}, test_case{45.94609255068201, 7.60320354814775, 2.0},
            test_case{23.163504610347342, 0.25318137444421307, -0.1389481133079018},
            test_case{46.49738252955772, -7.629461259574195, -1.0},
            test_case{63.071055091009974, 7.753408442909997, 2.0}, test_case{96.55498800547863, 5.051930923931616, 1.0},
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

    SECTION("von Mises CDF periodic extension") {
        // The CDF is extended periodically over x, adding one CDF period per
        // 2*pi shift.
        using test_case = std::tuple<double, double, int>;
        auto [k, x, periods] = GENERATE(
            test_case{0.25, 0.75, 1}, test_case{3.0, -1.25, -1}, test_case{49.999, 2.0, 2}, test_case{50.0, -0.5, -2}
        );
        const double expected = xsf::von_mises_cdf(k, x) + periods;
        const double output = xsf::von_mises_cdf(k, x + periods * 2.0 * M_PI);
        CAPTURE(k, x, periods, output, expected);
        REQUIRE(std::abs(output - expected) <= 1e-12);
    }

    SECTION("von Mises CDF scipy.stats periodic distribution behavior") {
        // Port of scipy.stats.tests.test_distributions.TestVonMises.
        //
        // @pytest.mark.parametrize('k', [0.1, 1, 101])
        // @pytest.mark.parametrize('x', [0, 1, np.pi, 10, 100])
        // def test_vonmises_periodic(self, k, x):
        //     def check_vonmises_cdf_periodic(k, L, s, x):
        //         vm = stats.vonmises(k, loc=L, scale=s)
        //         assert_almost_equal(vm.cdf(x) % 1,
        //                             vm.cdf(x % (2 * np.pi * s)) % 1)
        //
        //     check_vonmises_cdf_periodic(k, 0, 1, x)
        //     check_vonmises_cdf_periodic(k, 1, 1, x)
        //     check_vonmises_cdf_periodic(k, 0, 10, x)
        auto k = GENERATE(0.1, 1.0, 101.0);
        auto x = GENERATE(0.0, 1.0, M_PI, 10.0, 100.0);
        auto [loc, scale] = GENERATE(
            std::tuple<double, double>{0.0, 1.0}, std::tuple<double, double>{1.0, 1.0},
            std::tuple<double, double>{0.0, 10.0}
        );

        auto modulo_one = [](double value) {
            const double result = std::fmod(value, 1.0);
            return result < 0.0 ? result + 1.0 : result;
        };
        const double period = 2.0 * M_PI * scale;
        const double shifted_x = (x - loc) / scale;
        const double shifted_wrapped_x = (std::fmod(x, period) - loc) / scale;
        const double output = modulo_one(xsf::von_mises_cdf(k, shifted_x));
        const double expected = modulo_one(xsf::von_mises_cdf(k, shifted_wrapped_x));
        CAPTURE(k, x, loc, scale, output, expected);
        REQUIRE(std::abs(output - expected) <= 1e-7);
    }

    SECTION("von Mises CDF uniform distribution at k = 0") {
        // At zero concentration the von Mises distribution is uniform on the
        // circle.
        auto x = GENERATE(-M_PI_2, 0.0, M_PI_2);
        const double expected = 0.5 + x / (2.0 * M_PI);
        const double output = xsf::von_mises_cdf(0.0, x);
        CAPTURE(x, output, expected);
        REQUIRE(std::abs(output - expected) <= 1e-14);
    }

    SECTION("von Mises CDF branch boundary at x = 0") {
        // k = 50 switches from the series implementation to the normal
        // approximation.
        auto k = GENERATE(49.999, 50.0);
        const double output = xsf::von_mises_cdf(k, 0.0);
        CAPTURE(k, output);
        REQUIRE(std::abs(output - 0.5) <= 1e-14);
    }

    SECTION("von Mises CDF scipy.stats large-k numerical behavior") {
        // Port of scipy.stats.tests.test_distributions.TestVonMises.
        //
        // def test_vonmises_numerical(self):
        //     vm = stats.vonmises(800)
        //     assert_almost_equal(vm.cdf(0), 0.5)
        const double output = xsf::von_mises_cdf(800.0, 0.0);
        CAPTURE(output);
        REQUIRE(std::abs(output - 0.5) <= 1e-7);
    }

    SECTION("von Mises CDF float overload") {
        // The float overload delegates through the double implementation.
        auto [k, x] = GENERATE(
            std::tuple<float, float>{0.5F, -1.25F}, std::tuple<float, float>{10.0F, 2.5F},
            std::tuple<float, float>{50.0F, 0.5F}
        );
        const float output = xsf::von_mises_cdf(k, x);
        const float expected = static_cast<float>(xsf::von_mises_cdf(static_cast<double>(k), static_cast<double>(x)));
        CAPTURE(k, x, output, expected);
        REQUIRE(output == expected);
    }
}
