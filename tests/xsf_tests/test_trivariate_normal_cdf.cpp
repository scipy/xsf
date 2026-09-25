#include "../testing_utils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

#include <catch2/catch_approx.hpp>

#include <xsf/multivariate_normal.h>
#include <xsf/stats.h>

TEST_CASE("trivariate normal CDF test", "[trivariate_normal_cdf][xsf_tests]") {
    SECTION("trivariate normal CDF known SciPy value") {
        // Port of scipy.stats.tests.test_multivariate.TestMultivariateNormal.test_cdf_known
        // for ndim=3. For an equicorrelation matrix with rho=0.5,
        // P(X < 0, Y < 0, Z < 0) = 1 / (1 + ndim) = 0.25.
        const double output = xsf::trivariate_normal_cdf(0.0, 0.0, 0.0, 0.5, 0.5, 0.5);
        CAPTURE(output);
        REQUIRE(std::abs(output - 0.25) <= 1e-14);
    }

    SECTION("trivariate normal CDF infinite inputs") {
        // These expected values follow from the definition of a CDF. A -inf
        // upper limit makes the probability zero, three +inf limits give one,
        // and one +inf limit removes that variable. The remaining bivariate
        // CDF is evaluated as bivariate_normal_cdf(h1, h2, r).
        using test_case = std::tuple<double, double, double, double, double, double, double, double>;
        auto [dh, dk, dl, r_xy, r_xz, r_yz, expected, rtol] = GENERATE(
            test_case{-std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.2, -0.3, 0.4, 0.0, 1e-13},
            test_case{
                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(), 0.2, -0.3, 0.4, 1.0, 1e-13
            },
            test_case{
                std::numeric_limits<double>::infinity(), 1.0, 2.0, 0.2, -0.3, 0.4,
                xsf::bivariate_normal_cdf(1.0, 2.0, 0.4), 1e-13
            },
            test_case{
                1.0, std::numeric_limits<double>::infinity(), 2.0, 0.2, -0.3, 0.4,
                xsf::bivariate_normal_cdf(1.0, 2.0, -0.3), 1e-13
            },
            test_case{
                1.0, 2.0, std::numeric_limits<double>::infinity(), 0.2, -0.3, 0.4,
                xsf::bivariate_normal_cdf(1.0, 2.0, 0.2), 1e-13
            }
        );
        const double output = xsf::trivariate_normal_cdf(dh, dk, dl, r_xy, r_xz, r_yz);
        CAPTURE(dh, dk, dl, r_xy, r_xz, r_yz);
        REQUIRE(output == Catch::Approx(expected).epsilon(rtol).margin(1e-15));
    }

    SECTION("trivariate normal analytical value at zero") {
        // The expected values use the trivariate normal orthant formula
        //   P(X <= 0, Y <= 0, Z <= 0)
        //       = (1 + 2*(asin(r_xy) + asin(r_xz) + asin(r_yz))/pi)/8.
        using test_case = std::tuple<double, double, double, double>;
        auto [r_xy, r_xz, r_yz, expected] = GENERATE(
            test_case{0.9003169287616585, 0.9663883660249392, 0.8452039768406975, 0.39860596161542950},
            test_case{0.24761277161200018, 0.8986247023063125, 0.3305664724164368, 0.26057962506858890},
            test_case{-0.09236368828963112, 0.47195045831645027, 0.7576654495126486, 0.22516693852411820},
            test_case{-0.10055444141790104, 0.44079466125979155, -0.2947936519629224, 0.12949768309151660},
            test_case{0.1065473509416631, 0.2969282597497007, -0.4265042580475209, 0.12242167847560080},
            test_case{-0.9961097121761741, -0.8465939176143792, 0.8057658811126734, 0.00124546072254073},
            test_case{0.7616901423134714, -0.2812967181610399, -0.8112909809942859, 0.09590771689824270},
            test_case{0.3325635658924687, -0.4511411429898154, -0.8177101437857275, 0.03854347572563056},
            test_case{0.18346562422397583, 0.9873586596792112, 0.19912700304811634, 0.26796895124475840},
            test_case{-0.03493338924713792, 0.3200305655432352, -0.1970746854218671, 0.13235678304407310}
        );
        const double output = xsf::trivariate_normal_cdf(0.0, 0.0, 0.0, r_xy, r_xz, r_yz);
        CAPTURE(r_xy, r_xz, r_yz);
        REQUIRE(output == Catch::Approx(expected).epsilon(1e-13).margin(1e-14));
    }

    SECTION("trivariate normal independent variables") {
        // With zero correlations the variables are independent, so the joint
        // CDF is the product of the three univariate standard normal CDFs.
        const double h = GENERATE(-1.0, 0.0, 1.0);
        const double k = GENERATE(-0.5, 0.5);
        const double l = 0.25;
        const double cdf = xsf::trivariate_normal_cdf(h, k, l, 0.0, 0.0, 0.0);
        const double expected_cdf = xsf::ndtr(h) * xsf::ndtr(k) * xsf::ndtr(l);
        CAPTURE(h, k, l, cdf, expected_cdf);
        REQUIRE(xsf::extended_relative_error(cdf, expected_cdf) <= 1e-14);
    }

    SECTION("trivariate normal invalid inputs") {
        // NaN input, a correlation outside [-1, 1], and a correlation matrix
        // with negative determinant are outside the function's domain. The
        // expected value in each case is therefore NaN.
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        REQUIRE(std::isnan(xsf::trivariate_normal_cdf(nan, 0.0, 0.0, 0.0, 0.0, 0.0)));
        REQUIRE(std::isnan(xsf::trivariate_normal_cdf(0.0, 0.0, 0.0, 1.1, 0.0, 0.0)));
        // All pairwise correlations are in range, but the matrix is not
        // positive semidefinite.
        REQUIRE(std::isnan(xsf::trivariate_normal_cdf(0.0, 0.0, 0.0, 0.9, 0.9, -0.9)));
    }

    SECTION("trivariate normal singular correlations") {
        // Z = Y or Z = -Y reduces the event to a bivariate probability.
        // rho=0.3 also exercises roundoff in the determinant at singularity.
        const double rho = GENERATE(0.0, 0.3, -0.3, 1.0, -1.0);
        const double h1 = 0.4;
        const double h2 = 0.7;
        const double h3 = GENERATE(-0.9, -0.7, 0.2);
        const double sign = GENERATE(-1.0, 1.0);
        const double expected =
            sign > 0.0
                ? xsf::bivariate_normal_cdf(h1, std::min(h2, h3), rho)
                : (h2 > -h3 ? xsf::bivariate_normal_cdf(h1, h2, rho) - xsf::bivariate_normal_cdf(h1, -h3, rho) : 0.0);
        CAPTURE(rho, h3, sign);
        REQUIRE(xsf::trivariate_normal_cdf(h1, h2, h3, rho, sign * rho, sign) == Catch::Approx(expected).margin(1e-14));
        REQUIRE(xsf::trivariate_normal_cdf(h2, h1, h3, rho, sign, sign * rho) == Catch::Approx(expected).margin(1e-14));
        REQUIRE(xsf::trivariate_normal_cdf(h2, h3, h1, sign, rho, sign * rho) == Catch::Approx(expected).margin(1e-14));
    }

    SECTION("trivariate normal CDF reference values") {
        // Each expected value was computed in Mathematica by evaluating
        //   CDF[MultinormalDistribution[{0, 0, 0},
        //       {{1, r_xy, r_xz}, {r_xy, 1, r_yz}, {r_xz, r_yz, 1}}],
        //       {x, y, z}]
        // and recording the result to the displayed decimal precision.
        using test_case = std::tuple<double, double, double, double, double, double, double>;
        auto [r_xy, r_xz, r_yz, x, y, z, expected] = GENERATE(
            test_case{
                0.9003169287616585, 0.9663883660249392, 0.8452039768406975, 1.95324484, -0.53790441, -0.30868175,
                0.24865818416840870
            },
            test_case{
                0.24761277161200018, 0.8986247023063125, 0.3305664724164368, 0.61305867, -0.43716194, 1.24196053,
                0.26996099541815140
            },
            test_case{
                -0.09236368828963112, 0.47195045831645027, 0.7576654495126486, -0.9064342, -1.76571409, -0.13530755,
                0.00516301442987467
            },
            test_case{
                -0.10055444141790104, 0.44079466125979155, -0.2947936519629224, -0.52343401, 0.04762014, 0.34583066,
                0.10778223248220700
            },
            test_case{
                0.1065473509416631, 0.2969282597497007, -0.4265042580475209, -0.66472225, 0.41138206, 0.61340843,
                0.14361434787732490
            },
            test_case{
                -0.9961097121761741, -0.8465939176143792, 0.8057658811126734, -0.67022528, 0.28551754, -0.07723174,
                8.77076189453874e-14
            },
            test_case{
                0.7616901423134714, -0.2812967181610399, -0.8112909809942859, -2.9621344, 0.50887025, 2.51361749,
                0.00145384520124692
            },
            test_case{
                0.3325635658924687, -0.4511411429898154, -0.8177101437857275, -0.52819329, -0.08140761, 0.78034451,
                0.07738570033617177
            },
            test_case{
                0.18346562422397583, 0.9873586596792112, 0.19912700304811634, -1.39884246, 0.36286964, 0.56865498,
                0.06167444999671867
            },
            test_case{
                -0.03493338924713792, 0.3200305655432352, -0.1970746854218671, 1.34788001, 0.43214523, -0.07194596,
                0.27254103329345060
            }
        );
        const double output = xsf::trivariate_normal_cdf(x, y, z, r_xy, r_xz, r_yz, 1e-10);
        CAPTURE(x, y, z, r_xy, r_xz, r_yz);
        // epsi specifies absolute accuracy, including for tiny tail probabilities.
        REQUIRE(output == Catch::Approx(expected).epsilon(0.0).margin(1e-10));
        // All permutations describe the same event and exercise both sorting steps.
        REQUIRE(
            xsf::trivariate_normal_cdf(x, z, y, r_xz, r_xy, r_yz, 1e-10) ==
            Catch::Approx(expected).epsilon(0.0).margin(1e-10)
        );
        REQUIRE(
            xsf::trivariate_normal_cdf(y, x, z, r_xy, r_yz, r_xz, 1e-10) ==
            Catch::Approx(expected).epsilon(0.0).margin(1e-10)
        );
        REQUIRE(
            xsf::trivariate_normal_cdf(y, z, x, r_yz, r_xy, r_xz, 1e-10) ==
            Catch::Approx(expected).epsilon(0.0).margin(1e-10)
        );
        REQUIRE(
            xsf::trivariate_normal_cdf(z, x, y, r_xz, r_yz, r_xy, 1e-10) ==
            Catch::Approx(expected).epsilon(0.0).margin(1e-10)
        );
        REQUIRE(
            xsf::trivariate_normal_cdf(z, y, x, r_yz, r_xz, r_xy, 1e-10) ==
            Catch::Approx(expected).epsilon(0.0).margin(1e-10)
        );
    }
}
