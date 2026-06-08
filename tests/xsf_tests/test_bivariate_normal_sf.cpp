#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("bivariate normal SF test", "[bivariate_normal_sf][xsf_tests]") {

    SECTION("bivariate normal SF infinite inputs") {
        // test cases with infinite inputs
        using test_case = std::tuple<double, double, double, double, double>;
        auto [dh, dk, r, expected, rtol] = GENERATE(
            test_case{std::numeric_limits<double>::infinity(), 0.0, 0.5, 0.0, 1e-13},
            test_case{
                -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), 0.5, 1.0, 1e-13
            },
            test_case{-std::numeric_limits<double>::infinity(), 1.0, 0.5, xsf::ndtr(-1.0), 1e-13},
            test_case{1.0, -std::numeric_limits<double>::infinity(), 0.5, xsf::ndtr(-1.0), 1e-13}
        );
        const double output = xsf::bivariate_normal_sf(dh, dk, r);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("bivariate normal SF analytical dh=0, dk=0") {
        // bivariate_normal_sf(0, 0, r) = 0.25 + asin(r)/(2*pi)
        double dh = 0.0;
        double dk = 0.0;
        double r = GENERATE(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0);
        double expected = 0.25 + std::asin(r) / (2 * M_PI);
        const double output = xsf::bivariate_normal_sf(dh, dk, r);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        double rtol = 1e-13;
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("bivariate normal SF r=0 independence") {
        // bivariate_normal_sf(h, k, 0) = ndtr(-h) * ndtr(-dk)
        double dh = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double dk = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double r = 0.0;
        const double output = xsf::bivariate_normal_sf(dh, dk, r);
        const double expected = xsf::ndtr(-dh) * xsf::ndtr(-dk);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        double rtol = 1e-13;
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("bivariate normal SF scipy reference values") {
        // Reference values computed with scipy.stats._stats_pythran._bvnu using:
        //
        // import numpy as np
        // from scipy.stats._stats_pythran import _bvnu

        // rng = np.random.default_rng(12345)

        // dhs = rng.uniform(-5, 5, 30)
        // dks = rng.uniform(-5, 5, 30)
        // rs = rng.uniform(-1.0, 1.0, 30)

        // for dh, dk, r in zip(dhs, dks, rs):
        //     out = _bvnu(dh, dk, r)
        //     print(f"bvnu({dh:.2f}, {dk:.2f}, {r}) = {out:.15e}")
        using test_case = std::tuple<double, double, double, double>;
        auto [dh, dk, r, expected] = GENERATE(
            test_case{
                -2.72663977532830337e+00, 3.54741904374001304e+00, -4.89535292360999463e-01, 1.63048185234028855e-04
            },
            test_case{
                -1.83241660290247133e+00, 1.01621241693713116e+00, -4.04869463703903865e-01, 1.39709320362846334e-01
            },
            test_case{
                2.97365457332734096e+00, 4.31988361135983467e+00, -4.41865760372466720e-01, 1.28031806796334057e-13
            },
            test_case{
                1.76254670750974562e+00, 2.24781361092020049e+00, -4.78841575017404875e-01, 4.39780764859199922e-06
            },
            test_case{
                -1.08890449398091005e+00, 3.60551317393292381e+00, -3.44768144013685163e-02, 1.29386382694497590e-04
            },
            test_case{
                -1.67186072133615471e+00, 4.29337801575316291e+00, -5.76041927296978784e-01, 1.15778545700046103e-06
            },
            test_case{
                9.83087535871898233e-01, 4.61860090823530411e-01, -8.73880665391868483e-03, 5.16611180443857371e-02
            },
            test_case{
                -3.13265814396286668e+00, 4.37672958767756892e+00, -5.07477348338524870e-01, 4.95977845167354134e-06
            },
            test_case{
                1.72756044014621302e+00, -5.01205992117570442e-02, 6.76965304933889511e-01, 4.09487936281593343e-02
            },
            test_case{
                4.41802865269937151e+00, -2.26226817510012523e+00, -6.39738819809929860e-01, 9.24117593428671459e-07
            },
            test_case{
                -2.51754285370429010e+00, -4.82212925252393232e-01, 7.24312583018472900e-01, 6.85090983008043941e-01
            },
            test_case{
                4.48881151833318270e+00, 1.65038923399530302e+00, -6.43401110309625102e-01, 2.62689218848210853e-15
            },
            test_case{
                1.67237453100372413e+00, -1.69109069532945355e+00, 5.01062663874488168e-01, 4.71807844828924308e-02
            },
            test_case{
                -4.04102064405887873e+00, 4.03454006808239107e+00, 2.22240807661130413e-01, 2.73546726966612863e-05
            },
            test_case{
                -5.81603338321872165e-01, -2.42925824723465666e+00, -5.81689930142785361e-01, 7.12049615861800822e-01
            },
            test_case{
                3.86479919327517685e+00, -1.60171662389680169e+00, 5.19744842247990446e-01, 5.55899147446801687e-05
            },
            test_case{
                1.97453499882022099e+00, -2.41146601357072665e+00, -5.01478860930175019e-01, 2.21962018865216500e-02
            },
            test_case{
                -1.73527135929887888e+00, -1.44553520055713980e+00, -8.28856536026884028e-01, 8.84499818506475077e-01
            },
            test_case{
                2.33928163330066496e+00, -4.94977666282868256e+00, 2.36113444636181891e-01, 9.66043125988507757e-03
            },
            test_case{
                -2.79865044454513789e+00, 1.28604544099678697e+00, 7.39366620646710881e-02, 9.90483577651285751e-02
            },
            test_case{
                -4.18405430457791905e+00, -2.17617292574881738e+00, 2.69053422430551414e-01, 9.85216691184765403e-01
            },
            test_case{
                -3.40104398924952456e+00, -4.31912310512054276e+00, -6.51251782617223496e-01, 9.99656522357890154e-01
            },
            test_case{
                -1.59899815045294713e+00, 1.16828977256380528e+00, -5.03671020287095095e-01, 9.76198757306667997e-02
            },
            test_case{
                -3.48068462979490567e-01, -3.23673679718796592e+00, 3.69645969278798336e-01, 6.36011782596151320e-01
            },
            test_case{
                -2.33578971709229055e+00, -1.95611612780410393e+00, -8.38256707498185039e-01, 9.65023160244897871e-01
            },
            test_case{
                3.15776403424806951e+00, -5.91131891238819707e-01, 7.50147201512252337e-01, 7.94919761218776865e-04
            },
            test_case{
                -3.06705610710505505e+00, -3.49797658937299172e+00, -1.42611236920016315e-01, 9.98684739814051370e-01
            },
            test_case{
                -3.70530923822799707e+00, -2.82071136914564979e+00, 2.36788390794755665e-01, 9.97501342993959095e-01
            },
            test_case{
                -4.08335248455064104e+00, -2.56668846664555028e-01, -3.73788991629760314e-01, 6.01261092417943566e-01
            },
            test_case{
                9.85680136649131899e-01, -2.36311449188081335e-01, -6.42074289414264809e-01, 2.99037436714108401e-02
            }
        );

        constexpr double rtol = 1e-10;
        constexpr double atol = 1e-10;
        const double output = xsf::bivariate_normal_sf(dh, dk, r);
        CAPTURE(dh, dk, r, output, expected, rtol, atol);
        CHECK(output >= 0.0);
        CHECK(output <= 1.0);
        if (std::abs(expected) < 1e-10) {
            REQUIRE(std::abs(output - expected) <= atol);
        } else {
            REQUIRE(xsf::extended_relative_error(output, expected) <= rtol);
        }
    }

    SECTION("bivariate normal SF symmetry h,k") {
        // bivariate_normal_sf(h, k, r) == bivariate_normal_sf(k, h, r)
        REQUIRE(xsf::bivariate_normal_sf(1.0, 2.0, 0.5) == xsf::bivariate_normal_sf(2.0, 1.0, 0.5));
    }

    SECTION("bivariate normal SF complementarity") {
        // bivariate_normal_sf(h, k, r) + bivariate_normal_sf(h, -k, -r) == ndtr(-h)
        double dh = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double dk = GENERATE(-1.0, 0.0, 1.0, 2.0);
        double r = GENERATE(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0);
        double output = xsf::bivariate_normal_sf(dh, dk, r) + xsf::bivariate_normal_sf(dh, -dk, -r);
        double expected = xsf::ndtr(-dh);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        double rtol = 1e-13;
        CAPTURE(dh, dk, r, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
