#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("students t CDF test", "[students_t_cdf][xsf_tests]") {

    SECTION("students t CDF t = 0 special case") { REQUIRE(xsf::students_t_cdf(1.0, 0.0) == 0.5); }

    SECTION("students t CDF large t inputs") {
        // test cases with large inputs
        using test_case = std::tuple<double, double, double, double>;
        auto [nu, t, expected, rtol] = GENERATE(
            test_case{1.0, 1e200, 1.0, 1e-13}, test_case{1.0, -1e200, 0.0, 1e-13},
            test_case{1.0, std::numeric_limits<double>::infinity(), 1.0, 1e-13},
            test_case{1.0, -std::numeric_limits<double>::infinity(), 0.0, 1e-13}
        );
        const double output = xsf::students_t_cdf(nu, t);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(nu, t, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("students t CDF large degrees of freedom") {
        double nu = GENERATE(1e16, 1e100, std::numeric_limits<double>::infinity());
        double t = GENERATE(-5.0, -1.0, 1.0, 5.0);
        CAPTURE(nu, t);
        REQUIRE(xsf::students_t_cdf(nu, t) == xsf::ndtr(t));
    }

    SECTION("students t CDF scipy reference values") {
        // Reference values computed with scipy.special.stdtr using:
        //
        // import numpy as np
        // from scipy.special import stdtr

        // rng = np.random.default_rng(12345)

        // nus = rng.uniform(0, 30, 30)
        // ts = rng.uniform(-10, 10, 30)

        // for nu, t in zip(nus, ts):
        //     out = stdtr(nu, t)
        //     print(f"test_case{{{nu:.17e}, {t:.17e}, {out:.17e}}},")
        using test_case = std::tuple<double, double, double>;
        auto [nu, t, expected] = GENERATE(
            test_case{6.82008067401508988e+00, 7.09483808748002609e+00, 9.99890138723964594e-01},
            test_case{9.50275019129258602e+00, 2.03242483387426232e+00, 9.64503094798292659e-01},
            test_case{2.39209637199820229e+01, 8.63976722271966935e+00, 9.99999995954123477e-01},
            test_case{2.02876401225292362e+01, 4.49562722184040098e+00, 9.99892904079463074e-01},
            test_case{1.17332865180572696e+01, 7.21102634786584762e+00, 9.99993938039812624e-01},
            test_case{9.98441783599153609e+00, 8.58675603150632583e+00, 9.99996816400137156e-01},
            test_case{1.79492626076156938e+01, 9.23720181647060823e-01, 8.16058031053985111e-01},
            test_case{5.60202556811140084e+00, 8.75345917535513784e+00, 9.99910700863291702e-01},
            test_case{2.01826813204386397e+01, -1.00241198423514088e-01, 4.60570954832333901e-01},
            test_case{2.82540859580981163e+01, -4.52453635020025047e+00, 4.99188016405870876e-05},
            test_case{7.44737143888713060e+00, -9.64425850504786464e-01, 1.82549024492416806e-01},
            test_case{2.84664345549995481e+01, 3.30077846799060604e+00, 9.98699918568591527e-01},
            test_case{2.00171235930111706e+01, -3.38218139065890711e+00, 1.47890638351711313e-03},
            test_case{2.87693806782336248e+00, 8.06908013616478215e+00, 9.97687075268381518e-01},
            test_case{1.32551899850343844e+01, -4.85851649446931333e+00, 1.47786958685395566e-04},
            test_case{2.65943975798255323e+01, -3.20343324779360339e+00, 1.75504191225598532e-03},
            test_case{2.09236049964606643e+01, -4.82293202714145330e+00, 4.60009186664417870e-05},
            test_case{9.79418592210336314e+00, -2.89107040111427960e+00, 8.20406122949137416e-03},
            test_case{2.20178448999019949e+01, -9.89955332565736512e+00, 7.19158779738934737e-10},
            test_case{6.60404866636458721e+00, 2.57209088199357394e+00, 9.80604786750692625e-01},
            test_case{2.44783708626624374e+00, -4.35234585149763475e+00, 1.67989748330685286e-02},
            test_case{4.79686803225142544e+00, -8.63824621024108552e+00, 2.10799420311352358e-04},
            test_case{1.02030055486411584e+01, 2.33657954512761057e+00, 9.79448065462404593e-01},
            test_case{1.39557946110615276e+01, -6.47347359437593184e+00, 7.43521811249499339e-06},
            test_case{7.99263084872312923e+00, -3.91223225560820786e+00, 2.23751923309389172e-03},
            test_case{2.44732921027442103e+01, -1.18226378247763941e+00, 1.24230448248389705e-01},
            test_case{5.79883167868483440e+00, -6.99595317874598344e+00, 2.46391679471401283e-04},
            test_case{3.88407228531600790e+00, -5.64142273829129959e+00, 2.65310907635721380e-03},
            test_case{2.74994254634807778e+00, -5.13337693329110056e-01, 3.23041117900401020e-01},
            test_case{1.79570404099473961e+01, -4.72622898376162670e-01, 3.21087841259371398e-01}
        );

        constexpr double rtol = 1e-12;
        const double output = xsf::students_t_cdf(nu, t);
        CAPTURE(nu, t, output, expected, rtol);
        CHECK(output >= 0.0);
        CHECK(output <= 1.0);
        REQUIRE(xsf::extended_relative_error(output, expected) <= rtol);
    }
}
