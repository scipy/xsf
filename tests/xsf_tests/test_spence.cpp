#include "../testing_utils.h"

#include <xsf/spence.h>

constexpr double rtol = 1e-14;

TEST_CASE("cspence special values test", "[cspence][xsf_tests]") {
    double phi = (1 + std::sqrt(5)) / 2;
    std::vector<double> z = {
        1, 2, 0.5, 0, -1, (-1 + std::sqrt(5)) / 2, (3 - std::sqrt(5)) / 2, phi, (3 + std::sqrt(5)) / 2
    };
    std::vector<std::complex<double>> expected = {
        0,
        -M_PI * M_PI / 12,
        M_PI * M_PI / 12 - std::log(2) * std::log(2) / 2,
        M_PI * M_PI / 6,
        std::complex<double>(M_PI * M_PI / 4, -M_PI * std::log(2)),
        std::complex<double>(M_PI * M_PI / 15 - std::log(phi) * std::log(phi), 0),
        std::complex<double>(M_PI * M_PI / 10 - std::log(phi) * std::log(phi), 0),
        std::complex<double>(-M_PI * M_PI / 15 + std::log(phi) * std::log(phi) / 2, 0),
        std::complex<double>(-M_PI * M_PI / 10 - std::log(phi) * std::log(phi), 0)
    };
    for (std::size_t i = 0; i < z.size(); ++i) {
        const std::complex<double> res = xsf::cspence(std::complex<double>(z[i]));
        const double rel_error = xsf::extended_relative_error(res, expected[i]);
        CAPTURE(i, z[i], res, expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("cspence nan", "[cspence][xsf_tests]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto w1 = xsf::cspence(std::complex<double>{nan, 0.0});
    REQUIRE(std::isnan(w1.real()));
    REQUIRE(std::isnan(w1.imag()));
}

TEST_CASE("cspence circle about z=1", "[cspence][xsf_tests]") {
    /*
    import mpmath
    import numpy as np
    mpmath.mp.dps = 50

    def spence(z):
        return complex(mpmath.polylog(2, 1 - z))

    r = np.linspace(0.5, 1.5, 4)
    theta = np.linspace(0, 2 * np.pi, 4)

    z_values = (1 + np.outer(r, np.exp(1j * theta))).flatten()

    inputs = []
    outputs = []

    for z in z_values:
        zc = complex(float(z.real), float(z.imag))
        yc = spence(z)
        yc = complex(float(yc.real), float(yc.imag))

        inputs.append(zc)
        outputs.append(yc)

    print("std::vector<std::complex<double>> z = {")
    for z in inputs:
        print(f"    {{{z.real:.17g}, {z.imag:.17g}}},")
    print("};\n")

    print("std::vector<std::complex<double>> expected = {")
    for z in outputs:
        print(f"    {{{z.real:.17g}, {z.imag:.17g}}},")
    print("};")
    */

    std::vector<std::complex<double>> z = {
        {1.5, 0},
        {0.75000000000000011, 0.43301270189221935},
        {0.74999999999999978, -0.43301270189221919},
        {1.5, -1.2246467991473532e-16},
        {1.8333333333333333, 0},
        {0.58333333333333348, 0.72168783648703216},
        {0.58333333333333304, -0.72168783648703194},
        {1.8333333333333333, -2.0410779985789218e-16},
        {2.1666666666666665, 0},
        {0.41666666666666696, 1.0103629710818451},
        {0.41666666666666619, -1.0103629710818447},
        {2.1666666666666665, -2.8575091980104902e-16},
        {2.5, 0},
        {0.25000000000000033, 1.299038105676658},
        {0.24999999999999933, -1.2990381056766576},
        {2.5, -3.6739403974420594e-16},
    };

    std::vector<std::complex<double>> expected = {
        {-0.4484142069236462, 0},
        {0.20399099867470188, -0.48285365695744431},
        {0.2039909986747023, 0.48285365695744425},
        {-0.4484142069236462, 9.9310309362119811e-17},
        {-0.7041493333253307, 6.5253044679985245e-54},
        {0.26684721128838046, -0.8383299706430849},
        {0.26684721128838101, 0.83832997064308479},
        {-0.7041493333253307, 1.4846045433819914e-16},
        {-0.93540928477883845, 3.2626522339992623e-55},
        {0.2676081604653901, -1.1866168564157165},
        {0.26760816046539104, 1.1866168564157167},
        {-0.93540928477883845, 1.8937690435164661e-16},
        {-1.1473806603755707, 0},
        {0.22237199660222873, -1.5088901672932802},
        {0.2223719966022297, 1.5088901672932804},
        {-1.1473806603755707, 2.2442650237561391e-16},
    };

    for (std::size_t i = 0; i < z.size(); ++i) {
        const std::complex<double> res = xsf::cspence(z[i]);
        const double rel_error = xsf::extended_relative_error(res, expected[i]);
        CAPTURE(i, z[i], res, expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
