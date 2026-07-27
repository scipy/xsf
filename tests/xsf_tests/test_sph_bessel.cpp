#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/bessel.h>
#include <xsf/sph_bessel.h>

TEST_CASE("spherical_j reflection complex", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, std::complex<double>, std::complex<double>, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_j, rtol] = GENERATE(
        // Even n => j_n(-z) = j_n(z)   https://dlmf.nist.gov/10.47#E14
        test_case{10, {5., 0}, {0.00040734424424946041, 0}, 1e-14},
        test_case{10, {-5., 0}, {0.00040734424424946041, 0}, 1e-14},
        // Odd n => j_n(-z) = -j_n(z)   https://dlmf.nist.gov/10.47#E14
        test_case{7, {5., 0}, {0.017902778177989527, 0}, 1e-14},
        test_case{7, {-5., 0}, {-0.017902778177989527, 0}, 1e-14}
    );

    std::complex<double> result_spherical_j = xsf::sph_bessel_j(n, z);
    double rel_err_spherical_j = xsf::extended_relative_error(result_spherical_j, ref_spherical_j);

    CAPTURE(n, z, result_spherical_j, ref_spherical_j, rel_err_spherical_j, rtol);
    REQUIRE(rel_err_spherical_j <= rtol);
}

TEST_CASE("spherical_j_jac reflection derivative complex", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, std::complex<double>, std::complex<double>, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_j_jac, rtol] = GENERATE(
        // Even n => j_n'(-z) = -j_n'(z)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{10, {5., 0}, {0.00072194237812414795, 0}, 1e-14},
        test_case{10, {-5., 0}, {-0.00072194237812414795, 0}, 1e-14},
        // Odd n => j_n'(-z) = j_n'(z)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{7, {5., 0}, {0.01932245477463754818, 0}, 1e-14},
        test_case{7, {-5., 0}, {0.01932245477463754818, 0}, 1e-14}
    );

    std::complex<double> result_spherical_j_jac = xsf::sph_bessel_j_jac(n, z);
    double rel_err_spherical_j_jac = xsf::extended_relative_error(result_spherical_j_jac, ref_spherical_j_jac);

    CAPTURE(n, z, result_spherical_j_jac, ref_spherical_j_jac, rel_err_spherical_j_jac, rtol);
    REQUIRE(rel_err_spherical_j_jac <= rtol);
}

TEST_CASE("spherical_j real reflection", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, double, double, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_j, rtol] = GENERATE(
        // Even n => j_n(-x) = j_n(x)   https://dlmf.nist.gov/10.47#E14
        test_case{10, 5., 0.00040734424424946041, 1e-14}, test_case{10, -5., 0.00040734424424946041, 1e-14},
        // Odd n => j_n(-x) = -j_n(x)   https://dlmf.nist.gov/10.47#E14
        test_case{7, 5., 0.017902778177989527, 1e-14}, test_case{7, -5., -0.017902778177989527, 1e-14}
    );

    double result_spherical_j = xsf::sph_bessel_j(n, z);
    double rel_err_spherical_j = xsf::extended_relative_error(result_spherical_j, ref_spherical_j);

    CAPTURE(n, z, result_spherical_j, ref_spherical_j, rel_err_spherical_j, rtol);
    REQUIRE(rel_err_spherical_j <= rtol);
}

TEST_CASE("spherical_j_jac reflection derivative", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, double, double, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_j_jac, rtol] = GENERATE(
        // Even n => j_n'(-x) = -j_n'(x)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{10, 5., 0.00072194237812414795, 1e-14}, test_case{10, -5., -0.00072194237812414795, 1e-14},
        // Odd n => j_n'(-x) = j_n'(x)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{7, 5., 0.01932245477463754818, 1e-14}, test_case{7, -5., 0.01932245477463754818, 1e-14}
    );

    double result_spherical_j_jac = xsf::sph_bessel_j_jac(n, z);
    double rel_err_spherical_j_jac = xsf::extended_relative_error(result_spherical_j_jac, ref_spherical_j_jac);

    CAPTURE(n, z, result_spherical_j_jac, ref_spherical_j_jac, rel_err_spherical_j_jac, rtol);
    REQUIRE(rel_err_spherical_j_jac <= rtol);
}

TEST_CASE("spherical_y reflection complex", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, std::complex<double>, std::complex<double>, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_y, rtol] = GENERATE(
        // Even n => y_n(-z) = -y_n(z)   https://dlmf.nist.gov/10.47#E14
        test_case{10, {5., 0}, {-26.656114405718699, 0}, 1e-14},
        test_case{10, {-5., 0}, {26.656114405718699, 0}, 1e-14},
        // Odd n => y_n(-z) = y_n(z)   https://dlmf.nist.gov/10.47#E14
        test_case{7, {5., 0}, {-1.0273946388125984, 0}, 1e-14}, test_case{7, {-5., 0}, {-1.0273946388125984, 0}, 1e-14}
    );

    std::complex<double> result_spherical_y = xsf::sph_bessel_y(n, z);
    double rel_err_spherical_y = xsf::extended_relative_error(result_spherical_y, ref_spherical_y);

    CAPTURE(n, z, result_spherical_y, ref_spherical_y, rel_err_spherical_y, rtol);
    REQUIRE(rel_err_spherical_y <= rtol);
}

TEST_CASE("spherical_y_jac reflection derivative complex", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, std::complex<double>, std::complex<double>, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_y_jac, rtol] = GENERATE(
        // Even n => y_n'(-z) = y_n'(z)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{10, {5., 0}, {50.954006758163679, 0}, 1e-14}, test_case{10, {-5., 0}, {50.954006758163679, 0}, 1e-14},
        // Odd n => y_n'(-z) = -y_n'(z)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{7, {5., 0}, {1.1254238507300277, 0}, 1e-14}, test_case{7, {-5., 0}, {-1.1254238507300277, 0}, 1e-14}
    );

    std::complex<double> result_spherical_y_jac = xsf::sph_bessel_y_jac(n, z);
    double rel_err_spherical_y_jac = xsf::extended_relative_error(result_spherical_y_jac, ref_spherical_y_jac);

    CAPTURE(n, z, result_spherical_y_jac, ref_spherical_y_jac, rel_err_spherical_y_jac, rtol);
    REQUIRE(rel_err_spherical_y_jac <= rtol);
}

TEST_CASE("spherical_y real reflection", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, double, double, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_y, rtol] = GENERATE(
        // Even n => y_n(-x) = -y_n(x)   https://dlmf.nist.gov/10.47#E14
        test_case{10, 5., -26.656114405718699, 1e-14}, test_case{10, -5., 26.656114405718699, 1e-14},
        // Odd n => y_n(-x) = y_n(x)   https://dlmf.nist.gov/10.47#E14
        test_case{7, 5., -1.0273946388125984, 1e-14}, test_case{7, -5., -1.0273946388125984, 1e-14}
    );

    double result_spherical_y = xsf::sph_bessel_y(n, z);
    double rel_err_spherical_y = xsf::extended_relative_error(result_spherical_y, ref_spherical_y);

    CAPTURE(n, z, result_spherical_y, ref_spherical_y, rel_err_spherical_y, rtol);
    REQUIRE(rel_err_spherical_y <= rtol);
}

TEST_CASE("spherical_y_jac reflection derivative", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, double, double, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref_spherical_y_jac, rtol] = GENERATE(
        // Even n => y_n'(-x) = y_n'(x)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{10, 5., 50.954006758163679, 1e-14}, test_case{10, -5., 50.954006758163679, 1e-14},
        // Odd n => y_n'(-x) = -y_n'(x)   from DLMF 10.47.14 by taking the derivative with respect to z
        test_case{7, 5., 1.1254238507300277, 1e-14}, test_case{7, -5., -1.1254238507300277, 1e-14}
    );

    double result_spherical_y_jac = xsf::sph_bessel_y_jac(n, z);
    double rel_err_spherical_y_jac = xsf::extended_relative_error(result_spherical_y_jac, ref_spherical_y_jac);

    CAPTURE(n, z, result_spherical_y_jac, ref_spherical_y_jac, rel_err_spherical_y_jac, rtol);
    REQUIRE(rel_err_spherical_y_jac <= rtol);
}

TEST_CASE("spherical_i tiny inputs", "[spherical_bessel][xsf_tests]") {
    using test_case = std::tuple<long, double, double, double>;

    // Reference values computed with mpmath.
    auto [n, z, ref, rtol] = GENERATE(
        test_case{1, 1e-50, 3.333333333333333e-51, 1e-14}, test_case{100, 0.1, 7.463271153267418e-290, 1e-10},
        test_case{10, 5.2250558491838786e-27, 1.10313444760931e-273, 1e-14},
        test_case{3, 1e-50, 9.523809523809524e-153, 1e-14}, test_case{10, 1e-29, 7.273091945557419e-301, 1e-8},
        test_case{20, 1e-14, 7.625979004892142e-306, 1e-8}, test_case{200, 5, 3.168106195219311e-297, 1e-10}
    );

    double result = xsf::sph_bessel_i(n, z);
    double rel_err = xsf::extended_relative_error(result, ref);

    CAPTURE(n, z, result, ref, rel_err, rtol);
    REQUIRE(rel_err <= rtol);
}
