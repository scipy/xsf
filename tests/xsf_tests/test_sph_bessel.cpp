#include "../testing_utils.h"
#include <complex>
#include <tuple>
#include <xsf/bessel.h>
#include <xsf/sph_bessel.h>

TEST_CASE("spherical_j reflection complex", "[bessel][xsf_tests]") {
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

TEST_CASE("spherical_j_jac reflection derivative complex", "[bessel][xsf_tests]") {
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

TEST_CASE("spherical_j real reflection", "[bessel][xsf_tests]") {
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

TEST_CASE("spherical_j_jac reflection derivative", "[bessel][xsf_tests]") {
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
