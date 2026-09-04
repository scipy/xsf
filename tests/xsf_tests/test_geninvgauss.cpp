#include "../testing_utils.h"

#include <array>
#include <cmath>

#include <xsf/cpu/stats.h>

TEST_CASE("generalized inverse Gaussian density", "[geninvgauss][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.1/scipy/stats/tests/test_distributions.py#L941-L950
    constexpr std::array<double, 10> expected_pdf{
        2.081176820e-21, 4.488660034e-01, 3.747774338e-01, 2.693297528e-01, 1.905637275e-01,
        1.351476913e-01, 9.636538981e-02, 6.909040154e-02, 4.978006801e-02, 3.602084467e-02,
    };
    const std::vector<double> x = linspace(0.01, 5.0, 10);
    for (int i = 0; i < x.size(); ++i) {
        double pdf = xsf::cpu::geninvgauss_pdf(x[i], 0.5, 1.0);
        CAPTURE(i, x[i], pdf, expected_pdf[i]);
        REQUIRE(xsf::extended_relative_error(pdf, expected_pdf[i]) <= 1e-8);
    }
}

TEST_CASE("generalized inverse Gaussian PDF support", "[geninvgauss][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.1/scipy/stats/tests/test_distributions.py#L952-L957
    REQUIRE(xsf::cpu::geninvgauss_pdf(0.0, 0.5, 0.5) == 0.0);
    REQUIRE(xsf::cpu::geninvgauss_pdf(2e6, 50.0, 2.0) == 0.0);
}
