#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("von_mises_cdf test", "[von_mises_cdf][xsf_tests]") {
    // Reference values computed with scipy.stats.von_mises_cdf
    // import numpy as np
    // from scipy.stats._stats import von_mises_cdf

    // np.set_printoptions(precision=20)
    // k_obj = np.linspace(1e-3, 3.0, 10)
    // x_obj = np.linspace(-10.0, 10.0, 10)
    // von_mises_cdf(k_obj, x_obj)

    SECTION("vector inputs (SciPy reference values)") {
        const size_t n = 10;
        const std::vector<double> k_obj = linspace(1e-3, 3.0, n);
        const std::vector<double> x_obj = linspace(-10.0, 10.0, n);
        const std::vector<double> expected = {-1.0914628654411138,   -0.7904352686647403, -0.3085050816322099,
                                              -0.008913110513925071, 0.1450905974469251,  0.8857870521643215,
                                              1.0018330872501384,    1.1579435355869516,  1.9789394168440955,
                                              2.0011107464769506};

        const double rtol = 1e-8;
        const std::vector<double> result = xsf::von_mises_cdf(k_obj, x_obj);
        for (size_t i = 0; i < n; ++i) {
            const double ref = expected[i];
            const auto rel_error = xsf::extended_relative_error(result[i], ref);
            CAPTURE(i, k_obj[i], x_obj[i], result[i], ref, rtol, rel_error);
            REQUIRE(rel_error <= rtol);
        }
    }

    SECTION("broadcast scalar k") {
        std::vector<double> k = {1.0};
        std::vector<double> x = linspace(-5.0, 5.0, 10);

        auto res = xsf::von_mises_cdf(k, x);
        REQUIRE(res.size() == x.size());
    }
}
