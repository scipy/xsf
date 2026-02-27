#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("poisson_binom_pmf test", "[poisson_binom_pmf][xsf_tests]") {
    // Reference values computed from scipy.stats._stats_pythran
    // import numpy as np
    // from scipy.stats._stats_pythran import _poisson_binom_pmf
    // ps = np.linspace(1e-5, 1 - 1e-5, 51)
    // expected = _poisson_binom_pmf(ps)

    const int n_points = 51;
    const double start = 1e-5;
    const double end = 1.0 - 1e-5;
    std::vector<double> ps(n_points);
    for (int i = 0; i < n_points; ++i) {
        ps[i] = start + (end - start) * i / (n_points - 1);
    }

    const std::vector<double> expected = {
        3.42860361e-26, 3.43456678e-21, 6.00201383e-19, 4.62011334e-17, 2.11987726e-15, 6.59045259e-14, 1.49189489e-12,
        2.57589301e-11, 3.50401547e-10, 3.84627529e-09, 3.46984905e-08, 2.61005054e-07, 1.65616184e-06, 8.94950358e-06,
        4.15099219e-05, 1.66345550e-04, 5.79119575e-04, 1.75968880e-03, 4.68498249e-03, 1.09649150e-02, 2.26209315e-02,
        4.12283129e-02, 6.65032178e-02, 9.50736474e-02, 1.20585365e-01, 1.35781059e-01, 1.35781059e-01, 1.20585365e-01,
        9.50736474e-02, 6.65032178e-02, 4.12283129e-02, 2.26209315e-02, 1.09649150e-02, 4.68498249e-03, 1.75968880e-03,
        5.79119575e-04, 1.66345550e-04, 4.15099219e-05, 8.94950358e-06, 1.65616184e-06, 2.61005054e-07, 3.46984905e-08,
        3.84627529e-09, 3.50401547e-10, 2.57589301e-11, 1.49189489e-12, 6.59045259e-14, 2.11987726e-15, 4.62011334e-17,
        6.00201383e-19, 3.43456678e-21, 3.42860361e-26
    };

    const double rtol = 1e-8;
    const std::vector<double> results = xsf::poisson_binom_pmf(ps);

    for (int i = 0; i < n_points; ++i) {
        const double rel_error = xsf::extended_relative_error(results[i], expected[i]);
        CAPTURE(i, ps[i], results[i], expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("poisson_binom test", "[poisson_binom][xsf_tests]") {
    // Reference values computed from scipy.stats._stats_pythran
    // import numpy as np
    // from scipy.stats._stats_pythran import _poisson_binom_pmf
    // args = np.array([[0.1], [0.5]])
    // for mode in ["pdf", "cdf"]:
    //     for k in [0, 1, 2]:
    //         print(_poisson_binom(np.array([k]), args, mode))
    // args = np.array([[0.1, 0.2, 0.3], [0.5, 0.6, 0.7]])
    // k = np.array([0, 0, 2])
    // for mode in ["pdf", "cdf"]:
    //     print(_poisson_binom(k, args, mode))
    using test_case =
        std::tuple<std::vector<std::vector<double>>, std::vector<int>, std::vector<double>, std::string, double>;

    auto [args, k, expected, mode, rtol] = GENERATE(
        test_case{{{0.1}, {0.5}}, {0}, {0.45}, "pdf", 1e-12}, test_case{{{0.1}, {0.5}}, {0}, {0.45}, "cdf", 1e-12},
        test_case{{{0.1}, {0.5}}, {1}, {0.5}, "pdf", 1e-12}, test_case{{{0.1}, {0.5}}, {1}, {0.95}, "cdf", 1e-12},
        test_case{{{0.1}, {0.5}}, {2}, {0.05}, "pdf", 1e-12}, test_case{{{0.1}, {0.5}}, {2}, {1.0}, "cdf", 1e-12},
        test_case{{{0.1, 0.2, 0.3}, {0.5, 0.6, 0.7}}, {0, 0, 2}, {0.45, 0.32, 0.21}, "pdf", 1e-12},
        test_case{{{0.1, 0.2, 0.3}, {0.5, 0.6, 0.7}}, {0, 0, 2}, {0.45, 0.32, 1.0}, "cdf", 1e-12}
    );

    const std::vector<double> results = xsf::poisson_binom(k, args, mode);

    for (size_t i = 0; i < expected.size(); ++i) {
        const double rel_error = xsf::extended_relative_error(results[i], expected[i]);
        CAPTURE(i, results[i], expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
