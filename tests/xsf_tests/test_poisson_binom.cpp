#include "../testing_utils.h"
#include <xsf/stats.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

TEST_CASE("poisson_binom_pmf test", "[poisson_binom_pmf][xsf_tests]") {
    // Reference values computed from scipy.stats._stats_pythran
    // import numpy as np
    // from scipy.stats._stats_pythran import _poisson_binom_pmf
    // ps = np.linspace(1e-5, 1 - 1e-5, 51)
    // expected = _poisson_binom_pmf(ps)

    const int n_trials = 51;
    const double start = 1e-5;
    const double end = 1.0 - 1e-5;
    std::vector<double> ps = linspace(start, end, n_trials);

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
    std::vector<double> results = std::vector(n_trials + 1, 0.0);

    std::mdspan p_span(ps.data(), ps.size());
    std::mdspan res_span(results.data(), results.size());

    xsf::poisson_binom_pmf_all(p_span, res_span);

    for (int i = 0; i < n_trials + 1; ++i) {
        const double rel_error = xsf::extended_relative_error(results[i], expected[i]);
        CAPTURE(i, results[i], expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("poisson_binom_cdf test", "[poisson_binom_cdf][xsf_tests]") {
    // Reference values computed as follows
    // import numpy as np
    // import itertools as it
    // from mpmath import mp
    // from scipy.stats._stats_pythran import _poisson_binom_pmf
    // ps = np.linspace(1e-5, 1 - 1e-5, 51)
    // expected_pmf = _poisson_binom_pmf(ps)
    // expected_cdf = list(it.accumulate([mp.mpf(p) for p in expected_pmf]))
    // expected_cdf = np.clip(expected_cdf, 0, 1)
    const int n_trials = 51;
    const double start = 1e-5;
    const double end = 1.0 - 1e-5;
    std::vector<double> ps = linspace(start, end, n_trials);

    const std::vector<double> expected = {
        3.42860361e-26, 3.43460107e-21, 6.03635984e-19, 4.68047694e-17, 2.16668202e-15, 6.80712079e-14, 1.55996609e-12,
        2.73188962e-11, 3.77720443e-10, 4.22399573e-09, 3.89224863e-08, 2.99927541e-07, 1.95608939e-06, 1.09055930e-05,
        5.24155149e-05, 2.18761065e-04, 7.97880640e-04, 2.55756944e-03, 7.24255193e-03, 1.82074669e-02, 4.08283984e-02,
        8.20567113e-02, 1.48559929e-01, 2.43633576e-01, 3.64218941e-01, 5.00000000e-01, 6.35781059e-01, 7.56366424e-01,
        8.51440071e-01, 9.17943289e-01, 9.59171602e-01, 9.81792533e-01, 9.92757448e-01, 9.97442431e-01, 9.99202119e-01,
        9.99781239e-01, 9.99947584e-01, 9.99989094e-01, 9.99998044e-01, 9.99999700e-01, 9.99999961e-01, 9.99999996e-01,
        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
        1.00000000e+00, 1.00000000e+00, 1.00000000e+00
    };

    const double rtol = 1e-8;
    std::vector<double> results = std::vector(n_trials + 1, 0.0);

    std::mdspan p_span(ps.data(), ps.size());
    std::mdspan res_span(results.data(), results.size());

    xsf::poisson_binom_cdf_all(p_span, res_span);

    for (int i = 0; i < n_trials + 1; ++i) {
        const double rel_error = xsf::extended_relative_error(results[i], expected[i]);
        CAPTURE(i, results[i], expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
