#include "../testing_utils.h"
#include <xsf/stats.h>
#include <xsf/third_party/kokkos/mdspan.hpp>

#include <limits>

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

TEST_CASE("poisson_binom_sf test", "[poisson_binom_sf][xsf_tests]") {
    // Reference values computed exactly with rational arithmetic
    // from fractions import Fraction
    // import numpy as np
    // ps = np.linspace(1e-5, 1 - 1e-5, 51)
    // dist = [Fraction(1)]
    // for p in map(Fraction, ps):
    //     dist = [a * (1 - p) + b * p for a, b in zip(dist + [0], [0] + dist)]
    // expected = [float(sum(dist[k + 1:])) for k in range(52)]
    const int n_trials = 51;
    std::vector<double> ps = linspace(1e-5, 1.0 - 1e-5, n_trials);

    const std::vector<double> expected = {
        1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
        9.9999999999999778e-01, 9.9999999999993194e-01, 9.9999999999844003e-01, 9.9999999997268108e-01,
        9.9999999962227959e-01, 9.9999999577600429e-01, 9.9999996107751377e-01, 9.9999970007245931e-01,
        9.9999804391061464e-01, 9.9998909440703476e-01, 9.9994758448514331e-01, 9.9978123893485937e-01,
        9.9920211936023795e-01, 9.9744243055810211e-01, 9.9275744806906907e-01, 9.8179253310320613e-01,
        9.5917160163748016e-01, 9.1794328873788755e-01, 8.5144007092119833e-01, 7.5636642350749517e-01,
        6.3578105865804135e-01, 5.0000000000000033e-01, 3.6421894134195937e-01, 2.4363357649250547e-01,
        1.4855992907880214e-01, 8.2056711262112725e-02, 4.0828398362520030e-02, 1.8207466896793922e-02,
        7.2425519309309737e-03, 2.5575694418978781e-03, 7.9788063976207594e-04, 2.1876106514061484e-04,
        5.2415514856685807e-05, 1.0905592965264153e-05, 1.9560893853878383e-06, 2.9992754067288714e-07,
        3.8922486276778966e-08, 4.2239957318624452e-09, 3.7772044305291541e-10, 2.7318896193863000e-11,
        1.5599660949304700e-12, 6.8071207937338763e-14, 2.1666820248981403e-15, 4.6804769402516302e-17,
        6.0363598439614694e-19, 3.4346010668113848e-21, 3.4286036110438167e-26, 0.0000000000000000e+00
    };

    const double rtol = 1e-14;
    std::vector<double> results(n_trials + 1, 0.0);

    std::mdspan p_span(ps.data(), ps.size());
    std::mdspan res_span(results.data(), results.size());

    xsf::poisson_binom_sf_all(p_span, res_span);

    for (int i = 0; i < n_trials + 1; ++i) {
        const double rel_error = xsf::extended_relative_error(results[i], expected[i]);
        CAPTURE(i, results[i], expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("poisson_binom_logsf test", "[poisson_binom_logsf][xsf_tests]") {
    // Logarithms of the exact values in the sf test above, taken with
    // decimal.Decimal at 80 digits: float((Decimal(x.numerator) / Decimal(x.denominator)).ln())
    const int n_trials = 51;
    std::vector<double> ps = linspace(1e-5, 1.0 - 1e-5, n_trials);

    const std::vector<double> expected = {-3.4286036110281392e-26, -3.4346010668112833e-21,
                                          -6.0363598439613279e-19, -4.6804769402515328e-17,
                                          -2.1666820248981016e-15, -6.8071207937339899e-14,
                                          -1.5599660949316617e-12, -2.7318896194235758e-11,
                                          -3.7772044312424661e-10, -4.2239957407834616e-09,
                                          -3.8922487034258501e-08, -2.9992758565115773e-07,
                                          -1.9560912985331558e-06, -1.0905652431675361e-05,
                                          -5.2416888597787997e-05, -2.1878499683270253e-04,
                                          -7.9819911593497427e-04, -2.5608456098331045e-03,
                                          -7.2689065370015633e-03, -1.8375262700605008e-02,
                                          -4.1685282006648838e-02, -8.5619667243901118e-02,
                                          -1.6082616195376029e-01, -2.7922933304626785e-01,
                                          -4.5290102230176743e-01, -6.9314718055994451e-01,
                                          -1.0100001049631351e+00, -1.4120899181106070e+00,
                                          -1.9067668393287394e+00, -2.5003446700060059e+00,
                                          -3.1983774013972726e+00, -4.0059234999821127e+00,
                                          -4.9277816580418881e+00, -5.9686979082498732e+00,
                                          -7.1335515459347940e+00, -8.4275304505058983e+00,
                                          -9.8563079254105510e+00, -1.1426234784248937e+01,
                                          -1.3144563289344731e+01, -1.5019724922553918e+01,
                                          -1.7061693699942879e+01, -1.9282484300992191e+01,
                                          -2.1696866762567844e+01, -2.4323442484902241e+01,
                                          -2.7186357028922220e+01, -3.0318222062091326e+01,
                                          -3.3765579418287800e+01, -3.7600546565851545e+01,
                                          -4.1951315611462441e+01, -4.7120386171700211e+01,
                                          -5.8635059350042162e+01, -std::numeric_limits<double>::infinity()};

    const double rtol = 1e-14;
    std::vector<double> results(n_trials + 1, 0.0);

    std::mdspan p_span(ps.data(), ps.size());
    std::mdspan res_span(results.data(), results.size());

    xsf::poisson_binom_logsf_all(p_span, res_span);

    for (int i = 0; i < n_trials + 1; ++i) {
        const double rel_error = xsf::extended_relative_error(results[i], expected[i]);
        CAPTURE(i, results[i], expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("poisson_binom small upper tails", "[poisson_binom_sf][poisson_binom_logsf][xsf_tests]") {
    // scipy/scipy#26409. P(X > k) and its logarithm, computed as in the two tests above.
    struct Case {
        std::vector<double> p;
        int k;
        double sf;
        double logsf;
    };
    std::vector<Case> cases;
    std::vector<double> p1(32);
    for (int i = 0; i < 32; ++i) {
        p1[i] = (1 + i % 9) / 128.0;
    }
    cases.push_back({p1, 23, 1.1963531710514960e-29, -6.6595689791368116e+01});
    cases.push_back({std::vector<double>(256, 1.0 / 128), 249, 0.0, -1.1864216412494768e+03});
    cases.push_back({std::vector<double>(16, 1e-20), 0, 1.5999999999999999e-19, -4.3279113137641133e+01});

    const double rtol = 1e-14;
    for (auto &c : cases) {
        std::vector<double> sf(c.p.size() + 1), logsf(c.p.size() + 1);
        std::mdspan p_span(c.p.data(), c.p.size());
        xsf::poisson_binom_sf_all(p_span, std::mdspan(sf.data(), sf.size()));
        xsf::poisson_binom_logsf_all(p_span, std::mdspan(logsf.data(), logsf.size()));

        const double sf_error = xsf::extended_relative_error(sf[c.k], c.sf);
        const double logsf_error = xsf::extended_relative_error(logsf[c.k], c.logsf);
        CAPTURE(c.p.size(), c.k, sf[c.k], c.sf, logsf[c.k], c.logsf, sf_error, logsf_error);
        REQUIRE(sf_error <= rtol);
        REQUIRE(logsf_error <= rtol);
    }
}
