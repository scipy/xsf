#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("cvm_cdf_inf test", "[cvm_cdf_inf][xsf_tests]") {
    // Reference values computed with scipy.stats._hypotests._cdf_cvm_inf
    // from scipy.stats._hypotests import _cdf_cvm_inf
    // xs = np.linspace(2e-3, 1-2e-3, 51)
    // expected = _cdf_cvm_inf(xs)

    const int n_points = 51;
    const double start = 2e-3;
    const double end = 1.0 - 2e-3;
    std::vector<double> xs(n_points);
    for (int i = 0; i < n_points; ++i) {
        xs[i] = start + (end - start) * i / (n_points - 1);
    }

    const std::vector<double> expected = {
        1.14362132e-27, 5.17604145e-03, 7.65622444e-02, 1.97103480e-01, 3.17803769e-01, 4.22913972e-01, 5.10702567e-01,
        5.83253490e-01, 6.43255891e-01, 6.93126921e-01, 7.34842094e-01, 7.69965969e-01, 7.99727106e-01, 8.25091357e-01,
        8.46822292e-01, 8.65528152e-01, 8.81697515e-01, 8.95726201e-01, 9.07937577e-01, 9.18597931e-01, 9.27928197e-01,
        9.36112962e-01, 9.43307428e-01, 9.49642860e-01, 9.55230863e-01, 9.60166782e-01, 9.64532425e-01, 9.68398256e-01,
        9.71825180e-01, 9.74866008e-01, 9.77566657e-01, 9.79967157e-01, 9.82102482e-01, 9.84003253e-01, 9.85696330e-01,
        9.87205316e-01, 9.88550980e-01, 9.89751628e-01, 9.90823411e-01, 9.91780601e-01, 9.92635818e-01, 9.93400238e-01,
        9.94083761e-01, 9.94695172e-01, 9.95242266e-01, 9.95731968e-01, 9.96170434e-01, 9.96563142e-01, 9.96914964e-01,
        9.97230241e-01, 9.97512843e-01
    };

    const double rtol = 1e-8;

    for (int i = 0; i < n_points; ++i) {
        const double x = xs[i];
        const double ref = expected[i];
        const double result = xsf::cvm_cdf_inf(x);
        const auto rel_error = xsf::extended_relative_error(result, ref);

        CAPTURE(i, x, result, ref, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
