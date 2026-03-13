#include "../testing_utils.h"
#include <xsf/stats.h>

namespace {

void check_cdf(const std::vector<double> &xs, const std::vector<double> &expected, int n, double rtol = 1e-12) {
    for (int i = 0; i < xs.size(); ++i) {
        const double result = xsf::cdf_cvm(xs[i], n);
        const auto rel_error = xsf::extended_relative_error(result, expected[i]);
        CAPTURE(i, xs[i], n, result, expected[i], rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

// Reference values computed with scipy.stats._hypotests._cdf_cvm_inf
// import numpy as np
// from scipy.stats._hypotests import _cdf_cvm, _cdf_cvm_inf

// np.set_printoptions(precision=20)

// print("n=inf")
// xs_n_inf = np.linspace(2e-3, 1 - 2e-3, 51)
// expected_n_inf = _cdf_cvm(xs_n_inf)
// assert np.allclose(expected_n_inf, _cdf_cvm_inf(xs_n_inf))
// print("expected_n_inf:", expected_n_inf)

// print("\nn=10")
// n = 10  # support is [1/(12*n), n/3]
// lower = 1 / (12 * n)
// upper = n / 3
// eps_abs = 0.1 * (upper - lower)
// xs_n_10 = np.linspace(lower - eps_abs, upper + eps_abs, 51)
// expected_n_10 = np.clip(_cdf_cvm(xs_n_10, n=10), 0, 1)
// print("expected_n_10:", expected_n_10)

const std::vector<double> expected_n_inf = {
    1.1436213247613867e-27, 5.1760414467472869e-03, 7.6562244376509977e-02, 1.9710348009439807e-01,
    3.1780376868441079e-01, 4.2291397179060869e-01, 5.1070256741610731e-01, 5.8325348963163892e-01,
    6.4325589114353210e-01, 6.9312692117440566e-01, 7.3484209398823697e-01, 7.6996596903666392e-01,
    7.9972710572561168e-01, 8.2509135670539901e-01, 8.4682229224958638e-01, 8.6552815214794321e-01,
    8.8169751461339307e-01, 8.9572620130472891e-01, 9.0793757738532332e-01, 9.1859793103629028e-01,
    9.2792819736657439e-01, 9.3611296152331869e-01, 9.4330742791372113e-01, 9.4964286011501597e-01,
    9.5523086307243110e-01, 9.6016678243439135e-01, 9.6453242541194095e-01, 9.6839825608258112e-01,
    9.7182518031147813e-01, 9.7486600764169828e-01, 9.7756665688682531e-01, 9.7996715678835211e-01,
    9.8210248156904612e-01, 9.8400325250598897e-01, 9.8569633002701529e-01, 9.8720531576628501e-01,
    9.8855098010733300e-01, 9.8975162770810998e-01, 9.9082341112961814e-01, 9.9178060082134489e-01,
    9.9263581823463787e-01, 9.9340023765193020e-01, 9.9408376136868304e-01, 9.9469517209560188e-01,
    9.9524226582262543e-01, 9.9573196787365759e-01, 9.9617043445916942e-01, 9.9656314168469795e-01,
    9.9691496368300525e-01, 9.9723024129513893e-01, 9.9751284252216244e-01
};

const std::vector<double> expected_n_10 = {
    0.,
    0.,
    0.,
    0.,
    0.,
    0.26671805597982823,
    0.6194519411038519,
    0.7894518664504734,
    0.87762684010206,
    0.926884872096052,
    0.9556668307472429,
    0.9729196588715867,
    0.9834106045043685,
    0.9898395691044443,
    0.9937938275270536,
    0.996228007757172,
    0.9977243562086301,
    0.9986410813136646,
    0.9991996874665492,
    0.9995374987804108,
    0.9997397024352754,
    0.9998590871289252,
    0.9999282800630996,
    0.9999673639639962,
    0.9999886286430856,
    0.9999995369701493,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.,
    1.
};

const int n_points = 51;
// Case: n=inf
const std::vector<double> xs_n_inf = linspace(2e-3, 1.0 - 2e-3, n_points);
// Case: n=10
// For n=10, the finite-sample support is [1/(12*n), n/3] = [1/120, 10/3].
// This interval is extended symmetrically by +/- eps_abs, where
// eps_abs = 10% of the support width (upper - lower), so the test
// includes points outside the support where cdf_cvm should be 0 or 1.
const int n = 10;
const double lower = 1.0 / (12.0 * n);
const double upper = n / 3.0;
const double eps_abs = 0.1 * (upper - lower);
const std::vector<double> xs_n_10 = linspace(lower - eps_abs, upper + eps_abs, n_points);

} // namespace

TEST_CASE("cdf_cvm n=10", "[cdf_cvm][xsf_tests]") { check_cdf(xs_n_10, expected_n_10, 10); }
TEST_CASE("cdf_cvm n=inf", "[cdf_cvm][xsf_tests]") { check_cdf(xs_n_inf, expected_n_inf, -1); }
TEST_CASE("cdf_cvm default parameter (n=-1)", "[cdf_cvm][xsf_tests]") {
    for (int i = 0; i < n_points; ++i) {
        const double res_default = xsf::cdf_cvm(xs_n_inf[i]);
        const double res_explicit = xsf::cdf_cvm(xs_n_inf[i], -1);
        CAPTURE(i, xs_n_inf[i], res_default, res_explicit);
        REQUIRE(res_default == res_explicit);
    }
}
