#include "../testing_utils.h"
#include <xsf/stats.h>

namespace {

std::vector<double> make_xs(double start, double end, int n_points) {
    std::vector<double> xs(n_points);
    for (int i = 0; i < n_points; ++i) {
        xs[i] = start + (end - start) * i / (n_points - 1);
    }
    return xs;
}

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
//
// np.set_printoptions(precision=20)
// xs = np.linspace(2e-3, 1 - 2e-3, 51)
// expected_n_inf = _cdf_cvm(xs)
// expected_n_10 = _cdf_cvm(xs, n=10)
// assert np.allclose(expected_n_inf, _cdf_cvm_inf(xs))
// print("expected_n_inf:", expected_n_inf)
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
    0.30142064363602883,
    0.6171562716813203,
    0.7781501473708698,
    0.8657378548074194,
    0.9166296801414354,
    0.947471730015775,
    0.9666408872017962,
    0.9787337815065694,
    0.9864289403414447,
    0.9913489626856629,
    0.9945012589632174,
    0.9965211845209974,
    0.9978135194658495,
    0.9986378194232566,
    0.999161174297481,
    0.999491369374648,
    0.9996979717974336,
    0.9998258481904505,
    0.9999038787884099,
    0.9999505957510425,
    0.9999778393308453,
    0.9999931309483924,
    1.000001213309659,
    1.0000050487947476,
    1.0000064653628036,
    1.0000065731540428,
    1.000006032108419,
    1.0000052229847072,
    1.000004355929261,
    1.0000035388284674,
    1.0000028198915138,
    1.000002213823316,
    1.000001717629139,
    1.0000013199320752,
    1.0000010062820168,
    1.0000007620926568,
    1.000000573851733,
    1.0000004299797818,
    1.0000003207931465,
    1.0000002384261597,
    1.000000176611692,
    1.0000001304297244,
    1.0000000960629287,
    1.0000000705776748,
    1.000000051737207,
    1.,
    1.,
    1.,
    1.,
    1.
};

const int n_points = 51;
const std::vector<double> xs_n_inf = make_xs(2e-3, 1.0 - 2e-3, n_points);
// For n=10, the finite-sample support is [1/(12*n), n/3] = [1/120, 10/3].
// We sample 10% beyond both bounds to include points where the CDF is clamped to 0 and 1.
// For large values of x, expected CDF values exceed 1.0 due to the limitations of the approximation
// by Csörgő and Faraway (1996) implemented in xsf::cdf_cvm
const int n = 10;
const std::vector<double> xs_n_10 = make_xs((1.0 / (12.0 * n)) * (1.0 - 1e-1), (n / 3.0) * (1.0 + 1e-1), n_points);

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
