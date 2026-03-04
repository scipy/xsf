#include "../testing_utils.h"
#include <xsf/stats.h>

namespace {

const int n_points = 51;
const double start = 2e-3;
const double end = 1.0 - 2e-3;

std::vector<double> make_xs() {
    std::vector<double> xs(n_points);
    for (int i = 0; i < n_points; ++i) {
        xs[i] = start + (end - start) * i / (n_points - 1);
    }
    return xs;
}

const std::vector<double> xs = make_xs();

void check_cdf(const std::vector<double> &expected, int n, double rtol = 1e-12) {
    for (int i = 0; i < n_points; ++i) {
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
    0.0030589048899929954,
    0.0673304065708938,
    0.18574849404342234,
    0.30720490368650677,
    0.413925988371887,
    0.5034242973475932,
    0.5775415452429734,
    0.6389054618066644,
    0.6899338489134077,
    0.7326228341757939,
    0.7685622230100523,
    0.7990048296623964,
    0.8249378564805038,
    0.8471428982610355,
    0.8662430863316606,
    0.8827391798704515,
    0.8970369467842458,
    0.9094679107388989,
    0.9203051143603324,
    0.9297751544920294,
    0.9380674282130427,
    0.9453412862104374,
    0.9517316094440662,
    0.9573531915590892,
    0.9623042112467177,
    0.9666690064825132,
    0.9705203093627187,
    0.9739210610063278,
    0.9769258969528326,
    0.9795823719303957,
    0.9819319768033516,
    0.9840109884777481,
    0.9858511844883349,
    0.9874804471351779,
    0.9889232768160711,
    0.9902012301965315,
    0.991333295767238,
    0.9923362169336412,
    0.9932247708982962,
    0.9940120101088166,
    0.994709471861316,
    0.9953273607013794,
    0.9958747074997525,
    0.9963595084585039,
    0.9967888467950952,
    0.9971689994334287,
    0.9975055306845169,
    0.997803374610999,
    0.9980669075283305,
    0.9983000118924547
};

} // namespace

TEST_CASE("cdf_cvm n=10", "[cdf_cvm][xsf_tests]") { check_cdf(expected_n_10, 10); }
TEST_CASE("cdf_cvm n=inf", "[cdf_cvm][xsf_tests]") { check_cdf(expected_n_inf, -1); }
TEST_CASE("cdf_cvm default parameter (n=-1)", "[cdf_cvm][xsf_tests]") {
    for (int i = 0; i < n_points; ++i) {
        const double res_default = xsf::cdf_cvm(xs[i]);
        const double res_explicit = xsf::cdf_cvm(xs[i], -1);
        CAPTURE(i, xs[i], res_default, res_explicit);
        REQUIRE(res_default == res_explicit);
    }
}
