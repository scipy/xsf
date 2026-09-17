#include "../testing_utils.h"
#include <xsf/stats.h>

TEST_CASE("Tukey lambda cdf goes to 1 for large x", "[tukey][xsf_tests]") {
    // Check that survival function goes to 0 as x gets large

    std::vector<double> xs = {1e4, 1e6, 1e8, 1e10, 1e12, 1e24};
    double lambda = -0.5;

    // Check that cdf is 1 for largest x
    double cdf = xsf::tukeylambdacdf(xs.back(), lambda);
    REQUIRE(cdf == 1.0);

    double cdf_prev;
    for (unsigned i = 1; i < xs.size(); i++) {
        // Check that the cdf is always increasing
        cdf = xsf::tukeylambdacdf(xs[i], lambda);
        cdf_prev = xsf::tukeylambdacdf(xs[i - 1], lambda);
        REQUIRE(cdf - cdf_prev >= 0);
    }
}
