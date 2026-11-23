#include "../testing_utils.h"

#include <limits>

#include <xsf/cephes/pdtr.h>

TEST_CASE("cephes::pdtri y = 1", "[cephes::pdtri][xsf_tests]") {
    // pdtr(k, 0) = 1 exactly, so pdtri(k, 1) = 0. It used to return nan.
    for (int k : {0, 1, 10, std::numeric_limits<int>::max()}) {
        CAPTURE(k);
        REQUIRE(xsf::cephes::pdtri(k, 1.0) == 0.0);
    }
}
