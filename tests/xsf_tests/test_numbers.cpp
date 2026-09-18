#include "../testing_utils.h"
#include <tuple>
#include <xsf/numbers.h>

TEST_CASE("numbers.h", "[numbers][xsf_tests]") {
    REQUIRE(xsf::numbers::pi_v<float> == float(M_PI));
    REQUIRE(xsf::numbers::pi_v<double> == double(M_PI));

    // does not compile on all platforms, , is_floating_point<_Float16>::value is false
    if (sizeof(long double) > sizeof(double)) {
        REQUIRE(xsf::numbers::pi_v<long double> != double(M_PI));
    }

    // does not compile, is_floating_point<_Float16>::value is false
    // REQUIRE(xsf::cephes::pi_v<_Float16> != _Float16(M_PI));
}
