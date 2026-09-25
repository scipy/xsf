#include "../testing_utils.h"
#include <tuple>
#include <xsf/numbers.h>

TEST_CASE("numbers.h", "[numbers][xsf_tests]") {
    REQUIRE(xsf::numbers::pi_v<float> == float(M_PI));
    REQUIRE(xsf::numbers::pi_v<double> == double(M_PI));
#if __STDCPP_FLOAT16_T__
    REQUIRE(xsf::numbers::pi_v<_Float16> == _Float16(M_PI));
#endif

    // Check that long double pi_v is more precise than double precision,
    // but only if the two types have different sizes.
    if constexpr (sizeof(long double) > sizeof(double)) {
        REQUIRE(xsf::numbers::pi_v<long double> != double(M_PI));
    }
}
