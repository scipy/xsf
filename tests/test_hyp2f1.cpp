#include <catch2/catch_test_macros.hpp>

#include <xsf/hyp2f1.h>
#include <xsf/fp_error_metrics.h>

#include "testing_utils.h"


TEST_CASE("hyp2f1 example text", "[hyp2f1_example]") {
    REQUIRE ( xsf::extended_relative_error(xsf::hyp2f1(1.0, 0.9, 0.8, 0.2), 1.2848765105679898) < 1e-12 );
}
