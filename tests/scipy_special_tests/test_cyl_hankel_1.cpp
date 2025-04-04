#include "../testing_utils.h"

#include <xsf/bessel.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "cyl_hankel_1"};

TEST_CASE("cyl_hankel_1 dD->D scipy_special_tests", "[cyl_hankel_1][dD->D][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(
        xsf_test_cases<std::tuple<double, std::complex<double>>, std::tuple<std::complex<double>, bool>, double>(
            tables_path / "In_d_cd-cd.parquet", tables_path / "Out_d_cd-cd.parquet",
            tables_path / ("Err_d_cd-cd_" + get_platform_str() + ".parquet")
        )
    );

    auto [v, z] = input;
    auto [desired, fallback] = output;
    auto out = xsf::cyl_hankel_1(v, z);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(v, z, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
