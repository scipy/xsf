#include "../testing_utils.h"

#include <xsf/bessel.h>

namespace fs = std::filesystem;

fs::path tables_path{fs::path(XSREF_TABLES_PATH) / "scipy_special_tests" / "cyl_bessel_i1e"};

TEST_CASE("cyl_bessel_i1e f->f scipy_special_tests", "[cyl_bessel_i1e][f->f][scipy_special_tests]") {
    SET_FP_FORMAT()
    auto [input, output, tol] = GENERATE(xsf_test_cases<float, std::tuple<float, bool>, float>(
        tables_path / "In_f-f.parquet", tables_path / "Out_f-f.parquet",
        tables_path / ("Err_f-f_" + get_platform_str() + ".parquet")
    ));

    auto z = input;
    auto [desired, fallback] = output;
    auto out = xsf::cyl_bessel_i1e(z);
    auto error = xsf::extended_relative_error(out, desired);
    tol = adjust_tolerance(tol);
    CAPTURE(z, out, desired, error, tol, fallback);
    REQUIRE(error <= tol);
}
