#include "../testing_utils.h"
#include <tuple>
#include <xsf/hyp2f1.h>

TEST_CASE("hyp2f1(-1, b, c, z)", "[hyp2f1][xsf_tests]") {
    constexpr double rtol = 1e-14;
    constexpr double atol = 1e-14; // absolute tolerance

    constexpr int size = 201;
    std::vector<double> b_values(size), c_values(size), z_values(size);

    for (int i = 0; i < size; ++i) {
        b_values[i] = -10.0 + i * 0.1;
        c_values[i] = -10.0 + i * 0.1;
        z_values[i] = -10.0 + i * 0.1;
    }

    for (double b : b_values) {
        for (double c : c_values) {
            if (c == 0.0)
                continue;

            for (double z : z_values) {
                double expected = 1 - (b / c) * z;
                double result = xsf::hyp2f1(-1.0, b, c, z);

                double error;
                if (std::abs(expected) < atol) {
                    error = std::abs(result - expected);
                } else {
                    error = xsf::extended_relative_error(result, expected);
                }

                CAPTURE(b, c, z, result, expected, rtol, atol, error);
                REQUIRE(error <= std::max(rtol * std::abs(expected), atol));
            }
        }
    }
}
