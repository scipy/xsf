#include "../testing_utils.h"

#include <xsf/cpu/specfun/specfun.h>
#include <xsf/zeta.h>

TEST_CASE("specfun::bernob", "[bernob][xsf_tests]") {
    // Check B_0, ..., B_269 against B_2k = (-1)**(k+1) * 2 * (2k)! * zeta(2k) / (2 pi)**(2k) (DLMF 25.6.2).
    // The zeta(n) sum for B_n used to be cut off too early (relative error 1.7e-12 for B_4).
    SET_FP_FORMAT()
    double bn[270] = {};
    xsf::specfun::bernob(269, bn);

    REQUIRE(bn[0] == 1.0);
    REQUIRE(bn[1] == -0.5);
    REQUIRE(bn[2] == 1.0 / 6.0);
    double f = 2.0 * 2.0 / (2.0 * M_PI) / (2.0 * M_PI); // 2 * (2k)! / (2 pi)**(2k) for k = 1
    for (int n = 3; n < 270; n++) {
        if (n % 2 == 1) {
            CAPTURE(n, bn[n]);
            REQUIRE(bn[n] == 0.0);
            continue;
        }
        f *= (n - 1) / (2.0 * M_PI) * n / (2.0 * M_PI);
        const double desired = ((n % 4 == 2) ? 1.0 : -1.0) * f * xsf::zeta(n, 1.0);
        if (n >= 260) {
            // |B_n| overflows
            CAPTURE(n, bn[n], desired);
            REQUIRE(bn[n] == desired);
            continue;
        }
        const double error = xsf::extended_relative_error(bn[n], desired);
        const double rtol = 5e-14;
        CAPTURE(n, bn[n], desired, error, rtol);
        REQUIRE(error <= rtol);
    }
}
