#include "../testing_utils.h"

#include <xsf/kelvin.h>

// Wronskian of the Kelvin equation, in the complex form the functions come in:
//
//     x * (Be*Kep - Bep*Ke) == -1
//
// with Be = ber + i*bei, Ke = ker + i*kei, and Bep, Kep their derivatives --
// exactly what xsf::kelvin() returns. Its real part is
// (ber*kerp - berp*ker) - (bei*keip - beip*kei) and its imaginary part is
// (ber*keip - berp*kei) + (bei*kerp - beip*ker), so one check covers all eight
// functions at once.

TEST_CASE("kelvin wronskian", "[ber][bei][ker][kei][kelvin][xsf_tests]") {
    // Loose enough to pass today. The worst residual measured on the current
    // implementation is 8.2e-11, at x = 9.9 -- just below the |x| < 10 branch
    // point, and about 111x worse than the first point above it.
    const double atol = 1e-9;

    const auto x = GENERATE(
        0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 9.5, 9.9, 10.0, 10.1, 10.5, 11.0, 15.0, 20.0, 30.0, 50.0
    );

    std::complex<double> Be, Ke, Bep, Kep;
    xsf::kelvin(x, Be, Ke, Bep, Kep);

    const auto wronskian = x * (Be * Kep - Bep * Ke);
    const auto residual = std::abs(wronskian + 1.0);
    CAPTURE(x, Be, Ke, Bep, Kep, residual, atol);
    REQUIRE(residual <= atol);
}
