/* Translated from Cython into C++ by SciPy developers in 2023.
 *
 * Original author: Josh Wilson, 2016.
 */

#pragma once

#include "config.h"

namespace xsf {
namespace detail {

    XSF_HOST_DEVICE inline cxx::complex<double> zlog1(cxx::complex<double> z) {
        /* Compute log, paying special attention to accuracy around 1. We
         * implement this ourselves because some systems (most notably the
         * Travis CI machines) are weak in this regime. */
        cxx::complex<double> coeff = -1.0;
        cxx::complex<double> res = 0.0;

        if (cxx::abs(z - 1.0) > 0.1) {
            return cxx::log(z);
        }

        z -= 1.0;
        for (int n = 1; n < 17; n++) {
            coeff *= -z;
            res += coeff / static_cast<double>(n);
            if (cxx::abs(res / coeff) < cxx::numeric_limits<double>::epsilon()) {
                break;
            }
        }
        return res;
    }
} // namespace detail
} // namespace xsf
