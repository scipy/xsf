#pragma once

#include "cephes/polevl.h"
#include "config.h"
#include "evalpoly.h"
#include "trig.h"

#include <complex>
#include <type_traits>

namespace xsf {

// https://github.com/JuliaMath/SpecialFunctions.jl/blob/c9ff36085c9008b9bdaa457e745f459a1827e52e/src/gamma.jl#L65
template <typename T>
XSF_HOST_DEVICE inline T trigamma(T z) {
    double x = std::real(z);
    if (x <= 0) { // reflection formula
        T s = xsf::sinpi(z);
        T c = M_PI / s;
        return c * c - trigamma(T(1) - z);
    }

    T psi = 0;
    int N = 10;
    int n;
    if (x < N) {
        n = N - std::floor(x);
        psi += 1.0 / (z * z);
        for (int v = 1; v < n; v++) {
            T w = z + static_cast<T>(v);
            psi += 1.0 / (w * w);
        }
        z += n;
    }
    T t = 1.0 / z;
    T w = t * t;
    psi += t + 0.5 * w;
    double coeffs[] = {-7.092156862745098,   1.1666666666666667,   -0.2531135531135531,  0.07575757575757576,
                       -0.03333333333333333, 0.023809523809523808, -0.03333333333333333, 0.16666666666666666};
    using base_t = std::remove_cv_t<T>;
    if constexpr (std::is_same_v<base_t, std::complex<double>> || std::is_same_v<base_t, std::complex<float>>) {
        psi += t * w * xsf::cevalpoly(coeffs, 7, w);
    } else {
        psi += t * w * xsf::cephes::polevl(w, coeffs, 7);
    }
    return psi;
}

} // namespace xsf
