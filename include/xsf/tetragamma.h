/*
 * The MIT License (MIT)
 *
 * Copyright 2026 Steven T. Smith
 *
 * https://github.com/scipy/scipy/pull/17933#issuecomment-5235227979
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include "cephes/polevl.h"
#include "config.h"
#include "evalpoly.h"
#include "trig.h"

#include <complex>
#include <type_traits>

namespace xsf {

// Asymptotic coefficients d_k = B_{2k}*(2k+1), k=1..10 (A&S 6.4.11, m=2)
template <typename T>
XSF_HOST_DEVICE inline T tetragamma(T z) {
    using base_t = std::remove_cv_t<T>;
    constexpr bool is_complex_t =
        std::is_same_v<base_t, std::complex<double>> || std::is_same_v<base_t, std::complex<float>>;

    if constexpr (is_complex_t) {
        if (std::imag(z) == 0) {
            # also catches 0.+0.j → -Inf because M_PI / 0. is Inf in IEEE 754,
            # but M_PI+0j / 0+0j is Nan+Nan*1j in complex
            using real_t = typename base_t::value_type;
            return T(tetragamma(static_cast<real_t>(std::real(z))), 0);
        }
    }

    double x = std::real(z);
    if (x <= 0) {
        // omit numerically 0 term and avoid overflow at large imag(z)
        bool use_trig_term = true;
        if constexpr (is_complex_t) {
            using real_t = typename base_t::value_type;
            use_trig_term = std::abs(std::imag(z)) < std::log(std::numeric_limits<real_t>::max());
        }
        if (use_trig_term) {
            T s = sinpi(z);
            T c = M_PI / s;
            return tetragamma(T(1) - z) - T(2) * cospi(z) * c * c * c;
        } else {
            return tetragamma(T(1) - z);
        }
    }

    T psi = 0;
    int N = 10;
    int n;
    if (x < N) {
        n = N - std::floor(x);
        psi -= 2.0 / (z * z * z);
        for (int v = 1; v < n; v++) {
            T w = z + static_cast<T>(v);
            psi -= 2.0 / (w * w * w);
        }
        z += n;
    }
    T t = 1.0 / z;
    T w = t * t;
    psi += -w * (1.0 + t);
    double coeffs[] = {-11111.609090909091, 1044.452380952381, -120.56666666666666,
                        17.5, -3.2904761904761903, 0.83333333333333337,
                        -0.29999999999999999, 0.16666666666666666,
                        -0.16666666666666666, 0.5};
    if constexpr (is_complex_t) {
        psi -= w * w * cevalpoly(coeffs, 9, w);
    } else {
        psi -= w * w * cephes::polevl(w, coeffs, 9);
    }
    return psi;
}

template <>
XSF_HOST_DEVICE inline float tetragamma(float z) {
    return static_cast<float>(tetragamma(static_cast<double>(z)));
}

template <>
XSF_HOST_DEVICE inline std::complex<float> tetragamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(tetragamma(static_cast<std::complex<double>>(z)));
}

} // namespace xsf
