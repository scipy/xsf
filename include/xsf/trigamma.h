/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2017 Jeff Bezanson, Stefan Karpinski, Viral B. Shah, and others:
 *
 * https://github.com/JuliaMath/SpecialFunctions.jl/graphs/contributors
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

// https://github.com/JuliaMath/SpecialFunctions.jl/blob/c9ff36085c9008b9bdaa457e745f459a1827e52e/src/gamma.jl#L65
// Implementation based on differentiation of equation 12 of
// K.S. Kölbig, Programs for computing the logarithm of the gamma function, and the
// digamma function, for complex argument, Computer Physics Communications, 1972
// https://doi.org/10.1016/0010-4655(72)90012-4.
template <typename T>
XSF_HOST_DEVICE inline T trigamma(T z) {
    double x = std::real(z);
    if (x <= 0) {
        // reflection formula http://dlmf.nist.gov/5.15.E6
        T s = sinpi(z);
        T c = M_PI / s;
        return c * c - trigamma(T(1) - z);
    }

    T psi = 0;
    int N = 10;
    int n;
    if (x < N) {
        // recursion based on https://dlmf.nist.gov/5.15.E1
        n = N - std::floor(x);
        psi += 1.0 / (z * z);
        for (int v = 1; v < n; v++) {
            T w = z + static_cast<T>(v);
            psi += 1.0 / (w * w);
        }
        z += n;
    }

    // https://dlmf.nist.gov/5.15.E8
    T t = 1.0 / z;
    T w = t * t;
    psi += t + 0.5 * w;
    // order of the Bernoulli numbers is reversed due to `evalpoly` convention
    double coeffs[] = {-7.092156862745098,   1.1666666666666667,   -0.2531135531135531,  0.07575757575757576,
                       -0.03333333333333333, 0.023809523809523808, -0.03333333333333333, 0.16666666666666666};
    using base_t = std::remove_cv_t<T>;
    psi += t * w * evalpoly(coeffs, 7, w);
    return psi;
}

template <>
XSF_HOST_DEVICE inline float trigamma(float z) {
    return static_cast<float>(trigamma(static_cast<double>(z)));
}

template <>
XSF_HOST_DEVICE inline std::complex<float> trigamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(trigamma(static_cast<std::complex<double>>(z)));
}

} // namespace xsf
