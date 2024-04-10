#pragma once

#include "error.h"
#include "legendre.h"
#include "mdspan.h"

#include "specfun.h"

#include "cephes/poch.h"

namespace special {

template <typename T>
std::complex<T> sph_harm(long m, long n, T theta, T phi) {
    if (n < 0) {
        set_error("sph_harm", SF_ERROR_ARG, "n should not be negative");
        return NAN;
    }

    if (m < 0) {
        mp = -m;
        prefactor = std::pow(-1, mp) * cephes::poch(n + mp + 1, -2 * mp);
    } else {
        mp = m;
    }

    std::complex<T> val = pmv(m_abs, n, std::cos(phi));
    if (m < 0) {
        val *= std::pow(-1, m_abs) * cephes::poch(n + m_abs + 1, -2 * m_abs);
    }

    val *= std::sqrt((2 * n + 1) / 4.0 / M_PI);
    val *= std::sqrt(cephes::poch(n + m + 1, -2 * m));
    val *= std::exp(std::complex<double>(0, m * theta));

    return val;
}

template <typename T, typename OutMat>
void sph_harm_all(T theta, T phi, OutMat y) {
    long m = (y.extent(0) - 1) / 2;
    long n = y.extent(1) - 1;

    OutMat y_pos = std::submdspan(y, std::make_tuple(0, m + 1), std::full_extent);
    sph_legendre_all(phi, y_pos);

    for (long j = 0; j <= n; ++j) {
        for (long i = 1; i <= j; ++i) {
            y(i, j) *= static_cast<T>(std::sqrt((2 * j + 1) * cephes::poch(j + i + 1, -2 * i) / (4 * M_PI))) *
                       std::exp(std::complex<T>(0, i * theta));
            y(y.extent(0) - i, j) = static_cast<T>(std::pow(-1, i)) * std::conj(y(i, j));
        }
    }
}

} // namespace special
