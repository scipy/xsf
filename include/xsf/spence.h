/*
Computes the Spence function (alternative definition of the dilogarithm) for complex z.
Uses the fast converging series [1] (5) near z = 1, and uses argument transformation
formulae [1] (3), [1] (4) elsewhere. The series is evaluated using the Estrin scheme [2].

Sources:
[1] Vollinga, Weinzierl, "Numerical evaluation of multiple polylogarithms." https://arxiv.org/abs/hep-ph/0410259
[2] https://en.wikipedia.org/wiki/Estrin%27s_scheme
[3] https://en.wikipedia.org/wiki/Sterbenz_lemma

Author: Rudolf Rosendorf
Date: 20/10/2026 (dd/mm/yyyy)

*/

#pragma once

#include "cephes/spence.h"
#include "config.h"
#include "zlog1.h"


namespace xsf {
namespace detail {

    constexpr double PI_SQ_6 = 1.64493406684822644;  // π²/6

    // The algorithm splits complex numbers to their real and imag parts for work to execute FMA instructions on them.
    // FMA(a, b, c) = a * b + c
    // This should results in faster speeds and better errors.


    XSF_HOST_DEVICE inline void series(
        double * __restrict resr,
        double * __restrict resi,
        const double ur, 
        const double ui) 
        {

        // Computes the series S(u).
        // S(-log(1 - z)) ≈ Li2(z)
        // S(-log(z)) ≈ Li2(1 - z) = Spence(z)
        
        // v = u²
        double vr = cxx::fma(ur, ur, -ui * ui);
        double vi = ur * ui, t;
        vi += vi;

        const double wr = vr, wi = vi;  // save w = v for later use.
        
        // Compute the truncated series s [1] (5) using the Estrin scheme [2].
        // res = u - 1/4 u² + 1/36 u³ - 1/3600 u⁵ + 1/211680 u⁷ - ...
        // s = 1/36 - 1/3600 u² + 1/211680 u⁴ - ...
        // res = u + u²(-1/4 + u s)

        const double p01r = cxx::fma(-0.000277777777777777778, vr, 0.0277777777777777778);
        const double p01i = -0.000277777777777777778 * vi;

        const double p23r = cxx::fma(-9.18577307466196355e-8, vr, 4.72411186696900983e-6);
        const double p23i = -9.18577307466196355e-8 * vi;

        const double p45r = cxx::fma(-4.06476164514422553e-11, vr, 1.89788699889709991e-9);
        const double p45i = -4.06476164514422553e-11 * vi;

        const double p67r = cxx::fma(-1.99392958607210757e-14, vr, 8.92169102045645256e-13);
        const double p67i = -1.99392958607210757e-14 * vi;
        
        const double p89r = cxx::fma(-1.03565176121812470e-17, vr, 4.51898002961991819e-16);
        const double p89i = -1.03565176121812470e-17 * vi;
        
        // v = v²  (u⁴)
        t = cxx::fma(vr, vr, -vi * vi);
        vi = vr * vi;
        vi += vi;
        vr = t;
        
        // q03 = p01 + p23 v
        const double q03r = p01r + cxx::fma(p23r, vr, -p23i * vi);
        const double q03i = p01i + cxx::fma(p23r, vi, p23i * vr);

        // q47 = p45 + p67 v
        const double q47r = p45r + cxx::fma(p67r, vr, -p67i * vi);
        const double q47i = p45i + cxx::fma(p67r, vi, p67i * vr);
        
        // v = v²  (u⁸)
        t = cxx::fma(vr, vr, -vi * vi);
        vi = vr * vi;
        vi += vi;
        vr = t;
        
        // r07 = q03 + q47 v
        const double r07r = q03r + cxx::fma(q47r, vr, -q47i * vi);
        const double r07i = q03i + cxx::fma(q47r, vi, q47i * vr);
        
        // v = v²  (u¹⁶)
        t = cxx::fma(vr, vr, -vi * vi);
        vi = vr * vi;
        vi += vi;
        vr = t;
        
        // s = r07 + p89 v
        double sr = r07r + cxx::fma(p89r, vr, -p89i * vi);
        double si = r07i + cxx::fma(p89r, vi, p89i * vr);
        
        // s = u s - 1/4;
        t = cxx::fma(sr, ur, -si * ui) - 0.25;
        si = cxx::fma(sr, ui, si * ur);
        sr = t;

        // res = u + u² s
        *resr = ur + cxx::fma(sr, wr, -si * wi);
        *resi = ui + cxx::fma(sr, wi, si * wr);
    }


    XSF_HOST_DEVICE inline cxx::complex<double> spence_algorithm(const cxx::complex<double> z) {
        
        const double x = z.real(), y = z.imag();
        const double y2 = y * y;

        // Obtain |z|² and |1 - z|² - 1:
        // |z|² = x² + y²
        // |1 - z|² = (1 - x)² + (-y)² = 1 - 2x + x² + y²
        // |1 - z|² - 1 = |z|² - x - x
        const double v = cxx::fma(x, x, y2); // |z|²
        const double w = v - x - x; // |1 - z|² - 1
        double sr, si;

        // Domain is split into three regions:
        // Region 1: Re(z) > 1/2  &&  |1 - z| < 1  :  No transformation.
        // Region 2: Re(z) < 1/2  &&  |z| < 1      :  Li₂(1 - z) = π²/6 - log(z)log(1 - z) - Li₂(z)
        // Region 3: Elsewhere                     :  Li₂(1 - z) = -π²/6 - 1/2 log²(z - 1) - Li₂(1 / (1 - z))
        // The transformation formulae were derived from [1] (3) and [1] (4).
        
        if (x < 0.5) {
            if (v < 1.0) {  // Region 2

                if (v == 0.0) [[unlikely]] {
                    if (y == 0.0) {
                        return cxx::complex<double>{PI_SQ_6, 0.0};
                    }
                    return cxx::complex<double>{PI_SQ_6, y * (cxx::log(cxx::abs(y)) - 1.0)};
                }
                
                // -log(1 - z) = -1/2 log(|1 - z|²) - i Arg(1 - z) = -1/2 log1p(w) - i atan2(-y, 1 - x)
                // = -1/2 log1p(w) + i atan2(y, 1 - x)
                double lr = -0.5 * cxx::log1p(w), li = cxx::atan2(y, 1.0 - x);
                // log(z) = 1/2 log(|z|²) + i Arg(z) = 1/2 log(v) + i atan2(y, x)
                double mr = 0.5 * cxx::log(v), mi = cxx::atan2(y, x);

                series(&sr, &si, lr, li);  // S(-log(1 - z)) ≈ Li2(z)
                
                return cxx::complex<double>{PI_SQ_6 + cxx::fma(lr, mr, -li * mi) - sr, 
                                            cxx::fma(lr, mi, li * mr) - si};
            }
        } else {
            if (w < 0.0) {  // Region 1

                // To obtain accurate -log(z) here, we write:
                // u = 1 - x
                // In this region 1/2 < x < 2, which perfectly satisfies the Sterbenz Lemma [3]. u is accurate.
                // |z|² - 1 = u² + y² - 2u
                // Then we can evaluate: -log(z) = -1/2 log1p(u² + y² - 2u) - i atan2(y, x);

                const double u = 1.0 - x;
                series(&sr, &si, -0.5 * cxx::log1p(cxx::fma(u, u, y2 - u - u)), -cxx::atan2(y, x));
                return cxx::complex<double>{sr, si};
            }   
        }

        // Region 3

        if (cxx::isinf(v)) [[unlikely]] {
            // Li₂(1 - z) = -π²/6 - 1/2 log²(z - 1) - Li₂(1 / (1 - z)) ~ -π²/6 - 1/2 log²(z - 1) ~ -π²/6 - 1/2 log²(z)
            cxx::complex<double> l = cxx::log(z);
            return -PI_SQ_6 - 0.5 * l * l;
        }

        // log(z) = 1/2 log(|z|²) + i Arg(z) = 1/2 log(v) + i atan2(y, x)
        double mr = 0.5 * cxx::log(v), mi = cxx::atan2(y, x);
        // log(z - 1) = 1/2 log(|1 - z|²) + i Arg(z - 1) = 1/2 log1p(w) + i atan2(y, x - 1)
        double nr = 0.5 * cxx::log1p(w), ni = cxx::atan2(y, x - 1.0);
        
        // A trick to get rid of the division in this region.
        // Li₂(1 / (1 - z)) ≈ S(-log(1 - 1 / (1 - z))) = S(-log(z / (z - 1))) = S(log(z - 1) - log(z))
        series(&sr, &si, nr - mr, ni - mi);

        return -cxx::complex<double>{PI_SQ_6 + 0.5 * cxx::fma(nr, nr, -ni * ni) + sr,
                                        nr * ni + si};
    }

} // namespace detail

XSF_HOST_DEVICE inline cxx::complex<double> spence(cxx::complex<double> z) {
 
    if (cxx::isnan(z.real()) || cxx::isnan(z.imag())) [[unlikely]] {
        return cxx::complex<double>{cxx::numeric_limits<double>::quiet_NaN(), cxx::numeric_limits<double>::quiet_NaN()};
    } 
    return detail::spence_algorithm(z);
}

XSF_HOST_DEVICE inline cxx::complex<float> spence(cxx::complex<float> z) {
    return static_cast<cxx::complex<float>>(spence(static_cast<cxx::complex<double>>(z)));
}

XSF_HOST_DEVICE inline double spence(double z) {
    if (cxx::isnan(z)) {
        return cxx::numeric_limits<double>::quiet_NaN();
    }
    return xsf::cephes::spence(z);
}

XSF_HOST_DEVICE inline float spence(float z) { return static_cast<float>(spence(static_cast<double>(z))); }

} // namespace xsf
