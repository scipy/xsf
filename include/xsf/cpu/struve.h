#pragma once

#include "../cephes/struve.h"

namespace xsf {
namespace detail {

    // Integrals of H0(t) (th0) and H0(t)/t (tth0) from 0 to x, for 0.5 <= x < 40. This integrates
    //   H0(t) = 4/pi sum_k J_{2k+1}(t) / (2k+1) (DLMF 11.4.21)
    // term by term, with
    //   int_0^x J_n(t) dt = 2 sum_{k >= 0} J_{n+2k+1}(x) (DLMF 10.22.8),
    //   int_0^x J_n(t) / t dt = (2 sum_{k >= 0} J_{n+2k}(x) - J_n(x)) / n (from DLMF 10.6.1),
    // and J_n(x) from Miller's backward recurrence normalized by DLMF 10.12.4.
    inline void itsh0_neumann(double x, double *th0, double *tth0) {
        const double pi = 3.141592653589793;
        double jp = 0.0, j = 1.0, jn;    // J_{n+1}(x) and J_n(x), unnormalized
        double so = 0.0, se = 0.0;       // sums of J_k(x) over odd and even k >= n
        double vo = 0.0, wo = 0.0;       // sums over odd k >= n of int_0^x J_k(t) dt / (2k) and int_0^x J_k(t)/t dt / k
        int n = 2 * (int)(0.5 * x) + 51; // J_n(x) < 1e-22
        for (; n >= 1; n--) {
            if (n % 2 == 1) {
                so += j;
                vo += se / n;
                wo += (2.0 * so - j) / (n * n);
            } else {
                se += j;
            }
            jn = 2.0 * n / x * j - jp;
            jp = j;
            j = jn;
        }
        const double norm = j + 2.0 * se; // J_0(x) + 2 J_2(x) + 2 J_4(x) + ... = 1
        *th0 = 8.0 / pi * vo / norm;
        *tth0 = 4.0 / pi * wo / norm;
    }

    inline double itth0(double x) {

        // ===========================================================
        // Purpose: Evaluate the integral H0(t)/t with respect to t
        //          from x to infinity
        // Input :  x   --- Lower limit  ( x ≥ 0 )
        // Output:  TTH --- Integration of H0(t)/t from x to infinity
        // ===========================================================

        int k;
        double f0, g0, r, s, t, tth, tty, xt;
        const double pi = 3.141592653589793;
        if (x >= 0.5 && x < 40.0) {
            double th0, tth0;
            itsh0_neumann(x, &th0, &tth0);
            return pi / 2.0 - tth0;
        }
        s = 1.0;
        r = 1.0;
        if (x < 24.5) {
            for (k = 1; k < 61; k++) {
                r = -r * x * x * (2.0 * k - 1.0) / pow(2.0 * k + 1.0, 3);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            tth = pi / 2.0 - 2.0 / pi * x * s;
        } else {
            for (k = 1; k < 11; k++) {
                r = -r * pow(2.0 * k - 1.0, 3) / ((2.0 * k + 1.0) * x * x);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            tth = 2.0 / (pi * x) * s;
            t = 8.0 / x;
            xt = x + 0.25 * pi;
            f0 = (((((0.18118e-2 * t - 0.91909e-2) * t + 0.017033) * t - 0.9394e-3) * t - 0.051445) * t - 0.11e-5) * t +
                 0.7978846;
            g0 =
                (((((-0.23731e-2 * t + 0.59842e-2) * t + 0.24437e-2) * t - 0.0233178) * t + 0.595e-4) * t + 0.1620695) *
                t;
            tty = (f0 * sin(xt) - g0 * cos(xt)) / (sqrt(x) * x);
            tth = tth + tty;
        }
        return tth;
    }

    inline double itsh0(double x) {

        // ===================================================
        // Purpose: Evaluate the integral of Struve function
        //          H0(t) with respect to t from 0 and x
        // Input :  x   --- Upper limit  ( x ≥ 0 )
        // Output:  TH0 --- Integration of H0(t) from 0 and x
        // ===================================================

        int k;
        double a[25], a0, a1, af, bf, bg, r, rd, s, s0, th0, ty, xp;
        const double pi = 3.141592653589793;
        const double el = 0.57721566490153;

        if (x >= 0.5 && x < 40.0) {
            double tth0;
            itsh0_neumann(x, &th0, &tth0);
            return th0;
        }
        r = 1.0;
        if (x <= 30.0) {
            s = 0.5;
            for (k = 1; k < 101; k++) {
                rd = 1.0;
                if (k == 1) {
                    rd = 0.5;
                }
                r = -r * rd * k / (k + 1.0) * pow(x / (2.0 * k + 1.0), 2);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            th0 = 2.0 / pi * x * x * s;
        } else {
            s = 1.0;
            for (k = 1; k < 13; k++) {
                r = -r * k / (k + 1.0) * pow((2.0 * k + 1.0) / x, 2);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            s0 = s / (pi * x * x) + 2.0 / pi * (log(2.0 * x) + el);
            a0 = 1.0;
            a1 = 5.0 / 8.0;
            a[0] = a1;
            for (k = 1; k < 21; k++) {
                af = ((1.5 * (k + 0.5) * (k + 5.0 / 6.0) * a1 - 0.5 * (k + 0.5) * (k + 0.5) * (k - 0.5) * a0)) /
                     (k + 1.0);
                a[k] = af;
                a0 = a1;
                a1 = af;
            }
            bf = 1.0;
            r = 1.0;
            for (k = 1; k < 11; k++) {
                r = -r / (x * x);
                bf += a[2 * k - 1] * r;
            }
            bg = a[0] / x;
            r = 1.0 / x;
            for (k = 1; k <= 10; k++) {
                r = -r / (x * x);
                bg += a[2 * k] * r;
            }
            xp = x + 0.25 * pi;
            ty = sqrt(2.0 / (pi * x)) * (bg * cos(xp) - bf * sin(xp));
            th0 = ty + s0;
        }
        return th0;
    }

    inline double itsl0(double x) {

        // ===========================================================
        // Purpose: Evaluate the integral of modified Struve function
        //          L0(t) with respect to t from 0 to x
        // Input :  x   --- Upper limit  ( x ≥ 0 )
        // Output:  TL0 --- Integration of L0(t) from 0 to x
        // ===========================================================

        // Number of terms of the asymptotic expansion for x > 40, where the last one is at most 2.2e-16
        constexpr int na = 25;
        int k;
        double a[na], a0, a1, af, r, rd, s, s0, ti, tl0;
        const double pi = 3.141592653589793;
        const double el = 0.57721566490153;
        r = 1.0;
        if (x <= 40.0) {
            s = 0.5;
            for (k = 1; k < 101; k++) {
                rd = 1.0;
                if (k == 1) {
                    rd = 0.5;
                }
                r = r * rd * k / (k + 1.0) * pow(x / (2.0 * k + 1.0), 2);
                s += r;
                if (fabs(r / s) < std::numeric_limits<double>::epsilon()) {
                    break;
                }
            }
            tl0 = 2.0 / pi * x * x * s;
        } else {
            s = 1.0;
            for (k = 1; k < 11; k++) {
                r = r * k / (k + 1.0) * pow((2.0 * k + 1.0) / x, 2);
                s += r;
                if (fabs(r / s) < 1.0e-12) {
                    break;
                }
            }
            s0 = -s / (pi * x * x) + 2.0 / pi * (log(2.0 * x) + el);
            a0 = 1.0;
            a1 = 5.0 / 8.0;
            a[0] = a1;
            for (k = 1; k < na; k++) {
                af = ((1.5 * (k + .50) * (k + 5.0 / 6.0) * a1 - 0.5 * pow(k + 0.5, 2) * (k - 0.5) * a0)) / (k + 1.0);
                a[k] = af;
                a0 = a1;
                a1 = af;
            }
            ti = 1.0;
            r = 1.0;
            for (k = 1; k <= na; k++) {
                r = r / x;
                ti += a[k - 1] * r;
            }
            tl0 = ti / sqrt(2 * pi * x) * exp(x) + s0;
        }
        return tl0;
    }

} // namespace detail

template <typename T>
T itstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = detail::itsh0(x);
    SPECFUN_CONVINF("itstruve0", out);
    return out;
}

template <typename T>
T it2struve0(T x) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }

    T out = detail::itth0(x);
    SPECFUN_CONVINF("it2struve0", out);
    if (flag) {
        out = M_PI - out;
    }
    return out;
}

template <typename T>
T itmodstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = detail::itsl0(x);
    SPECFUN_CONVINF("itmodstruve0", out);
    return out;
}

inline double struve_h(double v, double z) { return cephes::struve_h(v, z); }

inline float struve_h(float v, float z) { return struve_h(static_cast<double>(v), static_cast<double>(z)); }

inline double struve_l(double v, double z) { return cephes::struve_l(v, z); }

inline float struve_l(float v, float z) { return struve_l(static_cast<double>(v), static_cast<double>(z)); }

} // namespace xsf
