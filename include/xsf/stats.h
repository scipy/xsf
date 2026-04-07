#pragma once

#include "xsf/bessel.h"
#include "xsf/binom.h"
#include "xsf/cephes/bdtr.h"
#include "xsf/cephes/chdtr.h"
#include "xsf/cephes/const.h"
#include "xsf/cephes/fdtr.h"
#include "xsf/cephes/gdtr.h"
#include "xsf/cephes/incbet.h"
#include "xsf/cephes/incbi.h"
#include "xsf/cephes/kolmogorov.h"
#include "xsf/cephes/nbdtr.h"
#include "xsf/cephes/ndtr.h"
#include "xsf/cephes/ndtri.h"
#include "xsf/cephes/owens_t.h"
#include "xsf/cephes/pdtr.h"
#include "xsf/cephes/tukey.h"
#include "xsf/erf.h"
#include "xsf/gamma.h"

namespace xsf {

inline double bdtr(double k, int n, double p) { return cephes::bdtr(k, n, p); }

inline double bdtri(double k, int n, double y) { return cephes::bdtri(k, n, y); }

inline double bdtrc(double k, int n, double p) { return cephes::bdtrc(k, n, p); }

inline double chdtr(double df, double x) { return cephes::chdtr(df, x); }

inline float chdtr(float df, float x) { return static_cast<float>(cephes::chdtr(df, x)); }

inline double chdtrc(double df, double x) { return cephes::chdtrc(df, x); }

inline float chdtrc(float df, float x) { return static_cast<float>(cephes::chdtrc(df, x)); }

inline double chdtri(double df, double y) { return cephes::chdtri(df, y); }

inline double cdf_cvm_inf(double x) {
    // CDF of the Cramer-von Mises test statistic (infinite sample limit).
    // Accurate for practical hypothesis testing but not expected to be accurate
    // for large values of x, e.g. x > 4, when the cdf is very close to 1.
    if (std::isnan(x)) {
        return x;
    }
    if (x <= 0) {
        return 0.0;
    }

    double total = 0.0;
    int k = 0;

    while (true) {
        double u = std::exp(gammaln(k + 0.5) - gammaln(k + 1.0)) / (std::pow(M_PI, 1.5) * std::sqrt(x));
        double y = 4.0 * k + 1.0;
        double q = (y * y) / (16.0 * x);
        double b = cyl_bessel_k(0.25, q);
        double z = u * std::sqrt(y) * std::exp(-q) * b;

        total += z;
        if (std::abs(z) < 1e-7) {
            break;
        }
        ++k;
    }

    return total;
}

// Helpers for psi1 computations
inline double psi1_G(double y) {
    // Helper implementing the G-type Bessel combination appearing in equation 1.10
    // in Csörgő, S. and Faraway, J. (1996).
    double z = y * y / 4.0;
    double b = cyl_bessel_k(0.25, z) + cyl_bessel_k(0.75, z);
    return std::exp(-z) * std::pow(y / 2.0, 1.5) * b / std::sqrt(M_PI);
}

inline double psi1_H(double y) {
    // Helper implementing the H-type Bessel combination appearing in equation 1.10
    // in Csörgő, S. and Faraway, J. (1996).
    double z = y * y / 4.0;
    double b = 2.0 * cyl_bessel_k(0.25, z) + 3.0 * cyl_bessel_k(0.75, z) - cyl_bessel_k(1.25, z);
    return std::exp(-z) * std::pow(y / 2.0, 2.5) * b / std::sqrt(M_PI);
}

inline double psi1_term(int k, double x) {
    // Compute the k-th term of the psi1 series expansion in equation 1.10
    // of Csörgő, S. and Faraway, J. (1996).
    double m = 2.0 * k + 1.0;
    double sx = 2.0 * std::sqrt(x);
    double y1 = std::pow(x, 0.75);
    double y2 = std::pow(x, 1.25);
    double gamma_kp1_2 = gamma(k + 0.5);
    double gamma_kp3_2 = gamma(k + 1.5);
    double e1 = m * gamma_kp1_2 * psi1_G((4.0 * k + 3.0) / sx) / (9.0 * y1);
    double e2 = gamma_kp1_2 * psi1_H((4.0 * k + 1.0) / sx) / (72.0 * y2);
    double e3 = 2.0 * (m + 2.0) * gamma_kp3_2 * psi1_H((4.0 * k + 5.0) / sx) / (12.0 * y2);
    double e4 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k + 1.0) / sx) / (144.0 * y1);
    double e5 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k + 5.0) / sx) / (144.0 * y1);
    return e1 + e2 + e3 + e4 + e5;
}

inline double psi1_mod(double x) {
    // psi1 is defined in equation 1.10 in Csörgő, S. and Faraway, J. (1996).
    // This implements a modified version by excluding the term V(x) / 12
    // (here: cdf_cvm_inf(x) / 12) to avoid evaluating cdf_cvm_inf(x) twice in cdf_cvm.
    //
    // Implementation based on MAPLE code of Julian Faraway and R code of the
    // function pCvM in the package goftest (v1.1.1), permission granted
    // by Adrian Baddeley. Main difference in the implementation: the code
    // here keeps adding terms of the series until the terms are small enough.
    int k = 0;
    double tot = 0.0;
    while (true) {
        double gamma_kp1 = gamma(k + 1.0);
        double z = -psi1_term(k, x) / (M_PI * gamma_kp1);
        tot += z;
        if (std::abs(z) < 1e-7) {
            break;
        }
        ++k;
    }
    return tot;
}

inline double cdf_cvm(double x, int n = -1) {
    // Calculate the CDF of the Cramér-von Mises statistic for a finite sample
    // size n. If n=-1, use the asymptotic CDF (n=inf).
    //
    // See equation 1.8 in Csörgő, S. and Faraway, J. (1996) for finite samples,
    // 1.2 for the asymptotic CDF.
    //
    // For finite n, the approximation is less accurate near the support boundaries
    // [1/(12*n), n/3] and for larger x values where the CDF is close to 1.
    // The returned value is clipped to [0, 1].
    if (n == -1) {
        return cdf_cvm_inf(x);
    }
    double lower = 1.0 / (12.0 * n);
    double upper = static_cast<double>(n) / 3.0;
    if (x <= lower) {
        return 0.0;
    }
    if (x >= upper) {
        return 1.0;
    }
    double cdf = cdf_cvm_inf(x) * (1.0 + lower) + psi1_mod(x) / n;
    return std::min(std::max(cdf, 0.0), 1.0);
}

inline double fdtr(double a, double b, double x) { return cephes::fdtr(a, b, x); }

inline double fdtrc(double a, double b, double x) { return cephes::fdtrc(a, b, x); }

inline double fdtri(double a, double b, double y) { return cephes::fdtri(a, b, y); }

inline double gdtr(double a, double b, double x) { return cephes::gdtr(a, b, x); }

inline double gdtrc(double a, double b, double x) { return cephes::gdtrc(a, b, x); }

inline double kolmogorov(double x) { return cephes::kolmogorov(x); }

inline double kolmogc(double x) { return cephes::kolmogc(x); }

inline double kolmogi(double x) { return cephes::kolmogi(x); }

inline double kolmogci(double x) { return cephes::kolmogci(x); }

inline double kolmogp(double x) { return cephes::kolmogp(x); }

XSF_HOST_DEVICE inline double ndtr(double x) { return cephes::ndtr(x); }

XSF_HOST_DEVICE inline float ndtr(float x) { return ndtr(static_cast<double>(x)); }

inline std::complex<double> ndtr(std::complex<double> z) { return 0.5 * erfc(-z * M_SQRT1_2); }

inline std::complex<float> ndtr(std::complex<float> z) {
    return static_cast<std::complex<float>>(ndtr(static_cast<std::complex<double>>(z)));
}

/*
 * Log of the CDF of the normal distribution for double x.
 *
 * Let F(x) be the CDF of the standard normal distribution.
 * This implementation of log(F(x)) is based on the identities
 *
 *   F(x) = erfc(-x/√2)/2
 *        = 1 - erfc(x/√2)/2
 *
 * We use the first formula for x < -1, with erfc(z) replaced
 * by erfcx(z)*exp(-z**2) to ensure high precision for large
 * negative values when we take the logarithm:
 *
 *   log F(x) = log(erfc(-x/√2)/2)
 *            = log(erfcx(-x/√2)/2)*exp(-x**2/2))
 *            = log(erfcx(-x/√2)/2) - x**2/2
 *
 * For x >= -1, we use the second formula for F(x):
 *
 *   log F(x) = log(1 - erfc(x/√2)/2)
 *            = log1p(-erfc(x/√2)/2)
 */
inline double log_ndtr(double x) {
    double t = x * M_SQRT1_2;
    if (x < -1.0) {
        return log(erfcx(-t) / 2) - t * t;
    } else {
        return log1p(-erfc(t) / 2);
    }
}

inline float log_ndtr(float x) { return log_ndtr(static_cast<double>(x)); }

/*
 * Log of the normal CDF for complex arguments.
 *
 * This is equivalent to log(ndtr(z)), but is more robust to overflow at $z\to\infty$.
 * This implementation uses $\erfc(z) = \exp(-z^2) w(iz)$ taking special care to select
 * the principal branch of the log function log( exp(-z^2) w(i z) )
 */
inline std::complex<double> log_ndtr(std::complex<double> z) {
    if (z.real() > 6) {
        // Underflow. Close to the real axis, expand the log in log(1 - ndtr(-z)).
        std::complex<double> w = -0.5 * erfc(z * M_SQRT1_2);
        if (std::abs(w) < 1e-8) {
            return w;
        }
    }

    z *= -M_SQRT1_2;
    double x = std::real(z);
    double y = std::imag(z);

    /* Compute the principal branch of $log(exp(-z^2))$, using the fact that
     * $log(e^t) = log|e^t| + i Arg(e^t)$, and that if $t = r + is$, then
     * $e^t = e^r (\cos(s) + i \sin(s))$.
     */
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2 * x * y;        // Im(-z^2)

    double im = fmod(mIm_z2, 2.0 * M_PI);
    if (im > M_PI) {
        im -= 2.0 * M_PI;
    }

    std::complex<double> val1 = std::complex<double>(mRe_z2, im);

    std::complex<double> val2 = log(xsf::wofz(complex<double>(-y, x)));
    std::complex<double> result = val1 + val2 - M_LN2;

    /* Again, select the principal branch: log(z) = log|z| + i arg(z), thus
     * the imaginary part of the result should belong to [-pi, pi].
     */
    im = imag(result);
    if (im >= M_PI) {
        im -= 2 * M_PI;
    }
    if (im < -M_PI) {
        im += 2 * M_PI;
    }

    return {result.real(), im};
}

inline std::complex<float> log_ndtr(std::complex<float> z) {
    return static_cast<std::complex<float>>(log_ndtr(static_cast<std::complex<double>>(z)));
}

inline double nbdtr(int k, int n, double p) { return cephes::nbdtr(k, n, p); }

XSF_HOST_DEVICE inline bool bivariate_normal_cdf_boundary(double dh, double dk, double r, double &p) {
    // Handles degenerate cases for bivariate_normal_cdf (infinite arguments or zero correlation).
    // Returns true and sets p if a boundary case applies, false otherwise.

    const double math_inf = std::numeric_limits<double>::infinity();
    if (dh == math_inf || dk == math_inf) {
        // if dh ==  inf | dk ==  inf:p = 0;
        p = 0.0;
        return true;
    }
    if (dh == -math_inf) {
        // elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
        p = (dk == -math_inf) ? 1.0 : ndtr(-dk);
        return true;
    }
    if (dk == -math_inf) {
        // elseif dk == -inf, p = phid(-dh);
        p = ndtr(-dh);
        return true;
    }
    if (r == 0.0) {
        // elseif r == 0, p = phid(-dh)*phid(-dk);
        p = ndtr(-dh) * ndtr(-dk);
        return true;
    }
    return false;
}

constexpr double bvn_w6[3] = {0.1713244923791705, 0.3607615730481384, 0.4679139345726904};
constexpr double bvn_x6[3] = {0.9324695142031522, 0.6612093864662647, 0.2386191860831970};

constexpr double bvn_w12[6] = {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                               0.2031674267230659,  0.2334925365383547, 0.2491470458134029};
constexpr double bvn_x12[6] = {0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
                               0.5873179542866171, 0.3678314989981802, 0.1252334085114692};

constexpr double bvn_w20[10] = {0.01761400713915212, 0.04060142980038694, 0.06267204833410906, 0.08327674157670475,
                                0.1019301198172404,  0.1181945319615184,  0.1316886384491766,  0.1420961093183821,
                                0.1491729864726037,  0.1527533871307259};
constexpr double bvn_x20[10] = {0.9931285991850949, 0.9639719272779138, 0.9122344282513259, 0.8391169718222188,
                                0.7463319064601508, 0.6360536807265150, 0.5108670019508271, 0.3737060887154196,
                                0.2277858511416451, 0.07652652113349733};

XSF_HOST_DEVICE inline double bivariate_normal_cdf(double dh, double dk, double r) {
    // CDF of the bivariate normal distribution with zero means, unit variances,
    // and correlation r, evaluated at (dh, dk).
    //
    // dh, dk: the values at which to evaluate the CDF
    // r: the correlation coefficient, must be in [-1, 1]
    //
    // Adapted from the MATLAB original implementation by Dr. Alan Genz;
    // see license information in _qmvnt.py
    // In the comments, phid is the CDF of the standard normal distribution.

    double p;
    if (bivariate_normal_cdf_boundary(dh, dk, r, p)) {
        return p;
    }
    // else, tp = 2*pi; h = dh; k = dk; hk = h*k; bvn = 0;
    double tp = 2 * M_PI;
    double h = dh;
    double k = dk;
    double hk = h * k;
    double bvn = 0.0;
    const double *w;
    const double *x;
    int n;
    if (std::abs(r) < 0.3) {
        // Gauss Legendre points and weights, n = 6
        w = bvn_w6;
        x = bvn_x6;
        n = 3;
    } else if (std::abs(r) < 0.75) {
        // Gauss Legendre points and weights, n = 12
        w = bvn_w12;
        x = bvn_x12;
        n = 6;
    } else {
        // Gauss Legendre points and weights, n = 20
        w = bvn_w20;
        x = bvn_x20;
        n = 10;
    }

    // if abs(r) < 0.925, hs = ( h*h + k*k )/2; asr = asin(r)/2;
    if (std::abs(r) < 0.925) {
        double hs = (h * h + k * k) / 2.0;
        double asr = std::asin(r) / 2.0;
        // sn = sin(asr*x); bvn = exp((sn*hk-hs)./(1-sn.^2))*w';
        for (int i = 0; i < n; ++i) {
            const double x_i[2] = {1.0 - x[i], 1.0 + x[i]};
            for (int j = 0; j < 2; ++j) {
                double sn = std::sin(asr * x_i[j]);
                bvn += std::exp((sn * hk - hs) / (1.0 - sn * sn)) * w[i];
            }
        }
        // bvn = bvn*asr/tp + phid(-h)*phid(-k);
        bvn = bvn * asr / tp + ndtr(-h) * ndtr(-k);
    } else {
        // else, if r < 0, k = -k; hk = -hk; end
        if (r < 0.0) {
            k = -k;
            hk = -hk;
        }
        if (std::abs(r) < 1.0) {
            // if abs(r) < 1, as = 1-r^2; a = sqrt(as); bs = (h-k)^2;
            double as_ = 1.0 - r * r;
            double a = std::sqrt(as_);
            double bs = (h - k) * (h - k);
            // asr = -( bs/as + hk )/2; c = (4-hk)/8 ; d = (12-hk)/80;
            double asr = -(bs / as_ + hk) / 2.0;
            double c = (4.0 - hk) / 8.0;
            double d = (12.0 - hk) / 80.0;
            if (asr > -100.0)
                // if asr > -100, bvn = a*exp(asr)*(1-c*(bs-as)*(1-d*bs)/3+c*d*as^2); end
                bvn = a * std::exp(asr) * (1.0 - c * (bs - as_) * (1.0 - d * bs) / 3.0 + c * d * as_ * as_);
            if (hk > -100.0) {
                // if hk  > -100, b = sqrt(bs); sp = sqrt(tp)*phid(-b/a);
                double b = std::sqrt(bs);
                double sp = std::sqrt(tp) * ndtr(-b / a);
                // bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1-d*bs)/3 );
                bvn = bvn - std::exp(-hk / 2.0) * sp * b * (1.0 - c * bs * (1.0 - d * bs) / 3.0);
            }
            // end, a = a/2; xs = (a*x).^2; asr = -( bs./xs + hk )/2;
            a = a / 2.0;

            // ix = find( asr > -100 ); xs = xs(ix); sp = ( 1 + c*xs.*(1+5*d*xs) );
            double tmp = 0.0;
            for (int i = 0; i < n; ++i) {
                const double x_i[2] = {1.0 - x[i], 1.0 + x[i]};
                for (int j = 0; j < 2; ++j) {
                    double xs_i = (a * x_i[j]) * (a * x_i[j]);
                    double asr_i = -(bs / xs_i + hk) / 2.0;
                    if (asr_i <= -100.0) {
                        continue;
                    }
                    double sp = 1.0 + c * xs_i * (1.0 + 5.0 * d * xs_i);
                    // rs = sqrt(1-xs); ep = exp( -(hk/2)*xs./(1+rs).^2 )./rs;
                    double rs = std::sqrt(1.0 - xs_i);
                    double ep = std::exp(-(hk / 2.0) * xs_i / ((1.0 + rs) * (1.0 + rs))) / rs;
                    tmp += w[i] * std::exp(asr_i) * (sp - ep);
                }
            }

            // bvn = ( a*( (exp(asr(ix)).*(sp-ep))*w(ix)' ) - bvn )/tp;
            bvn = (a * tmp - bvn) / tp;
        }
        // end
        // if r > 0, bvn =  bvn + phid( -max( h, k ) );
        if (r > 0.0) {
            bvn = bvn + ndtr(-(h > k ? h : k));
        } else if (h >= k) {
            // elseif h >= k, bvn = -bvn;
            bvn = -bvn;
        } else {
            // else, if h < 0, L = phid(k)-phid(h); else, L = phid(-h)-phid(-k); end
            double L;
            if (h < 0.0) {
                L = ndtr(k) - ndtr(h);
            } else {
                L = ndtr(-h) - ndtr(-k);
            }
            // bvn =  L - bvn;
            bvn = L - bvn;
        }
    }
    // end, p = max( 0, min( 1, bvn ) );
    return bvn < 0.0 ? 0.0 : (bvn > 1.0 ? 1.0 : bvn);
}

inline double nbdtrc(int k, int n, double p) { return cephes::nbdtrc(k, n, p); }

inline double nbdtri(int k, int n, double p) { return cephes::nbdtri(k, n, p); }

inline double ndtri(double x) { return cephes::ndtri(x); }

inline float ndtri(float x) { return static_cast<float>(cephes::ndtri(x)); }

XSF_HOST_DEVICE inline double nrdtrimn(double p, double std, double x) {
    if (std::isnan(std) || std <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isnan(p) || p <= 0 || p >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return x - std * cephes::ndtri(p);
}

XSF_HOST_DEVICE inline float nrdtrimn(float p, float std, float x) {
    return static_cast<float>(nrdtrimn(static_cast<double>(p), static_cast<double>(std), static_cast<double>(x)));
}

XSF_HOST_DEVICE inline double nrdtrisd(double mean, double p, double x) {
    if (std::isnan(mean) || std::isnan(p) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (p <= 0 || p >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (x - mean) / cephes::ndtri(p);
}

XSF_HOST_DEVICE inline float nrdtrisd(float mean, float p, float x) {
    return static_cast<float>(nrdtrisd(static_cast<double>(mean), static_cast<double>(p), static_cast<double>(x)));
}

inline double owens_t(double h, double a) { return cephes::owens_t(h, a); }

inline double pdtr(double k, double m) { return cephes::pdtr(k, m); }

inline double pdtrc(double k, double m) { return cephes::pdtrc(k, m); }

inline double pdtri(int k, double y) { return cephes::pdtri(k, y); }

namespace detail {

    template <typename FreqTable2D>
    XSF_HOST_DEVICE inline void
    cvm_freq_table_all(int64_t m, int64_t n, int64_t a, int64_t b, FreqTable2D gs, FreqTable2D next_gs) {
        using T = typename FreqTable2D::value_type;
        int64_t K = static_cast<int64_t>(gs.extent(1));

        // initialize gs to 0
        for (int64_t v = 0; v < m + 1; ++v)
            for (int64_t k = 0; k < K; ++k)
                gs(v, k) = T(0);
        // base case: gs(0, 0) = 1
        gs(0, 0) = T(1);

        for (int64_t u = 0; u < n + 1; ++u) {
            // v = 0: no next_gs(v-1, ...) term
            {
                int64_t d = -b * u;
                int64_t d2 = d * d;
                int64_t kstart = (d2 < K) ? d2 : K;
                // next_gs(0, k) = gs(0, k - d2) for k >= d2, else 0
                for (int64_t k = 0; k < kstart; ++k) {
                    next_gs(0, k) = T(0);
                }
                for (int64_t k = kstart; k < K; ++k) {
                    next_gs(0, k) = gs(0, k - d2);
                }
            }
            // v > 0: both terms contribute
            for (int64_t v = 1; v < m + 1; ++v) {
                int64_t d = a * v - b * u;
                int64_t d2 = d * d; // d^2 = (a*v - b*u)^2
                int64_t kstart = (d2 < K) ? d2 : K;
                for (int64_t k = 0; k < kstart; ++k) {
                    next_gs(v, k) = T(0);
                }
                for (int64_t k = kstart; k < K; ++k) {
                    next_gs(v, k) = next_gs(v - 1, k - d2) + gs(v, k - d2);
                }
            }
            FreqTable2D tmp = gs;
            gs = next_gs;
            next_gs = tmp;
        }
        // We swap `gs` and `next_gs` at each u-step, so buffer parity depends on n.
        // If n is even, the final table ends up in the original `next_gs` buffer;
        // copy it back so the caller can always read results from the original `gs`.
        if (n % 2 == 0) {
            for (int64_t v = 0; v < m + 1; ++v) {
                for (int64_t k = 0; k < K; ++k) {
                    next_gs(v, k) = gs(v, k);
                }
            }
        }
    }

} // namespace detail

template <typename FreqTable2D>
XSF_HOST_DEVICE inline double
pval_cvm_2samp_exact(double s, int64_t m, int64_t n, FreqTable2D gs, FreqTable2D next_gs) {
    /*
     * Compute the exact p-value of the Cramér-von Mises two-sample test
     * for a given value s of the test statistic and where m and n are the sizes
     * of the samples.
     *
     * [1] Y. Xiao, A. Gordon, and A. Yakovlev, "A C++ Program for
     *     the Cramér-Von Mises Two-Sample Test", J. Stat. Soft.,
     *     vol. 17, no. 8, pp. 1-15, Dec. 2006.
     * [2] T. W. Anderson "On the Distribution of the Two-Sample Cramér-von Mises
     *     Criterion," The Annals of Mathematical Statistics, Ann. Math. Statist.
     *     33(3), 1148-1159, (September, 1962)
     */
    if (m <= 0 || n <= 0) {
        set_error("pval_cvm_2samp_exact", SF_ERROR_DOMAIN, "m and n must be positive");
        return std::numeric_limits<double>::quiet_NaN();
    }
    // [1, p. 3]
    int64_t lcm = std::lcm(m, n);
    // [1, p. 4], below eq. 3
    int64_t a = lcm / m;
    int64_t b = lcm / n;
    // Combine Eq. 9 in [2] with Eq. 2 in [1] and solve for $\zeta$
    // Hint: `s` is $U$ in [2], and $T_2$ in [1] is $T$ in [2]
    int64_t mn = m * n;

    // Uses double floor division since s is double
    int64_t zeta =
        static_cast<int64_t>(std::floor((lcm * lcm * (m + n) * (6.0 * s - mn * (4.0 * mn - 1))) / (6.0 * mn * mn)));

    detail::cvm_freq_table_all(m, n, a, b, gs, next_gs);

    int64_t K = static_cast<int64_t>(gs.extent(1));

    // Clamp to prevent negative indexing when zeta < 0.
    int64_t k0 = (zeta < 0) ? 0 : zeta;

    int64_t sum_freq = 0;
    for (int64_t k = k0; k < K; ++k) {
        sum_freq += gs(m, k);
    }

    double combinations = xsf::binom(static_cast<double>(m + n), static_cast<double>(m));
    return sum_freq / combinations;
}

inline double smirnov(int n, double x) { return cephes::smirnov(n, x); }

inline double smirnovc(int n, double x) { return cephes::smirnovc(n, x); }

inline double smirnovi(int n, double x) { return cephes::smirnovi(n, x); }

inline double smirnovci(int n, double x) { return cephes::smirnovci(n, x); }

inline double smirnovp(int n, double x) { return cephes::smirnovp(n, x); }

inline double tukeylambdacdf(double x, double lmbda) { return cephes::tukeylambdacdf(x, lmbda); }

inline float tukeylambdacdf(float x, double lmbda) {
    return tukeylambdacdf(static_cast<double>(x), static_cast<double>(lmbda));
}

namespace detail {

    template <typename InputMat, typename OutputMat>
    XSF_HOST_DEVICE inline void poisson_binom_pmf_all_impl(InputMat p, OutputMat res) {
        using T = typename OutputMat::value_type;
        auto n = p.extent(0);

        // Initialize output array with zeros.
        for (decltype(n) i = 0; i < n + 1; ++i) {
            res(i) = T(0.0);
        }

        if (n == 0) {
            res(0) = T(1.0);
            return;
        }

        res(0) = T(1.0) - p(0);
        res(1) = p(0);

        for (decltype(n) i = 1; i < n; i++) {
            T p_i = p(i);
            T q_i = 1 - p_i;
            for (decltype(n) j = i + 1; j >= 1; j--) {
                T tmp = res(j - 1) * p_i;
                res(j - 1) *= q_i;
                res(j) += tmp;
            }
        }
    }

} // namespace detail

template <typename InputMat, typename OutputMat>
XSF_HOST_DEVICE inline void poisson_binom_pmf_all(InputMat p, OutputMat res) {
    /* Outputs the full pmf of a Poisson-Binomial distribution
     *
     * p should be an std::mdspan view of a 1d array of length n containing
     * the success probabilities for the n Bernoulli trials. res should be
     * a std::mdspan view of a 1d array of of length n + 1. It is up to
     * the caller to pass valid values for p and res.
     *
     * Upon completion, res(k) will contain the probability of observing k
     * successes for k from 0 to n.
     */
    auto n = p.extent(0);
    auto out_size = res.extent(0);

    if (out_size != n + 1) {
        set_error("_poisson_binom_pmf_all", SF_ERROR_MEMORY, "out.shape[-1] must be p.shape[-1] + 1");
        return;
    }

    detail::poisson_binom_pmf_all_impl(p, res);
}

template <typename InputMat, typename OutputMat>
XSF_HOST_DEVICE inline void poisson_binom_cdf_all(InputMat p, OutputMat res) {
    /* Outputs the full cdf of a Poisson-Binomial distribution */
    using T = typename OutputMat::value_type;
    auto n = p.extent(0);
    auto out_size = res.extent(0);

    if (out_size != n + 1) {
        set_error("_poisson_binom_cdf_all", SF_ERROR_MEMORY, "out.shape[-1] must be p.shape[-1] + 1");
        return;
    }

    detail::poisson_binom_pmf_all_impl(p, res);
    for (decltype(n) i = 1; i < n; i++) {
        res(i) = std::min(res(i) + res(i - 1), T(1.0));
    }
    res(n) = T(1.0);
}

template <typename InputMat>
XSF_HOST_DEVICE inline typename InputMat::value_type take_from_pmf(InputMat pmf, long long int k) {
    using T = typename InputMat::value_type;
    auto size = pmf.extent(0);
    if ((k < 0) || (k >= static_cast<long long int>(size))) {
        return T(0.0);
    }
    return pmf(k);
}

template <typename InputMat>
XSF_HOST_DEVICE inline typename InputMat::value_type take_from_discrete_cdf(InputMat cdf, long long int k) {
    using T = typename InputMat::value_type;
    auto size = cdf.extent(0);
    if (k < 0) {
        return T(0.0);
    }
    if (k >= static_cast<long long int>(size) - 1) {
        return T(1.0);
    }
    return cdf(k);
}

namespace detail {

    XSF_HOST_DEVICE inline double von_mises_cdf_series(double k, double x, unsigned int p) {
        double s, c, sn, cn, r, v;
        unsigned int n;
        s = std::sin(x);
        c = std::cos(x);
        sn = std::sin(p * x);
        cn = std::cos(p * x);
        r = 0;
        v = 0;
        for (n = p - 1; n > 0; --n) {
            double sn_new = sn * c - cn * s;
            double cn_new = cn * c + sn * s;
            sn = sn_new;
            cn = cn_new;
            r = k / (2 * n + k * r);
            v = r * (sn / n + v);
        }
        return 0.5 + x / (2.0 * M_PI) + v / M_PI;
    }

    XSF_HOST_DEVICE inline double von_mises_cdf_normalapprox(double k, double x) {
        double b = xsf::cephes::detail::SQRT2OPI / cephes::i0e(k);
        double z = b * std::sin(x / 2.0);
        return ndtr(z);
    }

} // namespace detail

XSF_HOST_DEVICE inline double von_mises_cdf(double k, double x) {
    // CDF of the von Mises distribution with concentration k, extended
    // periodically over x.
    //
    // For k < 50, this uses the backward-recursion series method of Hill [1],
    // Algorithm 518, with constants chosen for about 12 decimal digits. For k >= 50,
    // this switches to a large-concentration normal approximation.
    //
    // References:
    // [1] G. W. Hill, "Algorithm 518: Incomplete Bessel Function I0.
    //     The Von Mises Distribution [S14]", ACM Transactions on Mathematical
    //     Software, 3(3), 279-284, 1977.
    //     DOI: https://doi.org/10.1145/355744.355753

    if (k < 0) {
        set_error("von_mises_cdf", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }

    double ix = std::round(x / (2 * M_PI));
    x -= ix * 2.0 * M_PI;

    // These values should give 12 decimal digits
    const double ck = 50.0;
    const double a1 = 28.0, a2 = 0.5, a3 = 100.0, a4 = 5.0;
    double result;
    if (k < ck) {
        unsigned int p = static_cast<unsigned int>(1 + a1 + a2 * k - a3 / (k + a4));
        result = detail::von_mises_cdf_series(k, x, p);
        result = std::min(std::max(result, 0.0), 1.0);
    } else {
        result = detail::von_mises_cdf_normalapprox(k, x);
    }
    return result + ix;
}

XSF_HOST_DEVICE inline float von_mises_cdf(float k, float x) {
    return static_cast<float>(von_mises_cdf(static_cast<double>(k), static_cast<double>(x)));
}

} // namespace xsf
