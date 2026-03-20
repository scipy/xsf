#pragma once

#include "xsf/bessel.h"
#include "xsf/cephes/bdtr.h"
#include "xsf/cephes/chdtr.h"
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

inline double ndtr(double x) { return cephes::ndtr(x); }

inline float ndtr(float x) { return ndtr(static_cast<double>(x)); }

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

inline double bvnu(double dh, double dk, double r) {
    // CDF of the bivariate normal distribution with zero means, unit variances,
    // and correlation r, evaluated at (dh, dk).
    //
    // dh, dk: the values at which to evaluate the CDF
    // r: the correlation coefficient, must be in [-1, 1]
    //
    // Adapted from the MATLAB original implementation by Dr. Alan Genz;
    // see license information in _qmvnt.py
    // In the comments, phid is the CDF of the standard normal distribution.
    const double math_inf = std::numeric_limits<double>::infinity();
    double p;
    if (dh == math_inf || dk == math_inf) {
        // if dh ==  inf | dk ==  inf:p = 0;
        p = 0.0;
    } else if (dh == -math_inf) {
        // elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
        if (dk == -math_inf) {
            p = 1.0;
        } else {
            p = ndtr(-dk);
        }
    } else if (dk == -math_inf) {
        // elseif dk == -inf, p = phid(-dh);
        p = ndtr(-dh);
    } else if (r == 0.0) {
        // elseif r == 0, p = phid(-dh)*phid(-dk);
        p = ndtr(-dh) * ndtr(-dk);
    } else {
        // else, tp = 2*pi; h = dh; k = dk; hk = h*k; bvn = 0;
        double tp = 2 * M_PI;
        double h = dh;
        double k = dk;
        double hk = h * k;
        double bvn = 0.0;
        std::vector<double> w, x;

        if (std::abs(r) < 0.3) {
            // if abs(r) < 0.3      % Gauss Legendre points and weights, n =  6
            //     w(1:3) = [0.1713244923791705 0.3607615730481384 0.4679139345726904];
            //     x(1:3) = [0.9324695142031522 0.6612093864662647 0.2386191860831970];

            w = {0.1713244923791705, 0.3607615730481384, 0.4679139345726904};
            x = {0.9324695142031522, 0.6612093864662647, 0.2386191860831970};

        } else if (std::abs(r) < 0.75) {
            //  elseif abs(r) < 0.75 % Gauss Legendre points and weights, n = 12
            //      w(1:3) = [.04717533638651177 0.1069393259953183 0.1600783285433464];
            //      w(4:6) = [0.2031674267230659 0.2334925365383547 0.2491470458134029];
            //      x(1:3) = [0.9815606342467191 0.9041172563704750 0.7699026741943050];
            //      x(4:6) = [0.5873179542866171 0.3678314989981802 0.1252334085114692];

            w = {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                 0.2031674267230659,  0.2334925365383547, 0.2491470458134029};
            x = {0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
                 0.5873179542866171, 0.3678314989981802, 0.1252334085114692};
        } else {
            // else,                % Gauss Legendre points and weights, n = 20
            //      w(1:3) = [.01761400713915212 .04060142980038694 .06267204833410906];
            //      w(4:6) = [.08327674157670475 0.1019301198172404 0.1181945319615184];
            //      w(7:9) = [0.1316886384491766 0.1420961093183821 0.1491729864726037];
            //      w(10) =   0.1527533871307259;
            //      x(1:3) = [0.9931285991850949 0.9639719272779138 0.9122344282513259];
            //      x(4:6) = [0.8391169718222188 0.7463319064601508 0.6360536807265150];
            //      x(7:9) = [0.5108670019508271 0.3737060887154196 0.2277858511416451];
            //      x(10) =   0.07652652113349733;

            w = {0.01761400713915212, 0.04060142980038694, 0.06267204833410906, 0.08327674157670475,
                 0.1019301198172404,  0.1181945319615184,  0.1316886384491766,  0.1420961093183821,
                 0.1491729864726037,  0.1527533871307259};
            x = {0.9931285991850949, 0.9639719272779138, 0.9122344282513259, 0.8391169718222188, 0.7463319064601508,
                 0.6360536807265150, 0.5108670019508271, 0.3737060887154196, 0.2277858511416451, 0.07652652113349733};
        }
        // end, w = [w  w]; x = [1-x 1+x];

        std::vector<double> w_new = w;
        for (size_t i = 0; i < w.size(); ++i) {
            w_new.push_back(w[i]);
        }
        w = w_new;

        std::vector<double> x_new;
        for (size_t i = 0; i < x.size(); ++i) {
            x_new.push_back(1.0 - x[i]);
        }
        for (size_t i = 0; i < x.size(); ++i) {
            x_new.push_back(1.0 + x[i]);
        }
        x = x_new;

        // if abs(r) < 0.925, hs = ( h*h + k*k )/2; asr = asin(r)/2;
        if (std::abs(r) < 0.925) {
            double hs = (h * h + k * k) / 2.0;
            double asr = std::asin(r) / 2.0;
            // sn = sin(asr*x); bvn = exp((sn*hk-hs)./(1-sn.^2))*w';
            std::vector<double> sn(x.size());
            for (size_t i = 0; i < x.size(); ++i) {
                sn[i] = std::sin(asr * x[i]);
                bvn += std::exp((sn[i] * hk - hs) / (1.0 - sn[i] * sn[i])) * w[i];
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
                std::vector<double> xs_sel, w_sel, asr_sel;
                for (size_t i = 0; i < x.size(); ++i) {
                    double xs_i = (a * x[i]) * (a * x[i]);
                    double asr_i = -(bs / xs_i + hk) / 2.0;
                    if (asr_i > -100.0) {
                        xs_sel.push_back(xs_i);
                        w_sel.push_back(w[i]);
                        asr_sel.push_back(asr_i);
                    }
                }

                double sp, rs, ep, tmp = 0.0;
                for (size_t i = 0; i < xs_sel.size(); ++i) {
                    sp = 1.0 + c * xs_sel[i] * (1.0 + 5.0 * d * xs_sel[i]);
                    // rs = sqrt(1-xs); ep = exp( -(hk/2)*xs./(1+rs).^2 )./rs;
                    rs = std::sqrt(1.0 - xs_sel[i]);
                    ep = std::exp(-(hk / 2.0) * xs_sel[i] / ((1.0 + rs) * (1.0 + rs))) / rs;
                    tmp += w_sel[i] * std::exp(asr_sel[i]) * (sp - ep);
                }

                // bvn = ( a*( (exp(asr(ix)).*(sp-ep))*w(ix)' ) - bvn )/tp;
                bvn = (a * tmp - bvn) / tp;
            }
            // end
            // if r > 0, bvn =  bvn + phid( -max( h, k ) );
            if (r > 0.0) {
                bvn = bvn + ndtr(-std::max(h, k));
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
        // end
        // end, p = max( 0, min( 1, bvn ) );
        p = std::max(0.0, std::min(1.0, bvn));
    }
    // end
    return p;
}

inline double nbdtrc(int k, int n, double p) { return cephes::nbdtrc(k, n, p); }

inline double nbdtri(int k, int n, double p) { return cephes::nbdtri(k, n, p); }

inline double ndtri(double x) { return cephes::ndtri(x); }

inline float ndtri(float x) { return static_cast<float>(cephes::ndtri(x)); }

inline double owens_t(double h, double a) { return cephes::owens_t(h, a); }

inline double pdtr(double k, double m) { return cephes::pdtr(k, m); }

inline double pdtrc(double k, double m) { return cephes::pdtrc(k, m); }

inline double pdtri(int k, double y) { return cephes::pdtri(k, y); }

inline double smirnov(int n, double x) { return cephes::smirnov(n, x); }

inline double smirnovc(int n, double x) { return cephes::smirnovc(n, x); }

inline double smirnovi(int n, double x) { return cephes::smirnovi(n, x); }

inline double smirnovci(int n, double x) { return cephes::smirnovci(n, x); }

inline double smirnovp(int n, double x) { return cephes::smirnovp(n, x); }

inline double tukeylambdacdf(double x, double lmbda) { return cephes::tukeylambdacdf(x, lmbda); }

inline float tukeylambdacdf(float x, double lmbda) {
    return tukeylambdacdf(static_cast<double>(x), static_cast<double>(lmbda));
}

} // namespace xsf
