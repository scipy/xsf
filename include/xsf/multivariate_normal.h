#pragma once

#include "cephes/ndtr.h"
#include "config.h"

namespace xsf {

namespace detail {

    XSF_HOST_DEVICE inline bool bivariate_normal_sf_boundary(double dh, double dk, double r, double &p) {
        // Handles degenerate cases for bivariate_normal_sf (infinite arguments or zero correlation).
        // Returns true and sets p if a boundary case applies, false otherwise.

        const double math_inf = std::numeric_limits<double>::infinity();
        if (dh == math_inf || dk == math_inf) {
            // if dh ==  inf | dk ==  inf:p = 0;
            p = 0.0;
            return true;
        }
        if (dh == -math_inf) {
            // elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
            p = (dk == -math_inf) ? 1.0 : cephes::ndtr(-dk);
            return true;
        }
        if (dk == -math_inf) {
            // elseif dk == -inf, p = phid(-dh);
            p = cephes::ndtr(-dh);
            return true;
        }
        if (r == 0.0) {
            // elseif r == 0, p = phid(-dh)*phid(-dk);
            p = cephes::ndtr(-dh) * cephes::ndtr(-dk);
            return true;
        }
        return false;
    }

    // Positive nodes and weights of symmetric Gauss-Legendre rules on [-1, 1].
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

} // namespace detail

XSF_HOST_DEVICE inline double bivariate_normal_sf(double dh, double dk, double r) {
    // Survival function of the standard bivariate normal distribution.
    //
    // Return P(X > dh, Y > dk) for a standard bivariate normal vector
    // (X, Y) with correlation r.
    //
    // dh, dk are the lower limits of the upper tail, and r must satisfy
    // -1 <= r <= 1.
    //
    // Adapted from the MATLAB original implementation by Dr. Alan Genz;
    // see license information in _qmvnt.py
    // In the comments, phid is the CDF of the standard normal distribution.

    double p;
    if (detail::bivariate_normal_sf_boundary(dh, dk, r, p)) {
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
        w = detail::bvn_w6;
        x = detail::bvn_x6;
        n = 3;
    } else if (std::abs(r) < 0.75) {
        // Gauss Legendre points and weights, n = 12
        w = detail::bvn_w12;
        x = detail::bvn_x12;
        n = 6;
    } else {
        // Gauss Legendre points and weights, n = 20
        w = detail::bvn_w20;
        x = detail::bvn_x20;
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
        bvn = bvn * asr / tp + cephes::ndtr(-h) * cephes::ndtr(-k);
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
                double sp = std::sqrt(tp) * cephes::ndtr(-b / a);
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
            bvn = bvn + cephes::ndtr(-(h > k ? h : k));
        } else if (h >= k) {
            // elseif h >= k, bvn = -bvn;
            bvn = -bvn;
        } else {
            // else, if h < 0, L = phid(k)-phid(h); else, L = phid(-h)-phid(-k); end
            double L;
            if (h < 0.0) {
                L = cephes::ndtr(k) - cephes::ndtr(h);
            } else {
                L = cephes::ndtr(-h) - cephes::ndtr(-k);
            }
            // bvn =  L - bvn;
            bvn = L - bvn;
        }
    }
    // end, p = max( 0, min( 1, bvn ) );
    return bvn < 0.0 ? 0.0 : (bvn > 1.0 ? 1.0 : bvn);
}

XSF_HOST_DEVICE inline float bivariate_normal_sf(float dh, float dk, float r) {
    return bivariate_normal_sf(static_cast<double>(dh), static_cast<double>(dk), static_cast<double>(r));
}

} // namespace xsf
