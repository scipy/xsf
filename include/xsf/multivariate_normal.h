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

XSF_HOST_DEVICE inline double bivariate_normal_cdf(double dh, double dk, double r) {
    // Cumulative distribution function of the standard bivariate normal distribution.
    //
    // Return P(X <= dh, Y <= dk) for a standard bivariate normal vector
    // (X, Y) with correlation r.
    //
    // dh, dk are the upper limits of the lower tail, and r must satisfy
    // -1 <= r <= 1.
    //
    // Adapted from the original MATLAB survival function implementation by Dr. Alan Genz;
    // see license information in LICENSES_bundled.txt.
    // In the comments, phid is the CDF of the standard normal distribution.

    // By symmetry, the CDF at (dh, dk) is the upper-tail probability at (-dh, -dk).
    dh = -dh;
    dk = -dk;
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

XSF_HOST_DEVICE inline float bivariate_normal_cdf(float dh, float dk, float r) {
    return bivariate_normal_cdf(static_cast<double>(dh), static_cast<double>(dk), static_cast<double>(r));
}

namespace detail {

    XSF_HOST_DEVICE inline void sin_cos2(double x, double &sx, double &cs) {
        // Computes sin(x) and cos(x)^2 with series approximation for |x| near pi/2
        double ee = (M_PI_2 - std::abs(x)) * (M_PI_2 - std::abs(x));
        if (ee < 5e-5) {
            cs = ee * (1.0 - ee * (1.0 - 2.0 * ee / 15.0) / 3.0);
            double sgn = (x > 0.0) - (x < 0.0);
            sx = (1.0 - ee * (1.0 - ee / 12.0) / 2.0) * sgn;
        } else {
            sx = std::sin(x);
            cs = 1.0 - sx * sx;
        }
    }

    XSF_HOST_DEVICE inline double
    plackett_term(double ba, double bb, double bc, double ra, double rb, double r, double rr) {
        // Evaluate one Plackett formula integrand term.
        double f = 0.0;
        double dt = rr * (rr - (ra - rb) * (ra - rb) - 2.0 * ra * rb * (1.0 - r));
        if (dt > 0.0) {
            double bt = (bc * rr + ba * (r * rb - ra) + bb * (r * ra - rb)) / std::sqrt(dt);
            double ft = (ba - r * bb) * (ba - r * bb) / rr + bb * bb;
            if (bt > -10.0 && ft < 100.0) {
                f = std::exp(-ft / 2.0);
                if (bt < 10.0) {
                    f *= cephes::ndtr(bt);
                }
            }
        }
        return f;
    }

    XSF_HOST_DEVICE inline double
    tvn_integrand(double x, double h1, double h2, double h3, double r23, double a12, double a13) {
        // Combine the two Plackett integrand terms at integration point x.
        double f = 0.0;
        double r12, rr2, r13, rr3;
        sin_cos2(a12 * x, r12, rr2);
        sin_cos2(a13 * x, r13, rr3);
        if (std::abs(a12) > 0.0) {
            f += a12 * plackett_term(h1, h2, h3, r13, r23, r12, rr2);
        }
        if (std::abs(a13) > 0.0) {
            f += a13 * plackett_term(h1, h3, h2, r12, r23, r13, rr3);
        }
        return f;
    }

    XSF_HOST_DEVICE inline void tvn_gauss_kronrod(
        double a, double b, double h1, double h2, double h3, double r23, double a12, double a13, double &resk,
        double &err
    ) {
        // Apply a Gauss-Kronrod rule on one adaptive subinterval.
        constexpr double wg0 = 0.2729250867779007;
        constexpr double wg[5] = {
            0.05566856711617449, 0.1255803694649048, 0.1862902109277352, 0.2331937645919914, 0.2628045445102478
        };
        constexpr double xgk[11] = {0.9963696138895427, 0.9782286581460570, 0.9416771085780681, 0.8870625997680953,
                                    0.8160574566562211, 0.7301520055740492, 0.6305995201619651, 0.5190961292068118,
                                    0.3979441409523776, 0.2695431559523450, 0.1361130007993617};
        constexpr double wgk0 = 0.1365777947111183;
        constexpr double wgk[11] = {0.00976544104596129, 0.02715655468210443, 0.04582937856442671, 0.06309742475037484,
                                    0.07866457193222764, 0.09295309859690074, 0.1058720744813894,  0.1167395024610472,
                                    0.1251587991003195,  0.1312806842298057,  0.1351935727998845};

        double wid = (b - a) / 2.0;
        double cen = (b + a) / 2.0;
        double fc = tvn_integrand(cen, h1, h2, h3, r23, a12, a13);
        double resg = fc * wg0;
        resk = fc * wgk0;

        for (int j = 0; j < 5; ++j) {
            double t = wid * xgk[2 * j];
            fc = tvn_integrand(cen - t, h1, h2, h3, r23, a12, a13) + tvn_integrand(cen + t, h1, h2, h3, r23, a12, a13);
            resk += wgk[2 * j] * fc;

            t = wid * xgk[2 * j + 1];
            fc = tvn_integrand(cen - t, h1, h2, h3, r23, a12, a13) + tvn_integrand(cen + t, h1, h2, h3, r23, a12, a13);
            resk += wgk[2 * j + 1] * fc;
            resg += wg[j] * fc;
        }

        double t = wid * xgk[10];
        fc = tvn_integrand(cen - t, h1, h2, h3, r23, a12, a13) + tvn_integrand(cen + t, h1, h2, h3, r23, a12, a13);
        resk = wid * (resk + wgk[10] * fc);
        err = std::abs(resk - wid * resg);
    }

    XSF_HOST_DEVICE inline double
    tvn_adaptive_integral(double h1, double h2, double h3, double r23, double a12, double a13, double tol) {
        // Adaptively integrate the Plackett integrand over [0, 1].
        constexpr int nl = 100;
        double ai[nl];
        double bi[nl];
        double fi[nl];
        double ei[nl];
        ai[0] = 0.0;
        bi[0] = 1.0;
        fi[0] = 0.0;
        ei[0] = 0.0;

        int ip = 0;
        int im = 0;
        double err = 1.0;
        double fin = 0.0;
        while (4.0 * err > tol && im < nl - 1) {
            ++im;
            bi[im] = bi[ip];
            ai[im] = (ai[ip] + bi[ip]) / 2.0;
            bi[ip] = ai[im];
            tvn_gauss_kronrod(ai[ip], bi[ip], h1, h2, h3, r23, a12, a13, fi[ip], ei[ip]);
            tvn_gauss_kronrod(ai[im], bi[im], h1, h2, h3, r23, a12, a13, fi[im], ei[im]);

            fin = 0.0;
            double err2 = 0.0;
            ip = 0;
            for (int i = 0; i <= im; ++i) {
                fin += fi[i];
                err2 += ei[i] * ei[i];
                if (ei[i] > ei[ip]) {
                    ip = i;
                }
            }
            err = std::sqrt(err2);
        }
        return fin;
    }

} // namespace detail

XSF_HOST_DEVICE inline double
trivariate_normal_cdf(double h1, double h2, double h3, double r12, double r13, double r23, double epsi = 1e-7) {
    // Cumulative distribution function of the standard trivariate normal distribution.
    //
    // Return P(X <= h1, Y <= h2, Z <= h3) for a standard trivariate normal
    // vector (X, Y, Z) with correlation matrix
    //
    //     [ 1   r12  r13 ]
    //     [ r12  1   r23 ]
    //     [ r13 r23   1  ].
    //
    // h1, h2, h3 are the upper integration limits. The correlations must
    // describe a positive semidefinite correlation matrix. epsi is the requested
    // absolute accuracy (default 1e-7), with a lower bound of 1e-14.
    //
    // Port of tvnl in Alan Genz's MATLAB tvn.m. See LICENSES_bundled.txt for the license.
    // Original source:
    // https://web.archive.org/web/20200205123040/http://www.math.wsu.edu/faculty/genz/software/matlab/tvn.m
    const double math_inf = std::numeric_limits<double>::infinity();
    double epst = epsi > 1e-14 ? epsi : 1e-14;

    if (std::isnan(h1) || std::isnan(h2) || std::isnan(h3) || std::isnan(r12) || std::isnan(r13) || std::isnan(r23) ||
        std::isnan(epsi)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double det = 1.0 + 2.0 * r12 * r13 * r23 - r12 * r12 - r13 * r13 - r23 * r23;
    // Tolerate tiny negative determinants caused by floating-point roundoff
    // when the correlation matrix is singular.
    constexpr double det_tol = 8.0 * std::numeric_limits<double>::epsilon();
    if (std::abs(r12) > 1.0 || std::abs(r13) > 1.0 || std::abs(r23) > 1.0 || det < -det_tol) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (h1 == -math_inf || h2 == -math_inf || h3 == -math_inf) {
        return 0.0;
    }
    if (h1 == math_inf) {
        if (h2 == math_inf) {
            return h3 == math_inf ? 1.0 : cephes::ndtr(h3);
        }
        if (h3 == math_inf) {
            return cephes::ndtr(h2);
        }
        return bivariate_normal_cdf(h2, h3, r23);
    }
    if (h2 == math_inf) {
        if (h3 == math_inf) {
            return cephes::ndtr(h1);
        }
        return bivariate_normal_cdf(h1, h3, r13);
    }
    if (h3 == math_inf) {
        return bivariate_normal_cdf(h1, h2, r12);
    }

    // Sort correlations and check special cases.
    if (std::abs(r12) > std::abs(r13)) {
        double h2_old = h2;
        double r12_old = r12;
        h2 = h3;
        h3 = h2_old;
        r12 = r13;
        r13 = r12_old;
    }
    if (std::abs(r13) > std::abs(r23)) {
        double h1_old = h1;
        double r13_old = r13;
        double r23_old = r23;
        h1 = h2;
        h2 = h1_old;
        r23 = r13_old;
        r13 = r23_old;
    }

    double tvn;
    if (std::abs(h1) + std::abs(h2) + std::abs(h3) < epst) {
        tvn = (1.0 + 2.0 * (std::asin(r12) + std::asin(r13) + std::asin(r23)) / M_PI) / 8.0;
    } else if (std::abs(r12) + std::abs(r13) < epst) {
        tvn = cephes::ndtr(h1) * bivariate_normal_cdf(h2, h3, r23);
    } else if (std::abs(r13) + std::abs(r23) < epst) {
        tvn = cephes::ndtr(h3) * bivariate_normal_cdf(h1, h2, r12);
    } else if (std::abs(r12) + std::abs(r23) < epst) {
        tvn = cephes::ndtr(h2) * bivariate_normal_cdf(h1, h3, r13);
    } else if (1.0 - r23 < epst) {
        tvn = bivariate_normal_cdf(h1, (h2 < h3 ? h2 : h3), r12);
    } else if (r23 == -1.0 && h2 <= -h3) {
        // When Z = -Y, the constraints require -h3 <= Y <= h2.
        // For h2 <= -h3, this event has probability zero; avoid quadrature.
        tvn = 0.0;
    } else if (r23 + 1.0 < epst && h2 > -h3) {
        tvn = bivariate_normal_cdf(h1, h2, r12) - bivariate_normal_cdf(h1, -h3, r12);
    } else {
        // Use numerical integration to compute probability and add
        // singular values from the Plackett formula.
        double a12 = std::asin(r12);
        double a13 = std::asin(r13);
        tvn = detail::tvn_adaptive_integral(h1, h2, h3, r23, a12, a13, epst) / (2.0 * M_PI);
        tvn += bivariate_normal_cdf(h2, h3, r23) * cephes::ndtr(h1);
    }

    return tvn < 0.0 ? 0.0 : (tvn > 1.0 ? 1.0 : tvn);
}

XSF_HOST_DEVICE inline float
trivariate_normal_cdf(float h, float k, float l, float r12, float r13, float r23, float epsi = 1e-7F) {
    return static_cast<float>(trivariate_normal_cdf(
        static_cast<double>(h), static_cast<double>(k), static_cast<double>(l), static_cast<double>(r12),
        static_cast<double>(r13), static_cast<double>(r23), static_cast<double>(epsi)
    ));
}

} // namespace xsf
