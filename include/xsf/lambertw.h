/* Translated from Cython into C++ by SciPy developers in 2023.
 * Original header with Copyright information appears below.
 */

/* Implementation of the Lambert W function [1]. Based on MPMath
 *  Implementation [2], and documentation [3].
 *
 * Copyright: Yosef Meller, 2009
 * Author email: mellerf@netvision.net.il
 *
 * Distributed under the same license as SciPy
 *
 *
 * References:
 * [1] On the Lambert W function, Adv. Comp. Math. 5 (1996) 329-359,
 *     available online: https://web.archive.org/web/20230123211413/https://cs.uwaterloo.ca/research/tr/1993/03/W.pdf
 * [2] mpmath source code,
 https://github.com/mpmath/mpmath/blob/c5939823669e1bcce151d89261b802fe0d8978b4/mpmath/functions/functions.py#L435-L461
 * [3]
 https://web.archive.org/web/20230504171447/https://mpmath.org/doc/current/functions/powers.html#lambert-w-function
 *

 * TODO: use a series expansion when extremely close to the branch point
 * at `-1/e` and make sure that the proper branch is chosen there.
 */

#pragma once

#include "config.h"
#include "error.h"
#include "evalpoly.h"
#include "cephes/polevl.h"

namespace xsf {
constexpr double EXPN1 = 0.36787944117144232159553; // exp(-1)
constexpr double OMEGA = 0.56714329040978387299997; // W(1, 0)

namespace detail {
    XSF_HOST_DEVICE inline std::complex<double> lambertw_branchpt(std::complex<double> z) {
        // Series for W(z, 0) around the branch point; see 4.22 in [1].
        double coeffs[] = {-1.0 / 3.0, 1.0, -1.0};
        std::complex<double> p = std::sqrt(2.0 * (M_E * z + 1.0));

        return cevalpoly(coeffs, 2, p);
    }

    XSF_HOST_DEVICE inline std::complex<double> lambertw_pade0(std::complex<double> z) {
        // (3, 2) Pade approximation for W(z, 0) around 0.
        double num[] = {12.85106382978723404255, 12.34042553191489361902, 1.0};
        double denom[] = {32.53191489361702127660, 14.34042553191489361702, 1.0};

        /* This only gets evaluated close to 0, so we don't need a more
         * careful algorithm that avoids overflow in the numerator for
         * large z. */
        return z * cevalpoly(num, 2, z) / cevalpoly(denom, 2, z);
    }

    XSF_HOST_DEVICE inline std::complex<double> lambertw_asy(std::complex<double> z, long k) {
        /* Compute the W function using the first two terms of the
         * asymptotic series. See 4.20 in [1].
         */
        std::complex<double> w = std::log(z) + 2.0 * M_PI * k * std::complex<double>(0, 1);
        return w - std::log(w);
    }

    XSF_HOST_DEVICE inline std::complex<double> lambertw_complex(std::complex<double> z, long k, double tol) {
        double absz;
        std::complex<double> w;
        std::complex<double> ew, wew, wewz, wn;

        if (std::isnan(z.real()) || std::isnan(z.imag())) {
            return z;
        }
        if (z.real() == std::numeric_limits<double>::infinity()) {
            return z + 2.0 * M_PI * k * std::complex<double>(0, 1);
        }
        if (z.real() == -std::numeric_limits<double>::infinity()) {
            return -z + (2.0 * M_PI * k + M_PI) * std::complex<double>(0, 1);
        }
        if (z == 0.0) {
            if (k == 0) {
                return z;
            }
            set_error("lambertw", SF_ERROR_SINGULAR, NULL);
            return -std::numeric_limits<double>::infinity();
        }
        if (z == 1.0 && k == 0) {
            // Split out this case because the asymptotic series blows up
            return OMEGA;
        }

        absz = std::abs(z);
        // Get an initial guess for Halley's method
        if (k == 0) {
            if (std::abs(z + EXPN1) < 0.3) {
                w = detail::lambertw_branchpt(z);
            } else if (-1.0 < z.real() && z.real() < 1.5 && std::abs(z.imag()) < 1.0 &&
                    -2.5 * std::abs(z.imag()) - 0.2 < z.real()) {
                /* Empirically determined decision boundary where the Pade
                * approximation is more accurate. */
                w = detail::lambertw_pade0(z);
            } else {
                w = detail::lambertw_asy(z, k);
            }
        } else if (k == -1) {
            if (absz <= EXPN1 && z.imag() == 0.0 && z.real() < 0.0) {
                w = std::log(-z.real());
            } else {
                w = detail::lambertw_asy(z, k);
            }
        } else {
            w = detail::lambertw_asy(z, k);
        }

        // Halley's method; see 5.9 in [1]
        if (w.real() >= 0) {
            // Rearrange the formula to avoid overflow in exp
            for (int i = 0; i < 100; i++) {
                ew = std::exp(-w);
                wewz = w - z * ew;
                wn = w - wewz / (w + 1.0 - (w + 2.0) * wewz / (2.0 * w + 2.0));
                if (std::abs(wn - w) <= tol * std::abs(wn)) {
                    return wn;
                }
                w = wn;
            }
        } else {
            for (int i = 0; i < 100; i++) {
                ew = std::exp(w);
                wew = w * ew;
                wewz = wew - z;
                wn = w - wewz / (wew + ew - (w + 2.0) * wewz / (2.0 * w + 2.0));
                if (std::abs(wn - w) <= tol * std::abs(wn)) {
                    return wn;
                }
                w = wn;
            }
        }

        set_error("lambertw", SF_ERROR_SLOW, "iteration failed to converge: %g + %gj", z.real(), z.imag());
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
} // namespace detail

XSF_HOST_DEVICE inline std::complex<double> lambertw(std::complex<double> z, long k, double tol) {
    return detail::lambertw_complex(z, k, tol);
}

XSF_HOST_DEVICE inline std::complex<float> lambertw(std::complex<float> z, long k, float tol) {
    return static_cast<std::complex<float>>(detail::lambertw_complex(static_cast<std::complex<double>>(z), k, static_cast<double>(tol))
    );
}

namespace fukushima {
    
    double P1[] = {-9.9999999999999988900E-1, -2.7399668668203659304,
        +2.6164207726990399347E-2, +6.3709168078949009170, +7.1013286517854026680,
        +2.9800826783006852573, +4.8819596813789865043E-1, +2.3753035787333611915E-2,
        +7.7365760093772431000E-5
    };
    double P2[] = {-9.9997801800578916749E-1, -7.0415751590483602272E-1,
        +2.1232260832802529071, +2.3896760702935718341, +7.7765311805029175244E-1,
        +8.9686698993644741433E-2, +3.3062485753746403559E-3, +2.5106760479132851033E-5
    };
    double P3[] = {-9.8967420337273506393E-1, +5.9587680606394382748E-1,
        +1.4225083018151943148, +4.4882889168323809798E-1, +4.4504943332390033511E-2,
        +1.5218794835419578554E-3, +1.6072263556502220023E-5, +3.3723373020306510843E-8
    };
    double P4[] = {-7.7316491997206225517E-1, +1.1391333504296703783,
        +4.3116117255217074492E-1, +3.5773078319037507449E-2, +9.6441640580559092740E-4,
        +8.9723854598675864757E-6, +2.5623503144117723217E-8, +1.4348813778416631453E-11
    };
    double P5[] = {+1.2007101671553688430E-1, +8.3352640829912822896E-1,
        +7.0142775916948337582E-2, +1.4846357985475124849E-3, +1.0478757366110155290E-5,
        +2.5715892987071038527E-8, +1.9384214479606474749E-11, +2.8447049039139409652E-15
    };
    double P6[] = {+1.7221104439937710112, +3.9919594286484275605E-1,
        +7.9885540140685028937E-3, +4.2889742253257920541E-5, +7.8146828180529864981E-8,
        +4.9819638764354682359E-11, +9.7650889714265294606E-15, +3.7052997281721724439E-19
    };
    double P7[] = {+3.7529314023434544256, +1.5491342690357806525E-1,
        +7.5663140675900784505E-4, +1.0271609235969979059E-6, +4.7853247675930066150E-10,
        +7.8328040770275474410E-14, +3.9433033758391036653E-18, +3.8232862205660283978E-23
    };
    double P8[] = {+6.0196542055606555577, +5.3496672841797864762E-2,
        +6.4340849275316501519E-5, +2.1969090100095967485E-8, +2.5927988937033061070E-12,
        +1.0779198161801527308E-16, +1.3780424091017898301E-21, +3.3768973150742552802E-27
    };
    double P9[] = {+8.4280268500989701597, +1.7155758546279713315E-2,
        +5.0836620669829321508E-6, +4.3354903691832581802E-10, +1.2841017145645583385E-14,
        +1.3419106769745885927E-19, +4.3101698455492225750E-25, +2.6422433422088187549E-31
    };
    double P10[] = {+1.0931063230472498189E+1, +5.2224234540245532982E-3,
        +3.7996105711810129682E-7, +8.0305793533410355824E-12, +5.9139785627090605866E-17,
        +1.5382020359533028724E-22, +1.2288944126268109432E-28, +1.8665089270660122398E-35
    };
    double P11[] = {+1.3502943080893871412E+1, +1.5284636506346264572E-3,
        +2.7156967358262346166E-8, +1.4110394051242161772E-13, +2.5605734311219728461E-19,
        +1.6421293724425337463E-25, +3.2324944691435843553E-32, +1.2054662641251783155E-39
    };
    double P12[] = {+1.6128076167439014775E+1, +4.3360385176467069131E-4,
        +1.8696403871820916466E-9, +2.3691795766901486045E-15, +1.0503191826963154893E-21,
        +1.6461927573606764263E-28, +7.9138276083474522931E-36, +7.1845890343701668760E-44
    };
    double P13[] = {+1.8796301105534486604E+1, +1.1989443339646469157E-4,
        +1.2463377528676863250E-10, +3.8219456858010368172E-17, +4.1055693930252083265E-24,
        +1.5595231456048464246E-31, +1.8157173553077986962E-39, +3.9807997764326166245E-48
    };
    double P14[] = {+2.1500582830667332906E+1, +3.2441943237735273768E-5,
        +8.0764963416837559148E-12, +5.9488445506122883523E-19, +1.5364106187215861531E-26,
        +1.4033231297002386995E-34, +3.9259872712305770430E-43, +2.0629086382257737517E-52
    };
    double P15[] = {+2.4235812532416977267E+1, +8.6161505995776802509E-6,
        +5.1033431561868273692E-13, +8.9642393665849638164E-21, +5.5254364181097420777E-29,
        +1.2045072724050605792E-37, +8.0372997176526840184E-47, +1.0049140812146492611E-56
    };
    double P16[] = {+2.6998134347987436511E+1, +2.2512257767572285866E-6,
        +3.1521230759866963941E-14, +1.3114035719790631541E-22, +1.9156784033962366146E-31,
        +9.8967003053444799163E-41, +1.5640423898448433548E-50, +4.6216193040664872606E-61
    };
    double P17[] = {+2.9784546702831970770E+1, +5.7971764392171329944E-7,
        +1.9069872792601950808E-15, +1.8668700870858763312E-24, +6.4200510953370940075E-34,
        +7.8076624650818968559E-44, +2.9029638696956315654E-54, +2.0141870458566179853E-65
    };
    double P18[] = {+7.4413499460126776143E-1, +4.1403243618005911160E-1,
        +2.6012564166773416170E-1, +2.1450457095960295520E-2, +5.1872377264705907577E-4,
        +4.3574693568319975996E-6, +1.2363066058921706716E-8, +9.0194147766309957537E-12
    };
    double P19[] = {-6.1514412812729761526E-1, +6.7979310133630936580E-1,
        +8.9685353704585808963E-2, +1.5644941483989379249E-3, +7.7349901878176351162E-6,
        +1.2891647546699435229E-8, +7.0890325988973812656E-12, +9.8419790334279711453E-16
    };
    double Q1[] = {+1, +5.0716108484174280050,
        +9.9868388183545283370, +9.6607551922078869080, +4.7943728991336119052,
        +1.1629703477704522300, +1.1849462500733755233E-1, +3.4326525132402226488E-3
    };
    double Q2[] = {+1, +3.0356026828085410884,
        +3.1434530151286777057, +1.3723156566592447275, +2.5844697415744211142E-1,
        +1.9551162251819044265E-2, +4.8775933244530123101E-4, +2.3165116841073152717E-6
    };
    double Q3[] = {+1, +1.6959402394626198052,
        +8.0968573415500900896E-1, +1.4002034999817021955E-1, +9.3571878493790164480E-3,
        +2.3251487593389773464E-4, +1.8060170751502988645E-6, +2.5750667337015924224E-9
    };
    double Q4[] = {+1, +7.4657287456514418083E-1,
        +1.2629777033419350576E-1, +6.9741512959563184881E-3, +1.4089339244355354892E-4,
        +1.0257432883152943078E-6, +2.2902687190119230940E-9, +9.2794231013264501664E-13
    };
    double Q5[] = {+1, +2.5396738845619126630E-1,
        +1.2839238907330317393E-2, +2.0275375632510997371E-4, +1.1482956073449141384E-6,
        +2.3188370605674263647E-9, +1.4271994165742563419E-12, +1.5884836942394796961E-16
    };
    double Q6[] = {+1, +7.4007438118020543008E-2,
        +1.0333501506697740545E-3, +4.4360858035727508506E-6, +6.7822912316371041570E-9,
        +3.6834356707639492021E-12, +6.0836159560266041168E-16, +1.8149869335981225316E-20
    };
    double Q7[] = {+1, +2.0112985338854443555E-2,
        +7.4712286154830141768E-5, +8.4800598003693837469E-8, +3.4182424130376911762E-11,
        +4.8866259139690957899E-15, +2.1223373626834634178E-19, +1.6642985671260582515E-24
    };
    double Q8[] = {+1, +5.2809683704233371675E-3,
        +5.1020501219389558082E-6, +1.5018312292270832103E-9, +1.5677706636413188379E-13,
        +5.7992041238911878361E-18, +6.5133170770320780259E-23, +1.3205080139213406071E-28
    };
    double Q9[] = {+1, +1.3572006754595300315E-3,
        +3.3535243481426203694E-7, +2.5206969246421264128E-11, +6.7136226273060530496E-16,
        +6.3324226680854686574E-21, +1.8128167400013774194E-26, +9.3662030058136796889E-33
    };
    double Q10[] = {+1, +3.4328702551197577797E-4,
        +2.1395351518538844476E-8, +4.0524170186631594159E-13, +2.7181424315335710420E-18,
        +6.4538986638355490894E-24, +4.6494613785888987942E-30, +6.0442024367299387616E-37
    };
    double Q11[] = {+1, +8.5701512879089462255E-5,
        +1.3311244435752691563E-9, +6.2788924440385347269E-15, +1.0483788152252204824E-20,
        +6.1943499966249160886E-27, +1.1101567860340917294E-33, +3.5897381128308962590E-41
    };
    double Q12[] = {+1, +2.1154255263102938752E-5,
        +8.1006115442323280538E-11, +9.4155986022169905738E-17, +3.8725127902295302254E-23,
        +5.6344651115570565066E-30, +2.4860951084210029191E-37, +1.9788304737427787405E-45
    };
    double Q13[] = {+1, +5.1691031988359922329E-6,
        +4.8325571823313711932E-12, +1.3707888746916928107E-18, +1.3754560850024480337E-25,
        +4.8811882975661805184E-33, +5.2518641828170201894E-41, +1.0192119593134756440E-49
    };
    double Q14[] = {+1, +1.2515317642433850197E-6,
        +2.8310314214817074806E-13, +1.9423666416123637998E-20, +4.7128616004157359714E-28,
        +4.0433347391839945960E-36, +1.0515141443831187271E-44, +4.9316490935436927307E-54
    };
    double Q15[] = {+1, +3.0046761844749477987E-7,
        +1.6309104270855463223E-14, +2.6842271030298931329E-22, +1.5619672632458881195E-30,
        +3.2131689030397984274E-39, +2.0032396245307684134E-48, +2.2520274554676331938E-58
    };
    double Q16[] = {+1, +7.1572676370907573898E-8,
        +9.2500506091115760826E-16, +3.6239819582787573031E-24, +5.0187712493800424118E-33,
        +2.4565861988218069039E-42, +3.6435658433991660284E-52, +9.7432490640155346004E-63
    };
    double Q17[] = {+1, +1.6924463180469706372E-8,
        +5.1703934311254540111E-17, +4.7871532721560069095E-26, +1.5664405832545149368E-35,
        +1.8113137982381331398E-45, +6.3454150289495419529E-56, +4.0072964025244397967E-67
    };
    double Q18[] = {+1, +3.3487811067467010907E-1,
        +2.3756834394570626395E-2, +5.4225633008907735160E-4, +4.4378980052579623037E-6,
        +1.2436585497668099330E-8, +9.0225825867631852215E-12, -4.2057836270109716654E-19
    };
    double Q19[] = {+1, +9.7300263710401439315E-2,
        +1.6103672748442058651E-3, +7.8247741003077000012E-6, +1.2949261308971345209E-8,
        +7.0986911219342827130E-12, +9.8426285042227044979E-16, -1.5960147252606055352E-24
    };

    XSF_HOST_DEVICE inline double rational_function(double x, double A[], int nA, double B[], int nB) {
        double numerator = cephes::polevl(x, A, nA);
        double denominator = cephes::polevl(x, B, nB);
        return numerator / denominator;
    }

    constexpr double z0 = -0.36787944117144232160;

    XSF_HOST_DEVICE inline double lambertw(double z, long k) {
        if(z < -0.36787944117144232160){ return std::numeric_limits<double>::quiet_NaN(); }
        if(z < +2.1820144653320312500){ return rational_function(std::sqrt(z - z0), P1, 9, Q1, 8); }
        if(z < +4.3246045021497925573E+1){ return rational_function(std::sqrt(z - z0), P2, 8, Q2, 8); }
        if(z < +5.9808565427761132714E+2){ return rational_function(std::sqrt(z - z0), P3, 8, Q3, 8); }
        if(z < +8.0491241056345904686E+3){ return rational_function(std::sqrt(z - z0), P4, 8, Q4, 8); }
        if(z < +1.1112458624177664276E+5){ return rational_function(std::sqrt(z - z0), P5, 8, Q5, 8); }
        if(z < +1.5870426133287885398E+6){ return rational_function(std::sqrt(z - z0), P6, 8, Q6, 8); }
        if(z < +2.3414708033996018338E+7){ return rational_function(std::sqrt(z - z0), P7, 8, Q7, 8); }
        if(z < +3.5576474271222021108E+8){ return rational_function(std::sqrt(z - z0), P8, 8, Q8, 8); }
        if(z < +5.5501716292484833443E+9){ return rational_function(std::sqrt(z - z0), P9, 8, Q9, 8); }
        if(z < +8.8674704839289895890E+10){ return rational_function(std::sqrt(z - z0), P10, 8, Q10, 8); }
        if(z < +1.4477791865269224022E+12){ return rational_function(std::sqrt(z - z0), P11, 8, Q11, 8); }
        if(z < +2.4111458632511484051E+13){ return rational_function(std::sqrt(z - z0), P12, 8, Q12, 8); }
        if(z < +4.0897036442600808776E+14){ return rational_function(std::sqrt(z - z0), P13, 8, Q13, 8); }
        if(z < +7.0555901476789968723E+15){ return rational_function(std::sqrt(z - z0), P14, 8, Q14, 8); }
        if(z < +1.2366607557976727250E+17){ return rational_function(std::sqrt(z - z0), P15, 8, Q15, 8); }
        if(z < +2.1999373487930999771E+18){ return rational_function(std::sqrt(z - z0), P16, 8, Q16, 8); }
        if(z < +3.9685392198344016155E+19){ return rational_function(std::sqrt(z - z0), P17, 8, Q17, 8); }
        if(z < +1.4127075145274652069E+104){ return rational_function(std::log(z), P18, 8, Q18, 8); }
        if(z < +2.8134195736211426913E+618){ return rational_function(std::log(z), P19, 8, Q19, 8); }
        if(std::isnan(z)){ return std::numeric_limits<double>::quiet_NaN(); }
        if(z == std::numeric_limits<double>::infinity()){ return std::numeric_limits<double>::infinity(); }
        if(z == -std::numeric_limits<double>::infinity()){ return -std::numeric_limits<double>::infinity(); }
    }

} // namespace fukushima

XSF_HOST_DEVICE inline double lambertw(double z, long k, double tol) {
    return fukushima::lambertw(z, k);
}

} // namespace xsf
