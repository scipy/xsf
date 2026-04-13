/* 
 * Calculation of the trigamma function for real and complex inputs
 * Author: Lorenzo Peri
 */

#pragma once

#include "cephes/psi.h"
#include "cephes/zeta.h"
#include "config.h"
#include "error.h"
#include "trig.h"

namespace xsf {
namespace detail {

    /* There exist no real roots of the trigamma function.
     * All the roots of the trigamma function lay in the negative real semi-plane (Re[z] < 0)
     * and come in an infinite number of complex pairs (z_n, \bar{z}_n).
     * This is the location of the two roots closer to the origin.
     */

    constexpr double firstroot_re = -0.41213454795193723424;
    constexpr double firstroot_im = +0.59781194232059657221;

    constexpr double secondroot_re = -1.4455692007504016647;
    constexpr double secondroot_im = +0.69926087106525701723;

    /* NOTE: 
     * For better accuracy, one may be able to use the fact that 
     *    trigamma(z) + trigamma(1-z) = pi**2 / sin(pi * z)**2
     * So for large negative Re[z], the roots must asymptotically approach
     * real part -n + 1/2, with slowly-growing imaginary part. Hence, for values around
     * z = -n + 1/2 +- 0.75j we could use the Taylor expansion via the polygamma, 
     * if more precision is needed for root-finding. 
     * However, testing suggests this is not the case, so this has not been implemented.
     */

    XSF_HOST_DEVICE inline std::complex<double> 
    trigamma_root_series(std::complex<double> const z, int const root_index, double const imsign) {
        
        /********************************** Series expansion constants **********************************/
        /* Because the trigamma does not have roots in the real axis, we cannot piggy-back on cephes::zeta.
         * However, psi(n+1, z)/(n)! can still be evaluated via mpmath and stored as a constant array
        */
    
        static constexpr std::array<double, 2> roots_location_re = {firstroot_re, secondroot_re};
        
        static constexpr std::array<double, 2> roots_location_im = {firstroot_im, secondroot_im};

        static constexpr std::array<std::array<double, 100>, 2> taylor_coeffs_re = {{
        { // First root
            9.1080172479838895865e-17, 
            -2.9833604957805599511e+0, -14.018950391053259352e+0, 9.0244567682892409266e+0, 30.874956766351548509e+0,
            34.541752961058271865e+0, 39.009566174081896861e+0, -132.61045192741005394e+0, -218.62357592235318293e+0,
            -64.411882496106031226e+0, 203.10382610710644258e+0, 843.63339672989093287e+0, 662.46913770761409523e+0,
            -772.94967564645764924e+0, -2195.3718826638200881e+0, -2930.2957335752762447e+0, 645.15496382139008347e+0,
            7427.2379542621074506e+0, 9439.1355798311033141e+0, 2265.3656112362482418e+0, -17976.059265797157423e+0,
            -34557.165315543585166e+0, -15488.046929912889027e+0, 40467.071329236867314e+0, 102237.07007987287943e+0,
            84308.324421938203159e+0, -81451.135771210465464e+0, -290705.554540852434e+0, -314668.7054179974366e+0,
            78288.649124173025484e+0, 785651.99318701238371e+0, 1078335.3382942241151e+0, 173941.59359342456446e+0,
            -1883168.9550851269159e+0, -3430120.4253031355329e+0, -1663440.5002393615432e+0, 4158107.3425906528719e+0,
            9981119.5116241183132e+0, 7815988.8380806706846e+0, -7551943.8107375940308e+0, -27587495.006474528462e+0,
            -29076752.908828314394e+0, 8077666.9182099932805e+0, 71044589.308923915029e+0, 97353898.125684320927e+0,
            14547456.322059284896e+0, -169035489.75971624255e+0, -299438887.83402949572e+0, -143212831.15478402376e+0,
            363534700.00453239679e+0, 861125488.17706620693e+0, 657213911.34497916698e+0, -650972657.23167860508e+0,
            -2332338481.4032406807e+0, -2429691855.3164582253e+0, 716433626.61413443089e+0, 5917764178.0552072525e+0,
            8008883489.7560710907e+0, 1112483090.8655650616e+0, -13937926524.107784271e+0, -24338541099.401115417e+0,
            -11362293924.424409866e+0, 29639987087.049156189e+0, 69304135290.439422607e+0, 52046135316.926963806e+0,
            -52924464353.45640564e+0, -185777990501.79312134e+0, -191320182474.67971802e+0, 59363105651.442115784e+0,
            467673151063.63690186e+0, 625806128716.41210938e+0, 79708654667.238372803e+0, -1094025457975.3127441e+0,
            -1889131923265.0839844e+0, -861612477677.46520996e+0, 2315097770655.0478516e+0, 5344755178877.3457031e+0,
            3957704192960.6518555e+0, -4133617154330.4306641e+0, -14243561160457.984375e+0, -14507740357232.351563e+0,
            4730156227368.0966797e+0, 35680709502628.15625e+0, 47264902566753.390625e+0, 5466706316509.2929688e+0,
            -83134951089908.90625e+0, -142094372390481.15625e+0, -63354475700337.945313e+0, 175532075173260.25e+0,
            400403810133614.4375e+0, 292525586156113.625e+0, -313960697597505.6875e+0, -1063201394869874.0e+0,
            -1071636435602862.75e+0, 366735533339498.3125e+0, 2655187226255846.5e+0, 3484089660964081.0e+0,
            362143126064491.1875e+0, -6172561211494055.0e+0, -10448437646105654.0e+0
        },
        { // Second root
            5.4524533247549979198e-17, 
            -1.1129516609162204421e+0, -8.4123458270441080487e+0, 4.4344884179193835294e+0, 20.732204861034905008e+0,
            3.0429402875860867006e+0, 4.8856876287330912589e+0, -34.394886945411350609e+0, -73.718725710581068711e+0,
            34.397533695264762343e+0, 82.019204983273539256e+0, 87.311394512730743145e+0, 79.813123929215834096e+0,
            -229.88961911997415655e+0, -307.25357348392651602e+0, 17.768513154646232266e+0, 233.9103567437451261e+0,
            644.57447585988916217e+0, 400.31597610602159421e+0, -838.99527370542307381e+0, -1142.667274602767975e+0,
            -672.50854221176700776e+0, 699.33343286998206167e+0, 2952.2729899380997267e+0, 1804.4156337893034561e+0,
            -1939.0222428852139274e+0, -4390.9776154142909945e+0, -4765.9689718843501396e+0, 1782.0073007910705201e+0,
            10472.703050074791463e+0, 8540.7804724000234273e+0, -1795.7734006212097029e+0, -16338.586309481521312e+0,
            -22056.889110903408437e+0, 1042.6133896764329165e+0, 31820.13958963152254e+0, 38860.710514671307465e+0,
            10325.776471356321053e+0, -53936.772552904549229e+0, -85802.932425871578744e+0, -23282.380866077233804e+0,
            86119.543405023272499e+0, 158925.90507253524265e+0, 84995.099389624680043e+0, -148968.3693368526001e+0,
            -305690.65955470677e+0, -182464.18278398935217e+0, 200570.45908379976754e+0, 578305.62896160327364e+0,
            438149.73490933538415e+0, -311083.99629818135872e+0, -1022984.1769126829458e+0, -934179.80446570529602e+0,
            320098.97928918036632e+0, 1889408.096175869694e+0, 1919098.4924058187753e+0, -260137.09258919744752e+0,
            -3194477.7004772522487e+0, -3949612.7951777679846e+0, -268240.58420638204552e+0, 5568022.6676272731274e+0, 
            7594112.9608210055158e+0, 1866021.784239478875e+0, -9080501.0774496924132e+0, -14858851.839382385835e+0,
            -5612890.7087618308142e+0, 14571131.499402185902e+0, 27616979.832255661488e+0, 14860656.240491846576e+0,
            -22371271.646269556135e+0, -51332035.234632596374e+0, -33734206.887385599315e+0, 31637528.883603923023e+0,
            92778118.069112360477e+0, 74537842.063162609935e+0, -41599123.161676205695e+0, -164575435.80708390474e+0,
            -154751889.32015013695e+0, 41749336.101391732693e+0, 287615800.32811558247e+0, 313345749.76099699736e+0,
            -17522700.13627949357e+0, -487787122.48196721077e+0, -617073974.94401979446e+0, -81001990.910296738148e+0,
            814449173.36699914932e+0, 1183835029.963737011e+0, 344342132.1734226346e+0, -1310270797.1101074219e+0,
            -2235174639.0448431969e+0, -972870005.92697286606e+0, 2044440891.2078654766e+0, 4120640100.7729520798e+0,
            2371630467.9593176842e+0, -3026173792.9559259415e+0, -7486827096.7126846313e+0, -5297885493.9058942795e+0,
            4160137914.0421552658e+0, 13321713752.551465988e+0, 11278342754.757139206e+0
        }
        }};

        static constexpr std::array<std::array<double, 100>, 2> taylor_coeffs_im= {{
        { // First root
            1.9689476472890448436e-16, 
            4.044337981596503262e+0, -7.3941432193042802012e+0, -26.4623069876907131e+0, -1.5657568234258929518e+0,
            12.808397661710225535e+0, 88.027833394290681213e+0, 123.62940643728001078e+0, -106.32109009482918793e+0,
            -272.45867815419285307e+0, -410.27930694666247291e+0, -86.782597941426899979e+0, 1102.1387232076265263e+0,
            1466.0506454643184497e+0, 545.1882812919824346e+0, -2217.9064686095880461e+0, -5754.8727812528886716e+0,
            -3283.4424631773640613e+0, 5543.4392134415256805e+0, 15833.635370604708442e+0, 16682.880006384319131e+0,
            -9435.5858181943276577e+0, -47167.297127159035881e+0, -56722.714163294331229e+0, -854.28115876288268282e+0,
            124033.47128391949809e+0, 196500.84514077924541e+0, 63124.484428748924984e+0, -284327.58732096874155e+0,
            -608727.64227164373733e+0, -392247.35741802095436e+0, 604454.0479207362514e+0, 1734042.636362915393e+0,
            1621968.5883482831996e+0, -909822.08165691560134e+0, -4736249.7619706066325e+0, -5710757.7356921844184e+0,
            207452.20227611850714e+0, 11818016.939825890586e+0, 18576645.04214990139e+0, 6031726.6782435318455e+0,
            -27102249.135446421802e+0, -55521819.735950939357e+0, -35071050.387607783079e+0, 54519628.983876109123e+0,
            156158894.62619554996e+0, 141738560.80960172415e+0, -82503555.349761158228e+0, -412818986.95092844963e+0,
            -494446402.06373292208e+0, 25868788.41499729827e+0, 1015950908.8245180845e+0, 1569817888.5925683975e+0,
            497068007.28298789263e+0, -2297081599.2161045074e+0, -4632770027.5792331696e+0, -2861134230.4570889473e+0,
            4554388755.4394798279e+0, 12857647166.110887527e+0, 11502606847.624053955e+0, -6921551638.6264343262e+0,
            -33553555348.813869476e+0, -39739817902.058311462e+0, 2586764571.5984015465e+0, 81853762473.440658569e+0,
            124901288496.36500549e+0, 38189522192.376235962e+0, -183499975705.06417847e+0, -365651105247.84295654e+0,
            -221707664548.28396606e+0, 362081755311.74279785e+0, 1006612243979.449707e+0, 889379762919.08422852e+0,
            -552674484944.32897949e+0, -2608847120879.1108398e+0, -3056520000009.1757813e+0, 238946410783.76174927e+0,
            6328034955998.734375e+0, 9555065868385.1386719e+0, 2817455991885.3300781e+0, -14122437331019.056641e+0,
            -27828358863467.332031e+0, -16589801973106.900391e+0, 27819982471996.28125e+0, 76235974161174.5e+0,
            66565105419682.960938e+0, -42721272861724.007813e+0, -196745977498914.03125e+0, -228195034186489.65625e+0,
            20884379691198.875e+0, 475554441415672.3125e+0, 711150412786520.875e+0, 202086534694355.9375e+0,
            -1058822981010458.0e+0, -2064494960700603.5e+0, -1210578706349626.5e+0, 2085822266866922.5e+0,
            5638495776870608.0e+0, 4867352096330547.0e+0, -3226346588093619.5e+0,
        },
        { // Second root
            -2.0886187184514818594e-16,
            2.6352910681361323419e+0, -2.7465273316761171785e+0, -16.081603936774637731e+0, 3.6349150477985818775e+0,
            15.326822302901804917e+0, 16.963613417604033629e+0, 38.030023581648258357e+0, -45.080272523330833678e+0,
            -94.683126434869024024e+0, -9.8547858644971260134e+0, 23.13494927816051927e+0, 175.86022148586815206e+0,
            204.21597173963738214e+0, -191.75729307902071241e+0, -334.08125960690273359e+0, -287.7192524549465702e+0,
            -18.440040361671922398e+0, 894.44902390096001454e+0, 825.34530871214099079e+0, -320.73324130858628678e+0,
            -1160.3663328222710334e+0, -1920.3364016438486033e+0, -324.32309629687068764e+0, 3144.1574371753281412e+0,
            3371.5495145367981422e+0, 851.25003006116332926e+0, -4081.6620591933688047e+0, -8578.2344898673782154e+0,
            -2668.3694373718035422e+0, 8579.661814979142946e+0, 14009.367669937275423e+0, 9337.3016036521203205e+0,
            -12620.991744398463197e+0, -31502.710372905752592e+0, -17620.59367900545476e+0, 18235.149970555630716e+0,
            54404.217666152544552e+0, 48932.207386342364771e+0, -28758.606484943069518e+0, -103989.97488247773435e+0,
            -93196.077915846719407e+0, 21679.049627542692178e+0, 188424.71855478230282e+0, 209095.51932349728304e+0,
            -21708.310137294363813e+0, -318917.72407036350342e+0, -412294.83355757704703e+0, -65133.869098773830046e+0,
            573262.05199860793073e+0, 811407.75078429805581e+0, 223905.47458710509818e+0, -898392.3956043071812e+0,
            -1594942.9401338619646e+0, -692055.71907006460242e+0, 1502419.6342246450949e+0, 2936510.239670007024e+0,
            1722814.7067913736682e+0, -2194149.2624233397655e+0, -5568379.2814786368981e+0, -3861876.4619128271006e+0,
            3140539.447914307937e+0, 9908628.1384462881833e+0, 8556352.4623205848038e+0, -3885309.5242206109688e+0,
            -17843844.457143317908e+0, -17440540.353781145066e+0, 3359989.4226588201709e+0, 30885963.191879630089e+0,
            35667011.386512152851e+0, 106318.61343497049529e+0, -52540829.694197170436e+0, -69449093.001229658723e+0,
            -13246966.086865155026e+0, 87453795.719745591283e+0, 133842334.35044507682e+0, 44941680.012016184628e+0,
            -139416689.22975465655e+0, -251742409.82044029236e+0, -121905667.62976655364e+0, 216799719.98306891322e+0,
            463634394.69740325212e+0, 287096515.37800866365e+0, -313586151.84701472521e+0, -842334206.9316624403e+0,
            -633756729.02875030041e+0, 421789259.77867001295e+0, 1492942912.7488236427e+0, 1335806458.6442694664e+0,
            -474527237.98944181204e+0, -2610230524.3749175072e+0, -2706852540.8357629776e+0, 324012012.95333015919e+0,
            4450062790.8111820221e+0, 5356330811.5707168579e+0, 390413362.75755184889e+0, -7433131911.5818634033e+0,
            -10309588508.552108765e+0, -2484592309.8536400795e+0, 12070922242.493497849e+0
        }
        }};

        /************************************************ Series evaluation ************************************************/

        std::complex<double> const root(roots_location_re[root_index], imsign * roots_location_im[root_index]);
        std::complex<double> const zr = z - root;

        std::complex<double> res(taylor_coeffs_re[root_index][0], imsign * taylor_coeffs_im[root_index][0]);
        std::complex<double> coeff = 1.;
        std::complex<double> term;
        for (int n = 1; n < 100; n++) {
            coeff *= zr;
            term = coeff * std::complex<double>(taylor_coeffs_re[root_index][n], imsign * taylor_coeffs_im[root_index][n]);
            res += term;
            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(res)) {
                break;
            }
        }
        return res;
    }

    XSF_HOST_DEVICE inline std::complex<double>
    trigamma_forward_recurrence(std::complex<double> const z, std::complex<double> const psiz, int const n) {
        /* Compute trigamma(z + n) using trigamma(z) using the recurrence relation
         * trigamma(z + 1) = trigamma(z) - 1/z^2.
         * See http://dlmf.nist.gov/5.15.E5 */

        std::complex<double> zk, res = psiz;

        for (int k = 0; k < n; k++) {
            zk = (z + static_cast<double>(k));
            res -= 1./(zk * zk);
        }
        return res;
    }

    XSF_HOST_DEVICE inline std::complex<double>
    trigamma_backward_recurrence(std::complex<double> const z, std::complex<double> const psiz, int const n) {
        /* Compute trigamma(z + n) using trigamma(z) using the recurrence relation
         * trigamma(z + 1) = trigamma(z) - 1/z^2.
         * See http://dlmf.nist.gov/5.15.E5 */

        std::complex<double> zk, res = psiz;

        for (int k=1;k<n+1;k++) {
            zk = (z - static_cast<double>(k));
            res += 1./(zk * zk);
        }
        return res;
    }

    XSF_HOST_DEVICE inline std::complex<double> trigamma_asymptotic_series(std::complex<double> z) {
        /* Evaluate digamma using an asymptotic series. See
         * http://dlmf.nist.gov/5.15.E8
         * Higher order converge slower, so we need slightly more Bernoulli 
         * numbers than the digamma (Bernoulli numbers calculated via sympy).
         */
        static constexpr std::array<double, 24> bernoulli2k = {
            0.16666666666666665741e+0, -0.033333333333333332871e+0, 0.023809523809523808202e+0, 
            -0.033333333333333332871e+0, 0.075757575757575759678e+0, -0.25311355311355310249e+0,
            1.1666666666666667407e+0, -7.0921568627450977118e+0, 54.971177944862155584e+0, 
            -529.12424242424242493e+0, 6192.1231884057970092e+0, -86580.253113553117146e+0, 
            1425517.1666666667443e+0, -27298231.067816093564e+0, 601580873.90064239502e+0, 
            -15116315767.092157364e+0, 429614643061.16668701e+0, -13711655205088.332031e+0, 
            488332318973593.1875e+0, -19296579341940068.0e+0, 841693047573682560.0e+0, 
            -40338071854059454464.0e+0, 2.1150748638081992622e+21, -1.2086626522296526202e+23,
        };
        
        std::complex<double> const rz = 1.0 / z;
        std::complex<double> const rzz = rz*rz;
        std::complex<double> zfac = rz;

        std::complex<double> term;
        std::complex<double> res;

        if (!(std::isfinite(z.real()) && std::isfinite(z.imag()))) {
            /* Check for infinity (or nan) and return early.
             * Result of division by complex infinity is implementation dependent.
             * and has been observed to vary between C++ stdlib and CUDA stdlib.
             */
            return 0.; // Asymptotically psi(1, z) ~ 1/z, so just return zero?
        }

        res = rz + 0.5*rzz;

        for (int k = 1; k < bernoulli2k.size()+1; k++) {
            zfac *= rzz;
            term = bernoulli2k[k-1]*zfac;
            res += term;
            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(res)) {
                break;
            }
        }
        return res;
    }

} // namespace detail

XSF_HOST_DEVICE inline double trigamma(double z) {
    /* To compute the trigamma for real inputs we can leverage Cehpes via the identity
     *    psi(n, x) = (-1)**(n+1) * n! * zeta(n+1, x)
     * where n is a positive integer and zeta is the Hurwitz zeta function
     */
    return cephes::zeta(2., z);
}

XSF_HOST_DEVICE inline float trigamma(float z) { return static_cast<float>(trigamma(static_cast<double>(z))); }

XSF_HOST_DEVICE inline std::complex<double> trigamma(std::complex<double> z) {
    /*
     * Compute the trigamma function for complex arguments. The strategy is:
     *
     * - Around the two zeros closest to the origin
     * use a Taylor series with precomputed coefficient.
     * - If close to the origin, use a recurrence relation to step away
     * from the origin.
     * - If close to the negative real axis, use the reflection formula
     * to move to the right halfplane.
     * - If |z| is large (> 16), use the asymptotic series.
     * - If |z| is small, use a recurrence relation to make |z| large
     * enough to use the asymptotic series.
     */

    int n;
    double absz = std::abs(z);
    std::complex<double> res = 0.;
    /* The reflection formula flips the sign for odd polygammas.
     * We must keep track of that. */
    std::complex<double> sign = 1.;
    /* Use the asymptotic series for z away from the negative real axis
     * with abs(z) > smallabsz. */
    int const smallabsz = 16;

    if (z.real() <= 0.0 && std::ceil(z.real()) == z) {
        // Poles
        set_error("trigamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    if(std::abs(z.imag()) < std::numeric_limits<double>::epsilon()){
        // use the zeta, because it evaluates to better precision
        return static_cast<std::complex<double>>(cephes::zeta(2., z.real()));
    } 
    // Near the first two roots, use the asymptotic series
    if (std::abs(z - std::complex<double>(detail::firstroot_re, +detail::firstroot_im)) < 0.1){
        return detail::trigamma_root_series(z, 0, +1.);
    } else if (std::abs(z - std::complex<double>(detail::firstroot_re, -detail::firstroot_im)) < 0.1){
        return detail::trigamma_root_series(z, 0, -1.);
    } else if (std::abs(z - std::complex<double>(detail::secondroot_re, +detail::secondroot_im)) < 0.2){
        return detail::trigamma_root_series(z, 1, +1.);
    } else if (std::abs(z - std::complex<double>(detail::secondroot_re, -detail::secondroot_im)) < 0.2){
        return detail::trigamma_root_series(z, 1, -1.);
    }

    if (z.real() < 0 && std::abs(z.imag()) < smallabsz) {
        /* Reflection formula.
         * trigamma(z) + trigamma(1-z) = pi**2 / sin(pi * z)**2
         * https://dlmf.nist.gov/5.15#E5 
         * NOTE: sign change */
        std::complex<double> const ssp = sinpi(z);
        res += (M_PI * M_PI) / (ssp * ssp) ;
        sign *= -1;
        z = 1. - z;
        absz = std::abs(z);
    }

    if(absz < 0.5){
        // Use one step of the recurrence relation to step away from the pole
        res += sign/(z*z);
        z += 1;
        absz = std::abs(z);
    }

    // At this point we are sure to be far away from any root
    if (absz > smallabsz) {
        res += sign * detail::trigamma_asymptotic_series(z);
    } else if (z.real() >= 0.0) {
        double n = std::trunc(smallabsz - absz) + 1;
        std::complex<double> init = detail::trigamma_asymptotic_series(z + n);
        res += sign * detail::trigamma_backward_recurrence(z + n, init, n);
    } else {
        double n = std::trunc(smallabsz - absz) - 1;
        std::complex<double> init = detail::trigamma_asymptotic_series(z - n);
        res += sign * detail::trigamma_forward_recurrence(z - n, init, n);
    }

    return res;
}

XSF_HOST_DEVICE inline std::complex<float> trigamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(trigamma(static_cast<std::complex<double>>(z)));
}

} // namespace xsf
