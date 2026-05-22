#include "../testing_utils.h"
#include <cmath>
#include <limits>
#include <xsf/cdflib.h>

TEST_CASE("gdtria basic", "[gdtria][xsf_tests]") {

    SECTION("gdtria double precision invalid and edge case inputs") {
        // test cases with invalid and edge case  inputs that should return NaN
        using test_case = std::tuple<double, double, double, double>;
        auto [p, b, x, expected] = GENERATE(
            test_case{
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            }, // all NaN
            test_case{
                0.5, std::numeric_limits<double>::quiet_NaN(), 0.2, std::numeric_limits<double>::quiet_NaN()
            }, // NaN b
            test_case{
                0.5, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            }, // NaN x
            test_case{
                std::numeric_limits<double>::quiet_NaN(), 1.0, 0.2, std::numeric_limits<double>::quiet_NaN()
            },                                                                   // NaN p
            test_case{0.5, 0.0, 0.2, std::numeric_limits<double>::quiet_NaN()},  //  b == 0
            test_case{0.5, -1.0, 0.2, std::numeric_limits<double>::quiet_NaN()}, // b < 0
            test_case{0.0, 1.0, 0.2, 0.0},                                       // p == 0
            test_case{1.0, 1.0, 0.2, std::numeric_limits<double>::infinity()},   // p == 1
            test_case{-0.1, 1.0, 0.2, std::numeric_limits<double>::quiet_NaN()}, // p < 0
            test_case{1.1, 1.0, 0.2, std::numeric_limits<double>::quiet_NaN()},  // p > 1
            test_case{0.5, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()},  // x == 0
            test_case{0.5, 2.0, -1.0, std::numeric_limits<double>::quiet_NaN()}  // x < 0
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::gdtria(p, b, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, b, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("gdtria double precision scipy reference values") {
        // Reference values computed with scipy.special.gdtria using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // p_s = rng.uniform(0.01, 0.99, 30)
        // b_s = rng.uniform(0.1, 10.0, 30)
        // x_s = rng.uniform(0.1, 10.0, 30)

        // for p, b, x in zip(p_s, b_s, x_s):
        //     output = special.gdtria(p, b, x)
        //     print(f"test_case{{{p}, {b}, {x}, {output}}}")
        using test_case = std::tuple<double, double, double, double>;
        auto [p, b, x, expected] = GENERATE(
            test_case{0.23278930201782627, 8.561944853302613, 2.6268003028130527, 2.405084269521892},
            test_case{0.32042317291555783, 6.056050292767759, 3.045896154665676, 1.538564193169049},
            test_case{0.7914181481860795, 9.326684775246237, 2.8627644861562898, 4.067752217331768},
            test_case{0.6727295773359551, 7.275335474810999, 2.679734203663846, 3.0550427850207247},
            test_case{0.3932873595898708, 8.619458042193594, 4.879339768713225, 1.5444305277470047},
            test_case{0.33615764930905684, 9.30044423559563, 2.198592459879955, 3.531925470472},
            test_case{0.596342578515446, 5.507241489915295, 5.006742907063102, 1.1483671939698457},
            test_case{0.19299950189163909, 9.382962291800794, 2.537987125724302, 2.6379984475056584},
            test_case{0.6693009231343289, 5.00038060678036, 8.400978259422754, 0.6757339534261768},
            test_case{0.9329668079645385, 2.8103545066508766, 1.8832928419408475, 2.9806603319044016},
            test_case{0.2532808003369796, 4.5726092040001305, 8.63534728594144, 0.3504995274097965},
            test_case{0.9399035287966518, 6.68388534165535, 1.865164503967356, 5.947696666265312},
            test_case{0.663892704038365, 3.375820211623841, 7.530260186178716, 0.5101045393295068},
            test_case{0.10397997688222983, 9.044194667401566, 6.150091997922595, 0.8967495954687005},
            test_case{0.44300287284445655, 2.64503433523769, 2.1706348457932125, 0.9726387404457946},
            test_case{0.8787503209409673, 3.4643005423421664, 7.622736969127553, 0.7429765390556826},
            test_case{0.6935044298843817, 2.662648646564981, 2.5676796383956337, 1.2416836601809187},
            test_case{0.32994340678870987, 3.6189201514484317, 0.947160146666924, 2.705696430644921},
            test_case{0.7292496000634652, 0.1497211037996047, 6.2187615509491, 0.013187304666629089},
            test_case{0.2257322564345765, 6.323184986586819, 5.415986477220121, 0.8016101315277192},
            test_case{0.08996267815136395, 2.895588803508671, 6.38181444103123, 0.15503259007767223},
            test_case{0.16669768905354657, 0.7740681259306629, 1.8263036760447438, 0.05154580955696343},
            test_case{0.3432981812556112, 6.206606874838167, 2.5568284495788793, 1.9385221690218037},
            test_case{0.46588929062800993, 1.845630570783914, 6.879747547930052, 0.20695096221503495},
            test_case{0.2710926077249555, 3.1134450334739374, 0.900629297883984, 2.108079040310248},
            test_case{0.8094608753563108, 4.464779427673569, 8.763228647485649, 0.7038889193249693},
            test_case{0.1994285015037046, 1.5870031765207382, 4.344074377245919, 0.12775141436455284},
            test_case{0.13687969465365626, 2.257495744545807, 6.222102534434041, 0.12830374869142872},
            test_case{0.09983145651403719, 4.79589784180209, 3.1997444914326865, 0.7150705062507153},
            test_case{0.5965966533916149, 4.8160516653037995, 1.8717322673993892, 2.683021856447799}
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::gdtria(p, b, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, b, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("gdtrix basic", "[gdtrix][xsf_tests]") {

    SECTION("gdtrix double precision invalid and edge case inputs") {
        // test cases with invalid and edge case  inputs that should return NaN
        using test_case = std::tuple<double, double, double, double>;
        auto [a, b, p, expected] = GENERATE(
            test_case{
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            }, // all NaN
            test_case{
                0.5, std::numeric_limits<double>::quiet_NaN(), 0.2, std::numeric_limits<double>::quiet_NaN()
            }, // NaN b
            test_case{
                0.5, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            }, // NaN p
            test_case{
                std::numeric_limits<double>::quiet_NaN(), 1.0, 0.2, std::numeric_limits<double>::quiet_NaN()
            },                                                                   // NaN a
            test_case{0.0, 0.1, 0.2, std::numeric_limits<double>::infinity()},   //  a == 0
            test_case{-0.1, 1.0, 0.2, std::numeric_limits<double>::quiet_NaN()}, // a < 0
            test_case{0.5, 0.0, 0.2, std::numeric_limits<double>::quiet_NaN()},  // b== 0
            test_case{0.5, -0.1, 0.2, std::numeric_limits<double>::quiet_NaN()}, // b < 0
            test_case{1.0, 0.2, -0.1, std::numeric_limits<double>::quiet_NaN()}, // p < 0
            test_case{1.0, 0.2, 1.1, std::numeric_limits<double>::quiet_NaN()},  // p > 1
            test_case{0.1, 0.2, 0.0, 0.0},                                       // p == 0
            test_case{0.1, 0.2, 1.0, std::numeric_limits<double>::infinity()}    // p == 1
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::gdtrix(a, b, p);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(a, b, p, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("gdtrix double precision scipy reference values") {
        // Reference values computed with scipy.special.gdtrix using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // a_s = rng.uniform(0.1, 10.0, 30)
        // b_s = rng.uniform(0.1, 10.0, 30)
        // p_s = rng.uniform(0.01, 0.99, 30)

        // for a, b, p in zip(a_s, b_s, p_s):
        //     output = special.gdtrix(a, b, p)
        //     print(f"test_case{{{a}, {b}, {p}, {output}}}")
        using test_case = std::tuple<double, double, double, double>;
        auto [a, b, p, expected] = GENERATE(
            test_case{2.3506266224249797, 8.561944853302613, 0.26012770674311025, 2.7768649808487305},
            test_case{3.2359075631265535, 6.056050292767759, 0.3016139627850871, 1.4142751482065519},
            test_case{7.993918027594067, 9.326684775246237, 0.2834857774174913, 0.9245677330981279},
            test_case{6.794921240434648, 7.275335474810999, 0.26536762824147164, 0.7987280731149369},
            test_case{3.9719845509588994, 8.619458042193594, 0.4831063609433294, 2.0562951762532324},
            test_case{3.394857885877207, 9.30044423559563, 0.2177394556244804, 2.0143983080362546},
            test_case{6.023256660513179, 5.507241489915295, 0.49571798473957984, 0.8555973546319342},
            test_case{1.9486684374767624, 9.382962291800794, 0.2513360993141228, 3.6885953530599167},
            test_case{6.760284835744751, 5.00038060678036, 0.8317129994176059, 1.0432424645638272},
            test_case{9.423848366172377, 2.8103545066508766, 0.18652797829313436, 0.14269158874617255},
            test_case{2.5576325748327533, 4.5726092040001305, 0.8549131656790517, 2.6564332284783903},
            test_case{9.49392340314985, 6.68388534165535, 0.1847334559482837, 0.45928278226520325},
            test_case{6.7056507856936864, 3.375820211623841, 0.7455207052984992, 0.6472183910845487},
            test_case{1.0493895623817098, 9.044194667401566, 0.6088979957539539, 9.102551099940268},
            test_case{4.474212695061347, 2.64503433523769, 0.2149719342300352, 0.2968264901165944},
            test_case{8.876151201342426, 3.4643005423421664, 0.7546749727015153, 0.5079961463781446},
            test_case{7.004789648832019, 2.662648646564981, 0.25427535814421426, 0.21093548993098002},
            test_case{3.33208135429411, 3.6189201514484317, 0.09386029734682683, 0.43716056150799765},
            test_case{7.365888816967658, 0.1497211037996047, 0.6156955878717292, 0.003425034866074647},
            test_case{2.279336059900314, 6.323184986586819, 0.5362289644116889, 2.7277629765729774},
            test_case{0.9077862384678603, 2.895588803508671, 0.6318361769909702, 3.4609394813264274},
            test_case{1.6829664506429707, 0.7740681259306629, 0.1808866265175605, 0.06254468679910116},
            test_case{3.4669918310515824, 6.206606874838167, 0.2532012000593234, 1.273896295781778},
            test_case{4.705412221650304, 1.845630570783914, 0.6811265249466112, 0.46126266632639995},
            test_case{2.737568180078633, 3.1134450334739374, 0.08925421332588933, 0.4074812156303422},
            test_case{8.176186393905589, 4.464779427673569, 0.8675721287410036, 0.8338801248741946},
            test_case{2.0136144539659955, 1.5870031765207382, 0.430120493909192, 0.5388714456912779},
            test_case{1.3817438541542828, 2.257495744545807, 0.6160263114894303, 1.7135782486629199},
            test_case{1.0074810402948657, 4.79589784180209, 0.31684339410141743, 3.5136522999619384},
            test_case{6.025823335282641, 4.8160516653037995, 0.18538359818701025, 0.47437002825897256}
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::gdtrix(a, b, p);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(a, b, p, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
