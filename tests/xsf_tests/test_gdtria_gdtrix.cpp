#include "../testing_utils.h"
#include <cmath>
#include <limits>
#include <xsf/cdflib.h>

TEST_CASE("gdtria basic", "[gdtria][xsf_tests]") {
    // checks for both double and float precision
    SECTION("gdtria invalid and edge case inputs") {
        // test cases with invalid and edge case inputs
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
            test_case{0.5, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()}   // x == 0
        );
        const double rtol_double = 100 * std::numeric_limits<double>::epsilon();
        const auto output_double = xsf::gdtria(p, b, x);
        const auto rel_error_double = xsf::extended_relative_error(output_double, expected);
        CAPTURE(p, b, x, output_double, expected, rtol_double, rel_error_double);
        REQUIRE(rel_error_double <= rtol_double);

        const float rtol_float = 100 * std::numeric_limits<float>::epsilon();
        const auto output_float = xsf::gdtria(static_cast<float>(p), static_cast<float>(b), static_cast<float>(x));
        const auto rel_error_float = xsf::extended_relative_error(output_float, static_cast<float>(expected));
        CAPTURE(p, b, x, output_float, expected, rtol_float, rel_error_float);
        REQUIRE(rel_error_float <= rtol_float);
    }

    SECTION("gdtria scipy reference values") {
        // Reference values computed with scipy.special.gdtria using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // p_s = rng.uniform(0.01, 0.99, 30)
        // b_s = rng.uniform(0.1, 10.0, 30)
        // x_s = rng.uniform(-10.0, 10.0, 30)

        // for p, b, x in zip(p_s, b_s, x_s):
        //     output = special.gdtria(p, b, x)
        //     print(f"test_case{{{p}, {b}, {x}, {output}}}")
        using test_case = std::tuple<double, double, double, double>;
        auto [p, b, x, expected] = GENERATE(
            test_case{0.23278930201782627, 8.561944853302613, -4.895352923609995, -1.290545581913255},
            test_case{0.32042317291555783, 6.056050292767759, -4.048694637039039, -1.1574858515650308},
            test_case{0.7914181481860795, 9.326684775246237, -4.418657603724667, -2.6354195392837925},
            test_case{0.6727295773359551, 7.275335474810999, -4.788415750174049, -1.7096891898701405},
            test_case{0.3932873595898708, 8.619458042193594, -0.3447681440136847, -21.85759161598062},
            test_case{0.33615764930905684, 9.30044423559563, -5.760419272969788, -1.3480381097736163},
            test_case{0.596342578515446, 5.507241489915295, -0.08738806653918729, -65.79364358100551},
            test_case{0.19299950189163909, 9.382962291800794, -5.074773483385249, -1.3193113188933627},
            test_case{0.6693009231343289, 5.00038060678036, 6.769653049338896, 0.8385697480377493},
            test_case{0.9329668079645385, 2.8103545066508766, -6.397388198099298, -0.8774606282295611},
            test_case{0.2532808003369796, 4.5726092040001305, 7.243125830184731, 0.41787002099682563},
            test_case{0.9399035287966518, 6.68388534165535, -6.4340111030962515, -1.7241861296982413},
            test_case{0.663892704038365, 3.375820211623841, 5.010626638744881, 0.7666146732226315},
            test_case{0.10397997688222983, 9.044194667401566, 2.2224080766113037, 2.481584084072316},
            test_case{0.44300287284445655, 2.64503433523769, -5.816899301427854, -0.362949989844559},
            test_case{0.8787503209409673, 3.4643005423421664, 5.197448422479905, 1.089672137381579},
            test_case{0.6935044298843817, 2.662648646564981, -5.01478860930175, -0.6357687432051166},
            test_case{0.32994340678870987, 3.6189201514484317, -8.28856536026884, -0.309188347644603},
            test_case{0.7292496000634652, 0.1497211037996047, 2.3611344463618185, 0.03473275456543726},
            test_case{0.2257322564345765, 6.323184986586819, 0.7393666206467113, 5.871930800120954},
            test_case{0.08996267815136395, 2.895588803508671, 2.690534224305514, 0.367729655044082},
            test_case{0.16669768905354657, 0.7740681259306629, -6.512517826172235, -0.01445497793499551},
            test_case{0.3432981812556112, 6.206606874838167, -5.036710202870951, -0.9840686543905368},
            test_case{0.46588929062800993, 1.845630570783914, 3.6964596927879825, 0.38517135128472085},
            test_case{0.2710926077249555, 3.1134450334739374, -8.38256707498185, -0.22649359426243212},
            test_case{0.8094608753563108, 4.464779427673569, 7.501472015122523, 0.8222838837552521},
            test_case{0.1994285015037046, 1.5870031765207382, -1.4261123692001636, -0.38914300007735786},
            test_case{0.13687969465365626, 2.257495744545807, 2.367883907947556, 0.33714451845838045},
            test_case{0.09983145651403719, 4.79589784180209, -3.737889916297603, -0.6121215350365443},
            test_case{0.5965966533916149, 4.8160516653037995, -6.4207428941426485, -0.7821366881755081}
        );
        const double rtol_double = 100 * std::numeric_limits<double>::epsilon();
        const auto output_double = xsf::gdtria(p, b, x);
        const auto rel_error_double = xsf::extended_relative_error(output_double, expected);
        CAPTURE(p, b, x, output_double, expected, rtol_double, rel_error_double);
        REQUIRE(rel_error_double <= rtol_double);

        const float rtol_float = 100 * std::numeric_limits<float>::epsilon();
        const auto output_float = xsf::gdtria(static_cast<float>(p), static_cast<float>(b), static_cast<float>(x));
        const auto rel_error_float = xsf::extended_relative_error(output_float, static_cast<float>(expected));
        CAPTURE(p, b, x, output_float, expected, rtol_float, rel_error_float);
        REQUIRE(rel_error_float <= rtol_float);
    }
}

TEST_CASE("gdtrix basic", "[gdtrix][xsf_tests]") {
    // checks for both double and float precision
    SECTION("gdtrix double precision invalid and edge case inputs") {
        // test cases with invalid and edge case  inputs
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
            test_case{0.5, 0.0, 0.2, std::numeric_limits<double>::quiet_NaN()},  // b == 0
            test_case{0.5, -0.1, 0.2, std::numeric_limits<double>::quiet_NaN()}, // b < 0
            test_case{1.0, 0.2, -0.1, std::numeric_limits<double>::quiet_NaN()}, // p < 0
            test_case{1.0, 0.2, 1.1, std::numeric_limits<double>::quiet_NaN()},  // p > 1
            test_case{0.1, 0.2, 0.0, 0.0},                                       // p == 0
            test_case{0.1, 0.2, 1.0, std::numeric_limits<double>::infinity()}    // p == 1
        );
        const double rtol_double = 100 * std::numeric_limits<double>::epsilon();
        const auto output_double = xsf::gdtrix(a, b, p);
        const auto rel_error_double = xsf::extended_relative_error(output_double, expected);
        CAPTURE(a, b, p, output_double, expected, rtol_double, rel_error_double);
        REQUIRE(rel_error_double <= rtol_double);

        const float rtol_float = 100 * std::numeric_limits<float>::epsilon();
        const auto output_float = xsf::gdtrix(static_cast<float>(a), static_cast<float>(b), static_cast<float>(p));
        const auto rel_error_float = xsf::extended_relative_error(output_float, static_cast<float>(expected));
        CAPTURE(a, b, p, output_float, expected, rtol_float, rel_error_float);
        REQUIRE(rel_error_float <= rtol_float);
    }

    SECTION("gdtrix double precision scipy reference values") {
        // Reference values computed with scipy.special.gdtrix using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // a_s = rng.uniform(-10.0, 10.0, 30)
        // b_s = rng.uniform(0.1, 10.0, 30)
        // p_s = rng.uniform(0.01, 0.99, 30)

        // for a, b, p in zip(a_s, b_s, p_s):
        //     output = special.gdtrix(a, b, p)
        //     print(f"test_case{{{a}, {b}, {p}, {output}}}")
        using test_case = std::tuple<double, double, double, double>;
        auto [a, b, p, expected] = GENERATE(
            test_case{-5.453279550656607, 8.561944853302613, 0.26012770674311025, -1.196962798299369},
            test_case{-3.6648332058049427, 6.056050292767759, 0.3016139627850871, -1.2487508684364081},
            test_case{5.947309146654682, 9.326684775246237, 0.2834857774174913, 1.242733223898785},
            test_case{3.525093415019491, 7.275335474810999, 0.26536762824147164, 1.5396171704885484},
            test_case{-2.1778089879618197, 8.619458042193594, 0.4831063609433294, -3.750362275771973},
            test_case{-3.34372144267231, 9.30044423559563, 0.2177394556244804, -2.045205050295446},
            test_case{1.9661750717437965, 5.507241489915295, 0.49571798473957984, 2.6210699845938144},
            test_case{-6.265316287925733, 9.382962291800794, 0.2513360993141228, -1.147244450688539},
            test_case{3.4551208802924265, 5.00038060678036, 0.8317129994176059, 2.041206793493986},
            test_case{8.836057305398743, 2.8103545066508766, 0.18652797829313436, 0.15218370014990265},
            test_case{-5.03508570740858, 4.5726092040001305, 0.8549131656790517, -1.3493673301385065},
            test_case{8.977623036666365, 6.68388534165535, 0.1847334559482837, 0.48569599518744366},
            test_case{3.3447490620074483, 3.375820211623841, 0.7455207052984992, 1.2975623678287898},
            test_case{-8.082041288117757, 9.044194667401566, 0.6088979957539539, -1.1818947435182037},
            test_case{-1.1632066766437443, 2.64503433523769, 0.2149719342300352, -1.141727327547755},
            test_case{7.729598386550354, 3.4643005423421664, 0.7546749727015153, 0.5833486268571888},
            test_case{3.949069997640443, 2.662648646564981, 0.25427535814421426, 0.3741535949787358},
            test_case{-3.4705427185977573, 3.6189201514484317, 0.09386029734682683, -0.4197195291755668},
            test_case{4.67856326660133, 0.1497211037996047, 0.6156955878717292, 0.005392344739215286},
            test_case{-5.597300889090276, 6.323184986586819, 0.5362289644116889, -1.110801194819156},
            test_case{-8.368108609155838, 2.895588803508671, 0.6318361769909702, -0.3754484292759633},
            test_case{-6.802087978499049, 0.7740681259306629, 0.1808866265175605, -0.015474749794707358},
            test_case{-3.197996300905894, 6.206606874838167, 0.2532012000593234, -1.3810485177331855},
            test_case{-0.6961369259589816, 1.845630570783914, 0.6811265249466112, -3.1178219493719785},
            test_case{-4.671579434184581, 3.1134450334739374, 0.08925421332588933, -0.23878596641782146},
            test_case{6.315528068496139, 4.464779427673569, 0.8675721287410036, 1.0795549092964738},
            test_case{-6.13411221421011, 1.5870031765207382, 0.430120493909192, -0.17689264460468218},
            test_case{-7.410618476455994, 2.257495744545807, 0.6160263114894303, -0.31950453544800145},
            test_case{-8.166704969101282, 4.79589784180209, 0.31684339410141743, -0.43345977206149267},
            test_case{1.971360273298263, 4.8160516653037995, 0.18538359818701025, 1.4499987772702376}
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output_double = xsf::gdtrix(a, b, p);
        const auto rel_error_double = xsf::extended_relative_error(output_double, expected);
        CAPTURE(a, b, p, output_double, expected, rtol, rel_error_double);
        REQUIRE(rel_error_double <= rtol);

        const float rtol_float = 100 * std::numeric_limits<float>::epsilon();
        const auto output_float = xsf::gdtrix(static_cast<float>(a), static_cast<float>(b), static_cast<float>(p));
        const auto rel_error_float = xsf::extended_relative_error(output_float, static_cast<float>(expected));
        CAPTURE(a, b, p, output_float, expected, rtol_float, rel_error_float);
        REQUIRE(rel_error_float <= rtol_float);
    }
}
