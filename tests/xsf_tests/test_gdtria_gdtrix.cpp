#include "../testing_utils.h"
#include <cmath>
#include <limits>
#include <xsf/cdflib.h>

TEST_CASE("gdtria basic", "[gdtria][xsf_tests]") {

    SECTION("gdtria double precision invalid and edge case inputs") {
        // test cases with invalid and edge case  inputs
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
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::gdtria(p, b, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, b, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("gdtria float precision invalid and edge case inputs") {
        // test cases with invalid and edge case  inputs
        using test_case = std::tuple<float, float, float, float>;
        auto [p, b, x, expected] = GENERATE(
            test_case{
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(),
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            }, // all NaN
            test_case{
                0.5, std::numeric_limits<float>::quiet_NaN(), 0.2, std::numeric_limits<float>::quiet_NaN()
            }, // NaN b
            test_case{
                0.5, 1.0, std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            }, // NaN x
            test_case{
                std::numeric_limits<float>::quiet_NaN(), 1.0, 0.2, std::numeric_limits<float>::quiet_NaN()
            },                                                                  // NaN p
            test_case{0.5, 0.0, 0.2, std::numeric_limits<float>::quiet_NaN()},  //  b == 0
            test_case{0.5, -1.0, 0.2, std::numeric_limits<float>::quiet_NaN()}, // b < 0
            test_case{0.0, 1.0, 0.2, 0.0},                                      // p == 0
            test_case{1.0, 1.0, 0.2, std::numeric_limits<float>::infinity()},   // p == 1
            test_case{-0.1, 1.0, 0.2, std::numeric_limits<float>::quiet_NaN()}, // p < 0
            test_case{1.1, 1.0, 0.2, std::numeric_limits<float>::quiet_NaN()},  // p > 1
            test_case{0.5, 1.0, 0.0, std::numeric_limits<float>::quiet_NaN()}   // x == 0
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
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
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::gdtria(p, b, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, b, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("gdtria float precision scipy reference values") {
        // Reference values computed with scipy.special.gdtria using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // p_s = rng.uniform(0.01, 0.99, 30).astype(np.float32)
        // b_s = rng.uniform(0.1, 10.0, 30).astype(np.float32)
        // x_s = rng.uniform(-10.0, 10.0, 30).astype(np.float32)

        // for p, b, x in zip(p_s, b_s, x_s):
        //     output = special.gdtria(p, b, x)
        //     print(f"test_case{{{p}, {b}, {x}, {output}}}")
        using test_case = std::tuple<float, float, float, float>;
        auto [p, b, x, expected] = GENERATE(
            test_case{0.2327893078327179, 8.561944961547852, -4.895352840423584, -1.2905455827713013},
            test_case{0.3204231858253479, 6.0560503005981445, -4.048694610595703, -1.157485842704773},
            test_case{0.7914181351661682, 9.326684951782227, -4.4186577796936035, -2.6354193687438965},
            test_case{0.6727295517921448, 7.275335311889648, -4.788415908813477, -1.7096890211105347},
            test_case{0.3932873606681824, 8.619458198547363, -0.34476813673973083, -21.85759162902832},
            test_case{0.3361576497554779, 9.300444602966309, -5.7604193687438965, -1.3480381965637207},
            test_case{0.5963425636291504, 5.507241725921631, -0.0873880684375763, -65.79364013671875},
            test_case{0.19299949705600739, 9.382962226867676, -5.07477331161499, -1.3193113803863525},
            test_case{0.66930091381073, 5.000380516052246, 6.769652843475342, 0.8385697603225708},
            test_case{0.9329668283462524, 2.810354471206665, -6.397387981414795, -0.8774607181549072},
            test_case{0.25328078866004944, 4.5726094245910645, 7.243125915527344, 0.41787004470825195},
            test_case{0.9399035573005676, 6.68388557434082, -6.434010982513428, -1.7241863012313843},
            test_case{0.6638926863670349, 3.3758201599121094, 5.010626792907715, 0.7666146159172058},
            test_case{0.10397997498512268, 9.044194221496582, 2.2224080562591553, 2.481583833694458},
            test_case{0.4430028796195984, 2.6450343132019043, -5.816899299621582, -0.36294999718666077},
            test_case{0.8787503242492676, 3.4643006324768066, 5.197448253631592, 1.0896722078323364},
            test_case{0.6935044527053833, 2.6626486778259277, -5.014788627624512, -0.6357687711715698},
            test_case{0.3299434185028076, 3.618920087814331, -8.288565635681152, -0.3091883361339569},
            test_case{0.7292495965957642, 0.14972110092639923, 2.3611345291137695, 0.03473275154829025},
            test_case{0.2257322520017624, 6.323184967041016, 0.7393665909767151, 5.871931076049805},
            test_case{0.0899626761674881, 2.8955888748168945, 2.6905341148376465, 0.36772969365119934},
            test_case{0.166697695851326, 0.7740681171417236, -6.512517929077148, -0.014454977586865425},
            test_case{0.34329816699028015, 6.206606864929199, -5.036710262298584, -0.9840686321258545},
            test_case{0.4658893048763275, 1.8456305265426636, 3.6964597702026367, 0.3851713538169861},
            test_case{0.27109259366989136, 3.1134450435638428, -8.382567405700684, -0.22649358212947845},
            test_case{0.8094608783721924, 4.464779376983643, 7.501471996307373, 0.8222838640213013},
            test_case{0.1994284987449646, 1.587003231048584, -1.426112413406372, -0.3891430199146271},
            test_case{0.13687969744205475, 2.257495641708374, 2.3678839206695557, 0.33714449405670166},
            test_case{0.09983145445585251, 4.795897960662842, -3.7378900051116943, -0.6121215224266052},
            test_case{0.5965966582298279, 4.816051483154297, -6.420742988586426, -0.7821366786956787}
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto output = xsf::gdtria(p, b, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, b, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("gdtrix basic", "[gdtrix][xsf_tests]") {

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
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::gdtrix(a, b, p);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(a, b, p, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("gdtrix float precision invalid and edge case inputs") {
        // test cases with invalid and edge case  inputs
        using test_case = std::tuple<float, float, float, float>;
        auto [a, b, p, expected] = GENERATE(
            test_case{
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(),
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            }, // all NaN
            test_case{
                0.5, std::numeric_limits<float>::quiet_NaN(), 0.2, std::numeric_limits<float>::quiet_NaN()
            }, // NaN b
            test_case{
                0.5, 1.0, std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            }, // NaN p
            test_case{
                std::numeric_limits<float>::quiet_NaN(), 1.0, 0.2, std::numeric_limits<float>::quiet_NaN()
            },                                                                  // NaN a
            test_case{0.0, 0.1, 0.2, std::numeric_limits<float>::infinity()},   //  a == 0
            test_case{0.5, 0.0, 0.2, std::numeric_limits<float>::quiet_NaN()},  // b == 0
            test_case{0.5, -0.1, 0.2, std::numeric_limits<float>::quiet_NaN()}, // b < 0
            test_case{1.0, 0.2, -0.1, std::numeric_limits<float>::quiet_NaN()}, // p < 0
            test_case{1.0, 0.2, 1.1, std::numeric_limits<float>::quiet_NaN()},  // p > 1
            test_case{0.1, 0.2, 0.0, 0.0},                                      // p == 0
            test_case{0.1, 0.2, 1.0, std::numeric_limits<float>::infinity()}    // p == 1
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
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
        const auto output = xsf::gdtrix(a, b, p);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(a, b, p, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("gdtrix float precision scipy reference values") {
        // Reference values computed with scipy.special.gdtrix using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // a_s = rng.uniform(-10.0, 10.0, 30).astype(np.float32)
        // b_s = rng.uniform(0.1, 10.0, 30).astype(np.float32)
        // p_s = rng.uniform(0.01, 0.99, 30).astype(np.float32)

        // for a, b, p in zip(a_s, b_s, p_s):
        //     output = special.gdtrix(a, b, p)
        //     print(f"test_case{{{a}, {b}, {p}, {output}}}")
        using test_case = std::tuple<float, float, float, float>;
        auto [a, b, p, expected] = GENERATE(
            test_case{-9.805757522583008, 3.381485939025879, 0.7704927921295166, -0.4613720774650574},
            test_case{-5.799140930175781, 3.843750476837158, 0.8609054684638977, -1.024255633354187},
            test_case{7.400013446807861, 3.5338046550750732, 0.477615088224411, 0.4200052320957184},
            test_case{9.456596374511719, 5.210931301116943, 0.17547166347503662, 0.3288891315460205},
            test_case{-1.1641530990600586, 0.189040869474411, 0.5332977771759033, -0.02040097489953041},
            test_case{-2.4250102043151855, 4.284514427185059, 0.42731523513793945, -1.4867595434188843},
            test_case{-4.481058597564697, 8.788811683654785, 0.16976900398731232, -1.3353513479232788},
            test_case{9.32208251953125, 0.9653109908103943, 0.17256002128124237, 0.01868532784283161},
            test_case{-8.83594799041748, 4.892439842224121, 0.2629871666431427, -0.37878063321113586},
            test_case{-1.8253220319747925, 4.86415433883667, 0.9342297315597534, -4.668603897094727},
            test_case{-6.627423286437988, 7.8474578857421875, 0.9770708084106445, -2.168081045150757},
            test_case{-5.197118759155273, 9.649139404296875, 0.6995930671691895, -2.1170544624328613},
            test_case{5.600157260894775, 7.100254535675049, 0.718694269657135, 1.4969807863235474},
            test_case{-5.924648284912109, 2.8099935054779053, 0.7278478741645813, -0.599925696849823},
            test_case{1.041019082069397, 6.734121799468994, 0.5170475840568542, 6.25558614730835},
            test_case{-2.6601171493530273, 3.5405945777893066, 0.18387778103351593, -0.7025629281997681},
            test_case{0.14563442766666412, 7.704465389251709, 0.8700932860374451, 74.55491638183594},
            test_case{-3.3312439918518066, 6.790136814117432, 0.8727825284004211, -2.938950538635254},
            test_case{-4.345566749572754, 9.777566909790039, 0.7053438425064087, -2.577073335647583},
            test_case{-4.363393783569336, 8.680426597595215, 0.9247602224349976, -3.0261147022247314},
            test_case{-8.292374610900879, 0.5564693212509155, 0.9703981876373291, -0.30314576625823975},
            test_case{-0.3637268841266632, 2.9742047786712646, 0.503406286239624, -7.318968772888184},
            test_case{7.666857719421387, 8.63767147064209, 0.7524181008338928, 1.3595703840255737},
            test_case{8.944555282592773, 6.048393726348877, 0.6142861247062683, 0.7197775840759277},
            test_case{-9.452325820922852, 3.5081536769866943, 0.1225912868976593, -0.1635589599609375},
            test_case{8.355045318603516, 0.6504654884338379, 0.319668710231781, 0.019422415643930435},
            test_case{-7.569509506225586, 7.652439594268799, 0.0763792172074318, -0.5435218214988708},
            test_case{4.956955432891846, 0.1150359958410263, 0.8791337609291077, 0.04991769790649414},
            test_case{7.93041467666626, 1.0191510915756226, 0.03162899240851402, 0.004364134278148413},
            test_case{-6.641404151916504, 9.935690879821777, 0.5365409851074219, -1.489241123199463}
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto output = xsf::gdtrix(a, b, p);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(a, b, p, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
