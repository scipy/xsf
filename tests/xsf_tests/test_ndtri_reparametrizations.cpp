#include "../testing_utils.h"
#include <xsf/stats.h>

#include <cmath>
#include <limits>

TEST_CASE("nrdtrimn test", "[nrdtrimn][xsf_tests]") {

    SECTION("nrdtrimn double precision invalid inputs") {
        // test cases with invalid inputs that should return NaN
        using test_case = std::tuple<double, double, double, double>;
        auto [p, std, x, expected] = GENERATE(
            test_case{
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            }, // all NaN
            test_case{
                0.5, std::numeric_limits<double>::quiet_NaN(), 0.0, std::numeric_limits<double>::quiet_NaN()
            }, // NaN std
            test_case{
                0.5, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            }, // NaN x
            test_case(
                std::numeric_limits<double>::quiet_NaN(), 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()
            ),                                                                   // NaN p
            test_case{0.5, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN()},  // std == 0
            test_case{0.5, -1.0, 0.0, std::numeric_limits<double>::quiet_NaN()}, // std < 0
            test_case{0.0, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()},  // p == 0
            test_case{1.0, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()},  // p == 1
            test_case(-0.1, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()), // p < 0
            test_case(1.1, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN())   // p > 1
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::nrdtrimn(p, std, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, std, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("nrdtrimn float precision invalid inputs") {
        // test cases with invalid inputs that should return NaN
        using test_case = std::tuple<float, float, float, float>;
        auto [p, std, x, expected] = GENERATE(
            test_case{
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(),
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            }, // all NaN
            test_case{
                0.5, std::numeric_limits<float>::quiet_NaN(), 0.0, std::numeric_limits<float>::quiet_NaN()
            }, // NaN std
            test_case{
                0.5, 1.0, std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            }, // NaN x
            test_case(
                std::numeric_limits<float>::quiet_NaN(), 1.0, 0.0, std::numeric_limits<float>::quiet_NaN()
            ),                                                                  // NaN p
            test_case{0.5, 0.0, 0.0, std::numeric_limits<float>::quiet_NaN()},  // std == 0
            test_case{0.5, -1.0, 0.0, std::numeric_limits<float>::quiet_NaN()}, // std < 0
            test_case{0.0, 1.0, 0.0, std::numeric_limits<float>::quiet_NaN()},  // p == 0
            test_case{1.0, 1.0, 0.0, std::numeric_limits<float>::quiet_NaN()},  // p == 1
            test_case(-0.1, 1.0, 0.0, std::numeric_limits<float>::quiet_NaN()), // p < 0
            test_case(1.1, 1.0, 0.0, std::numeric_limits<float>::quiet_NaN())   // p > 1
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto output = xsf::nrdtrimn(p, std, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, std, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("nrdtrimn double precision scipy reference values") {
        // Reference values computed with scipy.special.nrdtrimn using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // p_s = rng.uniform(0.01, 0.99, 30)
        // std_s = rng.uniform(0.1, 10.0, 30)
        // x_s = rng.uniform(-10.0, 10.0, 30)

        // for p, std, x in zip(p_s, std_s, x_s):
        //     output = special.nrdtrimn(p, std, x)
        //     print(f"test_case{{{p}, {std}, {x}, {output}}}")
        using test_case = std::tuple<double, double, double, double>;
        auto [p, std, x, expected] = GENERATE(
            test_case{0.23278930201782627, 8.561944853302613, -4.895352923609995, 1.3522278947568704},
            test_case{0.32042317291555783, 6.056050292767759, -4.048694637039039, -1.2234515325366395},
            test_case{0.7914181481860795, 9.326684775246237, -4.418657603724667, -11.985879548048604},
            test_case{0.6727295773359551, 7.275335474810999, -4.788415750174049, -8.043858701021936},
            test_case{0.3932873595898708, 8.619458042193594, -0.3447681440136847, 1.9890464820009215},
            test_case{0.33615764930905684, 9.30044423559563, -5.760419272969788, -1.8265867842421923},
            test_case{0.596342578515446, 5.507241489915295, -0.08738806653918729, -1.430557521821209},
            test_case{0.19299950189163909, 9.382962291800794, -5.074773483385249, 3.05927885166947},
            test_case{0.6693009231343289, 5.00038060678036, 6.769653049338896, 4.5795682242834514},
            test_case{0.9329668079645385, 2.8103545066508766, -6.397388198099298, -10.608022671395918},
            test_case{0.2532808003369796, 4.5726092040001305, 7.243125830184731, 10.2802579566605},
            test_case{0.9399035287966518, 6.68388534165535, -6.4340111030962515, -16.82053003113355},
            test_case{0.663892704038365, 3.375820211623841, 5.010626638744881, 3.582281429837138},
            test_case{0.10397997688222983, 9.044194667401566, 2.2224080766113037, 13.610811612312046},
            test_case{0.44300287284445655, 2.64503433523769, -5.816899301427854, -5.437706756134887},
            test_case{0.8787503209409673, 3.4643005423421664, 5.197448422479905, 1.1485040763743886},
            test_case{0.6935044298843817, 2.662648646564981, -5.01478860930175, -6.3615787400132735},
            test_case{0.32994340678870987, 3.6189201514484317, -8.28856536026884, -6.695989189721861},
            test_case{0.7292496000634652, 0.1497211037996047, 2.3611344463618185, 2.2697229651590494},
            test_case{0.2257322564345765, 6.323184986586819, 0.7393666206467113, 5.500571300095614},
            test_case{0.08996267815136395, 2.895588803508671, 2.690534224305514, 6.573475076898129},
            test_case{0.16669768905354657, 0.7740681259306629, -6.512517826172235, -5.763763732811857},
            test_case{0.3432981812556112, 6.206606874838167, -5.036710202870951, -2.532478735787074},
            test_case{0.46588929062800993, 1.845630570783914, 3.6964596927879825, 3.8544591751368613},
            test_case{0.2710926077249555, 3.1134450334739374, -8.38256707498185, -6.484885408582596},
            test_case{0.8094608753563108, 4.464779427673569, 7.501472015122523, 3.59072120751771},
            test_case{0.1994285015037046, 1.5870031765207382, -1.4261123692001636, -0.08721439084858873},
            test_case{0.13687969465365626, 2.257495744545807, 2.367883907947556, 4.838591246064786},
            test_case{0.09983145651403719, 4.79589784180209, -3.737889916297603, 2.412909149035708},
            test_case{0.5965966533916149, 4.8160516653037995, -6.4207428941426485, -7.598497069341223}
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::nrdtrimn(p, std, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, std, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("nrdtrimn float precision scipy reference values") {
        // Reference values computed with scipy.special.nrdtrimn using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // p_s = rng.uniform(0.01, 0.99, 30).astype(np.float32)
        // std_s = rng.uniform(0.1, 10.0, 30).astype(np.float32)
        // x_s = rng.uniform(-10.0, 10.0, 30).astype(np.float32)

        // for p, std, x in zip(p_s, std_s, x_s):
        //     output = special.nrdtrimn(p, std, x)
        //     print(f"test_case{{{p}, {std}, {x}, {output}}}")
        using test_case = std::tuple<float, float, float, float>;
        auto [p, std, x, expected] = GENERATE(
            test_case{0.2327893078327179, 8.561944961547852, -4.895352840423584, 1.3522279262542725},
            test_case{0.3204231858253479, 6.0560503005981445, -4.048694610595703, -1.2234517335891724},
            test_case{0.7914181351661682, 9.326684951782227, -4.4186577796936035, -11.985879898071289},
            test_case{0.6727295517921448, 7.275335311889648, -4.788415908813477, -8.043858528137207},
            test_case{0.3932873606681824, 8.619458198547363, -0.34476813673973083, 1.9890464544296265},
            test_case{0.3361576497554779, 9.300444602966309, -5.7604193687438965, -1.8265867233276367},
            test_case{0.5963425636291504, 5.507241725921631, -0.0873880684375763, -1.430557370185852},
            test_case{0.19299949705600739, 9.382962226867676, -5.07477331161499, 3.059279203414917},
            test_case{0.66930091381073, 5.000380516052246, 6.769652843475342, 4.579568386077881},
            test_case{0.9329668283462524, 2.810354471206665, -6.397387981414795, -10.608022689819336},
            test_case{0.25328078866004944, 4.5726094245910645, 7.243125915527344, 10.280258178710938},
            test_case{0.9399035573005676, 6.68388557434082, -6.434010982513428, -16.820531845092773},
            test_case{0.6638926863670349, 3.3758201599121094, 5.010626792907715, 3.5822818279266357},
            test_case{0.10397997498512268, 9.044194221496582, 2.2224080562591553, 13.610811233520508},
            test_case{0.4430028796195984, 2.6450343132019043, -5.816899299621582, -5.43770694732666},
            test_case{0.8787503242492676, 3.4643006324768066, 5.197448253631592, 1.1485037803649902},
            test_case{0.6935044527053833, 2.6626486778259277, -5.014788627624512, -6.361578941345215},
            test_case{0.3299434185028076, 3.618920087814331, -8.288565635681152, -6.695989608764648},
            test_case{0.7292495965957642, 0.14972110092639923, 2.3611345291137695, 2.2697229385375977},
            test_case{0.2257322520017624, 6.323184967041016, 0.7393665909767151, 5.500571250915527},
            test_case{0.0899626761674881, 2.8955888748168945, 2.6905341148376465, 6.573474884033203},
            test_case{0.166697695851326, 0.7740681171417236, -6.512517929077148, -5.763763904571533},
            test_case{0.34329816699028015, 6.206606864929199, -5.036710262298584, -2.5324785709381104},
            test_case{0.4658893048763275, 1.8456305265426636, 3.6964597702026367, 3.854459285736084},
            test_case{0.27109259366989136, 3.1134450435638428, -8.382567405700684, -6.4848856925964355},
            test_case{0.8094608783721924, 4.464779376983643, 7.501471996307373, 3.5907211303710938},
            test_case{0.1994284987449646, 1.587003231048584, -1.426112413406372, -0.08721437305212021},
            test_case{0.13687969744205475, 2.257495641708374, 2.3678839206695557, 4.8385910987854},
            test_case{0.09983145445585251, 4.795897960662842, -3.7378900051116943, 2.4129092693328857},
            test_case{0.5965966582298279, 4.816051483154297, -6.420742988586426, -7.59849739074707}
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto output = xsf::nrdtrimn(p, std, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(p, std, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}

TEST_CASE("nrdtrisd test", "[nrdtrisd][xsf_tests]") {

    SECTION("nrdtrisd double precision invalid inputs") {
        // test cases with invalid inputs that should return NaN
        using test_case = std::tuple<double, double, double, double>;
        auto [mean, p, x, expected] = GENERATE(
            test_case{
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            }, // all NaN
            test_case{
                0.0, std::numeric_limits<double>::quiet_NaN(), 0.0, std::numeric_limits<double>::quiet_NaN()
            }, // NaN p
            test_case{
                std::numeric_limits<double>::quiet_NaN(), 0.5, 0.0, std::numeric_limits<double>::quiet_NaN()
            }, // NaN mean
            test_case{
                0.0, 0.5, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
            },                                                                   // NaN x
            test_case{0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN()},  // p == 0
            test_case{0.0, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN()},  // p == 1
            test_case{0.0, -0.1, 0.0, std::numeric_limits<double>::quiet_NaN()}, // p < 0
            test_case{0.0, 1.1, 0.0, std::numeric_limits<double>::quiet_NaN()}   // p > 1
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::nrdtrisd(mean, p, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(mean, p, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("nrdtrisd float precision invalid inputs") {
        // test cases with invalid inputs that should return NaN
        using test_case = std::tuple<float, float, float, float>;
        auto [mean, p, x, expected] = GENERATE(
            test_case{
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(),
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            }, // all NaN
            test_case{
                0.0, std::numeric_limits<float>::quiet_NaN(), 0.0, std::numeric_limits<float>::quiet_NaN()
            }, // NaN p
            test_case{
                std::numeric_limits<float>::quiet_NaN(), 0.5, 0.0, std::numeric_limits<float>::quiet_NaN()
            }, // NaN mean
            test_case{
                0.0, 0.5, std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()
            },                                                                  // NaN x
            test_case{0.0, 0.0, 0.0, std::numeric_limits<float>::quiet_NaN()},  // p == 0
            test_case{0.0, 1.0, 0.0, std::numeric_limits<float>::quiet_NaN()},  // p == 1
            test_case{0.0, -0.1, 0.0, std::numeric_limits<float>::quiet_NaN()}, // p < 0
            test_case{0.0, 1.1, 0.0, std::numeric_limits<float>::quiet_NaN()}   // p > 1
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto output = xsf::nrdtrisd(mean, p, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(mean, p, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("nrdtrisd double precision scipy reference values") {
        // Reference values computed with scipy.special.nrdtrisd using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // mean_s = rng.uniform(-10.0, 10.0, 30)
        // p_s = rng.uniform(0.01, 0.99, 30)
        // x_s = rng.uniform(-10.0, 10.0, 30)

        // for mean, p, x in zip(mean_s, p_s, x_s):
        //     output = special.nrdtrisd(mean, p, x)
        //     print(f"test_case{{{mean}, {p}, {x}, {output}}}")
        using test_case = std::tuple<double, double, double, double>;
        auto [mean, p, x, expected] = GENERATE(
            test_case{-5.453279550656607, 0.8476470662865213, -4.895352923609995, 0.5435793502613633},
            test_case{-3.6648332058049427, 0.5995888168598389, -4.048694637039039, -1.5215512261933377},
            test_case{5.947309146654682, 0.9233485939132638, -4.418657603724667, -7.259273810304403},
            test_case{3.525093415019491, 0.7202857338701797, -4.788415750174049, -14.24300820343293},
            test_case{-2.1778089879618197, 0.8533402910454265, -0.3447681440136847, 1.7443119387869397},
            test_case{-3.34372144267231, 0.9207510455438099, -5.760419272969788, -1.713798119419246},
            test_case{1.9661750717437965, 0.5452622889007059, -0.08738806653918729, -18.06120318476792},
            test_case{-6.265316287925733, 0.9289194995924017, -5.074773483385249, 0.8111119350581782},
            test_case{3.4551208802924265, 0.4950881812772478, 6.769653049338896, -269.20244578782393},
            test_case{8.836057305398743, 0.2782977188401878, -6.397388198099298, 25.911365737533085},
            test_case{-5.03508570740858, 0.45274313332526545, 7.243125830184731, -103.40958296344068},
            test_case{8.977623036666365, 0.6617381449315397, -6.4340111030962515, -36.939619181983346},
            test_case{3.3447490620074483, 0.3342731118577136, 5.010626638744881, -3.890927476568311},
            test_case{-8.082041288117757, 0.8953849266720743, 2.2224080766113037, 8.206236566723723},
            test_case{-1.1632066766437443, 0.26193269177100365, -5.816899301427854, 7.30107382310783},
            test_case{7.729598386550354, 0.3430317708581134, 5.197448422479905, 6.264552127929731},
            test_case{3.949069997640443, 0.2636763306700688, -5.01478860930175, 14.18214553813266},
            test_case{-3.4705427185977573, 0.3583375503454003, -8.28856536026884, 13.27622737982676},
            test_case{4.67856326660133, 0.014921887042789152, 2.3611344463618185, 1.0668787434364317},
            test_case{-5.597300889090276, 0.6260324532176851, 0.7393666206467113, 19.718080938103416},
            test_case{-8.368108609155838, 0.2867350532766159, 2.690534224305514, -19.64415441534154},
            test_case{-6.802087978499049, 0.07672593569818682, -6.512517826172235, -0.20285916144950464},
            test_case{-3.197996300905894, 0.614492397711253, -5.036710202870951, -6.3175793046662605},
            test_case{-0.6961369259589816, 0.18279979387557937, 3.6964596927879825, -4.855056732853751},
            test_case{-4.671579434184581, 0.3083006194751978, -8.38256707498185, 7.411997925181353},
            test_case{6.315528068496139, 0.4420690746585957, 7.501472015122523, -8.138208867965012},
            test_case{-6.13411221421011, 0.1571982942414468, -1.4261123692001636, -4.679736760057691},
            test_case{-7.410618476455994, 0.22357028582372632, 2.367883907947556, -12.86322195904312},
            test_case{-8.166704969101282, 0.4748464530268736, -3.737889916297603, -70.19566984245017},
            test_case{1.971360273298263, 0.47684147797956805, -6.4207428941426485, 144.48604671757812}
        );
        const double rtol = 100 * std::numeric_limits<double>::epsilon();
        const auto output = xsf::nrdtrisd(mean, p, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(mean, p, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }

    SECTION("nrdtrisd float precision scipy reference values") {
        // Reference values computed with scipy.special.nrdtrisd using:
        // import numpy as np
        // from scipy import special

        // rng= np.random.default_rng(12345)

        // mean_s = rng.uniform(-10.0, 10.0, 30).astype(np.float32)
        // p_s = rng.uniform(0.01, 0.99, 30).astype(np.float32)
        // x_s = rng.uniform(-10.0, 10.0, 30).astype(np.float32)

        // for mean, p, x in zip(mean_s, p_s, x_s):
        //     output = special.nrdtrisd(mean, p, x)
        //     print(f"test_case{{{mean}, {p}, {x}, {output}}}")
        using test_case = std::tuple<float, float, float, float>;
        auto [mean, p, x, expected] = GENERATE(
            test_case{-5.453279495239258, 0.8476470708847046, -4.895352840423584, 0.5435793399810791},
            test_case{-3.6648333072662354, 0.5995888113975525, -4.048694610595703, -1.5215507745742798},
            test_case{5.9473090171813965, 0.923348605632782, -4.4186577796936035, -7.259273529052734},
            test_case{3.5250933170318604, 0.7202857136726379, -4.788415908813477, -14.243009567260742},
            test_case{-2.177809000015259, 0.8533402681350708, -0.34476813673973083, 1.7443121671676636},
            test_case{-3.343721389770508, 0.9207510352134705, -5.7604193687438965, -1.7137982845306396},
            test_case{1.9661750793457031, 0.5452622771263123, -0.0873880684375763, -18.061208724975586},
            test_case{-6.265316486358643, 0.9289194941520691, -5.07477331161499, 0.8111122250556946},
            test_case{3.455120801925659, 0.4950881898403168, 6.769652843475342, -269.2029113769531},
            test_case{8.836057662963867, 0.2782977223396301, -6.397387981414795, 25.911365509033203},
            test_case{-5.035085678100586, 0.45274314284324646, 7.243125915527344, -103.40960693359375},
            test_case{8.977622985839844, 0.6617381572723389, -6.434010982513428, -36.93961715698242},
            test_case{3.3447489738464355, 0.334273099899292, 5.010626792907715, -3.890927791595459},
            test_case{-8.08204174041748, 0.8953849077224731, 2.2224080562591553, 8.20623779296875},
            test_case{-1.163206696510315, 0.26193270087242126, -5.816899299621582, 7.301074028015137},
            test_case{7.729598522186279, 0.34303176403045654, 5.197448253631592, 6.264552593231201},
            test_case{3.9490699768066406, 0.26367634534835815, -5.014788627624512, 14.182146072387695},
            test_case{-3.4705426692962646, 0.35833755135536194, -8.288565635681152, 13.276227951049805},
            test_case{4.678563117980957, 0.014921886846423149, 2.3611345291137695, 1.0668786764144897},
            test_case{-5.597301006317139, 0.6260324716567993, 0.7393665909767151, 19.71807861328125},
            test_case{-8.368108749389648, 0.28673505783081055, 2.6905341148376465, -19.644155502319336},
            test_case{-6.802087783813477, 0.07672593742609024, -6.512517929077148, -0.20285895466804504},
            test_case{-3.1979963779449463, 0.6144924163818359, -5.036710262298584, -6.317578315734863},
            test_case{-0.6961369514465332, 0.1827998012304306, 3.6964597702026367, -4.8550567626953125},
            test_case{-4.671579360961914, 0.30830061435699463, -8.382567405700684, 7.411998748779297},
            test_case{6.31552791595459, 0.44206908345222473, 7.501471996307373, -8.138211250305176},
            test_case{-6.134112358093262, 0.15719829499721527, -1.426112413406372, -4.679737091064453},
            test_case{-7.410618305206299, 0.22357028722763062, 2.3678839206695557, -12.863222122192383},
            test_case{-8.166705131530762, 0.4748464524745941, -3.7378900051116943, -70.19567108154297},
            test_case{1.9713603258132935, 0.4768414795398712, -6.420742988586426, 144.48605346679688}
        );
        const float rtol = 100 * std::numeric_limits<float>::epsilon();
        const auto output = xsf::nrdtrisd(mean, p, x);
        const auto rel_error = xsf::extended_relative_error(output, expected);
        CAPTURE(mean, p, x, output, expected, rtol, rel_error);
        REQUIRE(rel_error <= rtol);
    }
}
