#include "../testing_utils.h"

#include <array>
#include <cmath>
#include <limits>

#include <xsf/cpu/stats.h>

TEST_CASE("generalized hyperbolic density", "[genhyperbolic][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.1/scipy/stats/tests/test_distributions.py#L962-L982
    constexpr std::array<double, 10> expected_pdf{
        2.94895678275316e-13, 1.75746848647696e-10, 9.48149804073045e-08, 4.17862521692026e-05, 0.0103947630463822,
        0.240864958986839,    0.162833527161649,    0.0374609592899472,   0.00634894847327781,  0.000941920705790324,
    };
    constexpr double p = 2.0;
    constexpr double a = 3.0;
    constexpr double b = 1.5;
    constexpr double loc = 0.5;
    constexpr double scale = 1.5;

    const std::vector<double> x_values = linspace(-10.0, 10.0, expected_pdf.size());
    for (int i = 0; i < expected_pdf.size(); ++i) {
        double x = x_values[i];
        double standardized_x = (x - loc) / scale;
        double pdf = xsf::cpu::genhyperbolic_pdf(standardized_x, p, a, b) / scale;
        CAPTURE(i, x, standardized_x, pdf, expected_pdf[i]);
        REQUIRE(xsf::extended_relative_error(pdf, expected_pdf[i]) <= 1e-13);
    }
}

TEST_CASE("generalized hyperbolic Student's t limit", "[genhyperbolic][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.1/scipy/stats/tests/test_distributions.py#L1081-L1097
    constexpr std::array<double, 10> lower{
        -31.820519750798752, -3.640296435003648, -2.9487520061496664, -2.7322159231076126, -2.6271579957950233,
        -2.565225109345393,  -2.524412431436311, -2.495498782411984,  -2.4739457843399073, -2.4572615423796655,
    };
    constexpr std::array<double, 10> upper{
        31.820140377530542, 3.640296435001152,  2.9487520061496824, 2.732215923107273, 2.627157995794383,
        2.5652251093455574, 2.5244124314366156, 2.4954987824114236, 2.473945784340125, 2.4572615423800115,
    };
    constexpr double alpha_epsilon = std::numeric_limits<float>::epsilon();

    const std::vector<double> degrees_of_freedom = linspace(1.0, 30.0, lower.size());
    for (int j = 0; j < degrees_of_freedom.size(); ++j) {
        double df = degrees_of_freedom[j];
        double p = -df / 2.0;
        double a = df * df * alpha_epsilon;
        double scale = std::sqrt(df);

        for (const double x : linspace(lower[j], upper[j], 50)) {
            double pdf = xsf::cpu::genhyperbolic_pdf(x / scale, p, a, 0.0) / scale;
            double expected_pdf = std::tgamma((df + 1.0) / 2.0) / (std::sqrt(df * M_PI) * std::tgamma(df / 2.0)) *
                                  std::pow(1.0 + x * x / df, -(df + 1.0) / 2.0);
            CAPTURE(j, df, x, pdf, expected_pdf);
            REQUIRE(xsf::extended_relative_error(pdf, expected_pdf) <= 1e-6);
        }
    }
}

TEST_CASE("generalized hyperbolic Cauchy limit", "[genhyperbolic][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.1/scipy/stats/tests/test_distributions.py#L1099-L1114
    constexpr double p = -0.5;
    constexpr double a = std::numeric_limits<float>::epsilon();
    constexpr double lower = -31.820519750798752;
    constexpr double upper = 31.820140377530542;

    for (const double x : linspace(lower, upper, 50)) {
        double pdf = xsf::cpu::genhyperbolic_pdf(x, p, a, 0.0);
        double expected_pdf = 1.0 / (M_PI * (1.0 + x * x));
        CAPTURE(x, pdf, expected_pdf);
        REQUIRE(xsf::extended_relative_error(pdf, expected_pdf) <= 1e-6);
    }
}

TEST_CASE("generalized hyperbolic Laplace limit", "[genhyperbolic][xsf_tests]") {
    // Mirrors https://github.com/scipy/scipy/blob/v1.18.1/scipy/stats/tests/test_distributions.py#L1116-L1135
    constexpr double scale = std::numeric_limits<float>::epsilon();

    for (const double loc : linspace(-10.0, 10.0, 10)) {
        for (const double x : linspace(-20.0, 20.0, 50)) {
            double pdf = xsf::cpu::genhyperbolic_pdf((x - loc) / scale, 1.0, scale, 0.0) / scale;
            double expected_pdf = 0.5 * std::exp(-std::abs(x - loc));
            CAPTURE(loc, x, pdf, expected_pdf);
            REQUIRE(xsf::extended_relative_error(pdf, expected_pdf) <= 1e-11);
        }
    }
}
