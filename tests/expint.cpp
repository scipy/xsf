#include <limits>
#include <complex>
#include <gtest/gtest.h>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

#include <xsf/expint.h>
#include <xtest.hpp>
namespace xtest = xsf::test;


TEST(ExpInt_exp1, SpecialValue) {
    EXPECT_TRUE(std::isnan(xsf::exp1(xtest::NaN64)));
    // x < 0:   NaN
    EXPECT_TRUE(std::isnan(xsf::exp1(-1.0)));
    EXPECT_TRUE(std::isnan(xsf::exp1(-1000.0)));
    EXPECT_TRUE(std::isnan(xsf::exp1(-xtest::Inf64)));
    // x == 0:  +Inf
    EXPECT_TRUE(std::isinf(xsf::exp1(0.0)));
    // x == +Inf:   0
    EXPECT_EQ(0.0, xsf::exp1(xtest::Inf64));

    // Test float exp1(float x)
    EXPECT_TRUE(std::isnan(xsf::exp1(xtest::NaN32)));
    EXPECT_TRUE(std::isinf(xsf::exp1((float)0.0)));
}
