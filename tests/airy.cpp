#include <limits>
#include <complex>
#include <gtest/gtest.h>
#include <xsf/airy.h>

TEST(Airy, BasicAssertions) {
  const double nan64 = std::numeric_limits<double>::quiet_NaN();
  const std::complex<double> nan64c(std::nan(""), std::nan(""));
  double x, ai, aip, bi, bip;
  
  x = 0.0;
  xsf::airy(x, ai, aip, bi, bip);
  EXPECT_NE(ai, nan64);
  EXPECT_NE(aip, nan64);
  EXPECT_NE(bi, nan64);
  EXPECT_NE(bip, nan64);

  xsf::airye(x, ai, aip, bi, bip);
  EXPECT_NE(ai, nan64);
  EXPECT_NE(aip, nan64);
  EXPECT_NE(bi, nan64);
  EXPECT_NE(bip, nan64);

  double apt, bpt, ant, bnt;
  x = 1.0;
  xsf::itairy(x, apt, bpt, ant, bnt);
  EXPECT_NE(apt, nan64);
  EXPECT_NE(bpt, nan64);
  EXPECT_NE(ant, nan64);
  EXPECT_NE(bnt, nan64);

  // std::complex
  std::complex<double> z, cai, caip, cbi, cbip;
  z = 0.0;
  xsf::airy(z, cai, caip, cbi, cbip);
  EXPECT_NE(cai, nan64c);
  EXPECT_NE(caip, nan64c);
  EXPECT_NE(cbi, nan64c);
  EXPECT_NE(cbip, nan64c);

  xsf::airye(z, cai, caip, cbi, cbip);
  EXPECT_NE(cai, nan64c);
  EXPECT_NE(caip, nan64c);
  EXPECT_NE(cbi, nan64c);
  EXPECT_NE(cbip, nan64c);
}
