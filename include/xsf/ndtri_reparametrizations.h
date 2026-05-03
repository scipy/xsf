#pragma once

#include <cmath>
#include <limits>

#include "xsf/cephes/ndtri.h"

namespace xsf {
namespace detail {

    // Returns the mean of a normal distribution with standard deviation std
    // such that the CDF at x equals p.
    inline double nrdtrimn(double p, double std, double x) {
        if (std::isnan(std) || std <= 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (std::isnan(p) || p <= 0 || p >= 1) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (std::isnan(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return x - std * cephes::ndtri(p);
    }

} // namespace detail

inline double nrdtrimn(double p, double std, double x) { return detail::nrdtrimn(p, std, x); }

inline float nrdtrimn(float p, float std, float x) {
    return static_cast<float>(detail::nrdtrimn(static_cast<double>(p), static_cast<double>(std), static_cast<double>(x))
    );
}

} // namespace xsf
