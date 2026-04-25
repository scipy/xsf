#pragma once

#include <cmath>
#include <limits>

#include "cephes/ndtri.h"

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

    // Returns the std deviation of a normal distribution with given mean
    // such that the CDF at x equals p.
    inline double nrdtrisd(double mean, double p, double x) {
        if (std::isnan(mean) || std::isnan(p) || std::isnan(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (p <= 0 || p >= 1) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return (x - mean) / cephes::ndtri(p);
    }

} // namespace detail

inline double nrdtrimn(double p, double std, double x) { return detail::nrdtrimn(p, std, x); }
inline double nrdtrisd(double mean, double p, double x) { return detail::nrdtrisd(mean, p, x);}

inline float nrdtrimn(float p, float std, float x) {
    return static_cast<float>(detail::nrdtrimn(static_cast<double>(p), static_cast<double>(std), static_cast<double>(x)));
}
inline float nrdtrisd(float mean, float p, float x) {
    return static_cast<float>(detail::nrdtrisd(static_cast<double>(mean), static_cast<double>(p), static_cast<double>(x)));
} 

} // namespace xsf