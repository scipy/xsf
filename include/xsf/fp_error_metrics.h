#include "config.h"

namespace xsf {

template <typename T>
XSF_HOST_DEVICE typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type
extended_absolute_error(T actual, T desired) {
    if (actual == desired || (cxx::isnan(actual) && cxx::isnan(desired))) {
        return T(0);
    }
    if (cxx::isnan(desired) || cxx::isnan(actual)) {
        /* If expected nan but got non-NaN or expected non-NaN but got NaN
         * we consider this to be an infinite error. */
        return cxx::numeric_limits<T>::infinity();
    }
    if (cxx::isinf(actual)) {
        /* We don't want to penalize early overflow too harshly, so instead
         * compare with the mythical value nextafter(max_float). */
        T sgn = cxx::copysign(1.0, actual);
        T max_float = cxx::numeric_limits<T>::max();
        // max_float * 2**-(mantissa_bits + 1) = ulp(max_float)
        T ulp = cxx::pow(2, -cxx::numeric_limits<T>::digits) * max_float;
        return cxx::abs((sgn * cxx::numeric_limits<T>::max() - desired) + sgn * ulp);
    }
    if (cxx::isinf(desired)) {
        T sgn = cxx::copysign(1.0, desired);
        T max_float = cxx::numeric_limits<T>::max();
        // max_float * 2**-(mantissa_bits + 1) = ulp(max_float)
        T ulp = cxx::pow(2, -cxx::numeric_limits<T>::digits) * max_float;
        return cxx::abs((sgn * cxx::numeric_limits<T>::max() - actual) + sgn * ulp);
    }
    return cxx::abs(actual - desired);
}

template <typename T>
XSF_HOST_DEVICE T extended_absolute_error(cxx::complex<T> actual, cxx::complex<T> desired) {
    return cxx::hypot(
        extended_absolute_error(actual.real(), desired.real()), extended_absolute_error(actual.imag(), desired.imag())
    );
}

template <typename T>
XSF_HOST_DEVICE typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type
extended_relative_error(T actual, T desired) {
    T abs_error = extended_absolute_error(actual, desired);
    T abs_desired = cxx::abs(desired);
    if (desired == 0.0) {
        /* If the desired result is 0.0, normalize by smallest subnormal instead
         * of zero. */
        abs_desired = cxx::numeric_limits<T>::denorm_min();
    } else if (cxx::isinf(desired)) {
        abs_desired = cxx::numeric_limits<T>::max();
    } else if (cxx::isnan(desired)) {
        /* This ensures extended_relative_error(nan, nan) = 0 but
         * extended_relative_error(x0, x1) is infinite if one but not both of
         * x0 and x1 equals NaN */
        abs_desired = T(1);
    }
    return abs_error / abs_desired;
}

template <typename T>
XSF_HOST_DEVICE T extended_relative_error(cxx::complex<T> actual, cxx::complex<T> desired) {
    T abs_error = extended_absolute_error(actual, desired);

    if (desired.real() == 0.0) {
        desired.real(cxx::copysign(cxx::numeric_limits<T>::denorm_min(), desired.real()));
    } else if (cxx::isinf(desired.real())) {
        desired.real(cxx::copysign(cxx::numeric_limits<T>::max(), desired.real()));
    } else if (cxx::isnan(desired.real())) {
        /* In this case, the value used for desired doesn't matter. If desired.real() is NaN
         * but actual.real() isn't NaN, then the extended_absolute_error will be inf already
         * anyway. */
        desired.real(1.0);
    }

    if (desired.imag() == 0.0) {
        desired.imag(cxx::copysign(cxx::numeric_limits<T>::denorm_min(), desired.imag()));
    } else if (cxx::isinf(desired.imag())) {
        desired.imag(cxx::copysign(cxx::numeric_limits<T>::max(), desired.imag()));
    } else if (cxx::isnan(desired.imag())) {
        /* In this case, the value used for desired doesn't matter. If desired.imag() is NaN
         * but actual.imag() isn't NaN, then the extended_absolute_error will be inf already
         * anyway. */
        desired.imag(1.0);
    }

    if (!(cxx::isinf(desired.real()) || cxx::isinf(desired.imag())) && cxx::isinf(cxx::abs(desired))) {
        /* Rescale to avoid overflow */
        return (abs_error / 2.0) / (cxx::abs(desired / 2.0));
    }

    return abs_error / abs(desired);
}

} // namespace xsf
