#pragma once

#ifndef __CUDACC__

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <math.h>
#include <tuple>
#include <type_traits>
#include <utility>

#endif

// Define math constants if they are not available
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#ifndef M_LOG2E
#define M_LOG2E 1.44269504088896340736
#endif

#ifndef M_LOG10E
#define M_LOG10E 0.434294481903251827651
#endif

#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4 0.785398163397448309616
#endif

#ifndef M_1_PI
#define M_1_PI 0.318309886183790671538
#endif

#ifndef M_2_PI
#define M_2_PI 0.636619772367581343076
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524401
#endif

#if defined(__CUDACC__) || defined(__HIPCC__)
#if defined(__CUDACC__)
#include <cuda/std/array>
#include <cuda/std/cmath>
#include <cuda/std/cstddef>
#include <cuda/std/cstdint>
#include <cuda/std/limits>
#include <cuda/std/tuple>
#include <cuda/std/type_traits>
#include <cuda/std/utility>
#endif
#define XSF_HOST_DEVICE __host__ __device__
#else
#define XSF_HOST_DEVICE
#endif

// Define target platform macros
#if defined(__CUDA_ARCH__)
#define XSF_TARGET_CUDA
#elif defined(__HIP_DEVICE_COMPILE__)
#define XSF_TARGET_HIP
#else
#define XSF_TARGET_CPU
#endif

#if defined(XSF_TARGET_CUDA)
#include <cuda_runtime.h>
#endif

namespace xsf::cxx {

#if defined(__CUDACC__)

using cuda::std::ptrdiff_t;
using cuda::std::size_t;
using cuda::std::uint64_t;

template <typename T, size_t N>
using array = cuda::std::array<T, N>;

// Must use thrust for complex types in order to support CuPy
template <typename T>
using complex = thrust::complex<T>;

template <typename T1, typename T2>
using pair = cuda::std::pair<T1, T2>;

template <typename... Types>
using tuple = cuda::std::tuple<Types...>;

template <typename T>
using is_floating_point = cuda::std::is_floating_point<T>;

template <typename T>
using is_integral = cuda::std::is_integral<T>;

template <typename T>
inline constexpr bool is_integral_v = cuda::std::is_integral_v<T>;

template <typename T>
using is_signed = cuda::std::is_signed<T>;

template <typename T>
inline constexpr bool is_signed_v = cuda::std::is_signed_v<T>;

template <typename T1, typename T2>
using is_same = cuda::std::is_same<T1, T2>;

template <typename T1, typename T2>
inline constexpr bool is_same_v = cuda::std::is_same_v<T1, T2>;

template <typename T>
using make_unsigned = cuda::std::make_unsigned<T>;

template <typename T>
using make_unsigned_t = cuda::std::make_unsigned_t<T>;

template <bool Cond, typename T = void>
using enable_if = cuda::std::enable_if<Cond, T>;

template <typename T>
using decay = cuda::std::decay<T>;

template <typename F>
struct invoke_result {
    using type = decltype(cuda::std::declval<F>()());
};

template <typename F>
using invoke_result_t = typename invoke_result<F>::type;

#else

using std::ptrdiff_t;
using std::size_t;
using std::uint64_t;

template <typename T, std::size_t N>
using array = std::array<T, N>;

template <typename T>
using complex = std::complex<T>;

template <typename T1, typename T2>
using pair = std::pair<T1, T2>;

template <typename... Types>
using tuple = std::tuple<Types...>;

// Type traits
template <typename T>
using is_floating_point = std::is_floating_point<T>;

template <typename T>
using is_integral = std::is_integral<T>;

template <typename T>
inline constexpr bool is_integral_v = std::is_integral_v<T>;

template <typename T>
using is_signed = std::is_signed<T>;

template <typename T>
inline constexpr bool is_signed_v = std::is_signed_v<T>;

template <typename T1, typename T2>
using is_same = std::is_same<T1, T2>;

template <typename T1, typename T2>
inline constexpr bool is_same_v = std::is_same_v<T1, T2>;

template <typename T>
using make_unsigned = std::make_unsigned<T>;

template <typename T>
using make_unsigned_t = std::make_unsigned_t<T>;

template <bool Cond, typename T = void>
using enable_if = std::enable_if<Cond, T>;

template <typename T>
using decay = std::decay<T>;

template <typename F>
using invoke_result = std::invoke_result<F>;

template <typename F>
using invoke_result_t = std::invoke_result_t<F>;

#endif

#if defined(__CUDACC__)
template <typename T>
using numeric_limits = cuda::std::numeric_limits<T>;
#else
template <typename T>
using numeric_limits = std::numeric_limits<T>;
#endif

XSF_HOST_DEVICE inline double abs(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::abs(num);
#else
    return std::abs(num);
#endif
}

XSF_HOST_DEVICE inline double fabs(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::fabs(num);
#else
    return std::fabs(num);
#endif
}

XSF_HOST_DEVICE inline double exp(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::exp(num);
#else
    return std::exp(num);
#endif
}

XSF_HOST_DEVICE inline double expm1(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::expm1(num);
#else
    return std::expm1(num);
#endif
}

XSF_HOST_DEVICE inline double log(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::log(num);
#else
    return std::log(num);
#endif
}

XSF_HOST_DEVICE inline double sqrt(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::sqrt(num);
#else
    return std::sqrt(num);
#endif
}

XSF_HOST_DEVICE inline bool isinf(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::isinf(num);
#else
    return std::isinf(num);
#endif
}

XSF_HOST_DEVICE inline bool isnan(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::isnan(num);
#else
    return std::isnan(num);
#endif
}

XSF_HOST_DEVICE inline bool isfinite(double num) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::isfinite(num);
#else
    return std::isfinite(num);
#endif
}

XSF_HOST_DEVICE inline double pow(double x, double y) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::pow(x, y);
#else
    return std::pow(x, y);
#endif
}

XSF_HOST_DEVICE inline double sin(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::sin(x);
#else
    return std::sin(x);
#endif
}

XSF_HOST_DEVICE inline double cos(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::cos(x);
#else
    return std::cos(x);
#endif
}

XSF_HOST_DEVICE inline double tan(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::tan(x);
#else
    return std::tan(x);
#endif
}

XSF_HOST_DEVICE inline double atan(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::atan(x);
#else
    return std::atan(x);
#endif
}

XSF_HOST_DEVICE inline double asin(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::asin(x);
#else
    return std::asin(x);
#endif
}

XSF_HOST_DEVICE inline double acos(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::acos(x);
#else
    return std::acos(x);
#endif
}

XSF_HOST_DEVICE inline double sinh(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::sinh(x);
#else
    return std::sinh(x);
#endif
}

XSF_HOST_DEVICE inline double cosh(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::cosh(x);
#else
    return std::cosh(x);
#endif
}

XSF_HOST_DEVICE inline double asinh(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::asinh(x);
#else
    return std::asinh(x);
#endif
}

XSF_HOST_DEVICE inline double tanh(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::tanh(x);
#else
    return std::tanh(x);
#endif
}

XSF_HOST_DEVICE inline double atanh(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::atanh(x);
#else
    return std::atanh(x);
#endif
}

XSF_HOST_DEVICE inline bool signbit(double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::signbit(x);
#else
    return std::signbit(x);
#endif
}

XSF_HOST_DEVICE inline double hypot(double x, double y) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::hypot(x, y);
#else
    return std::hypot(x, y);
#endif
}

XSF_HOST_DEVICE inline double atan2(double y, double x) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::atan2(y, x);
#else
    return std::atan2(y, x);
#endif
}

// TODO: Check if separating for NVRTC compilation is necessary
XSF_HOST_DEVICE inline double ceil(double x) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::ceil(x);
#else
    return cuda::std::ceil(x);
#endif
#else
    return std::ceil(x);
#endif
}

XSF_HOST_DEVICE inline double floor(double x) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::floor(x);
#else
    return cuda::std::floor(x);
#endif
#else
    return std::floor(x);
#endif
}

XSF_HOST_DEVICE inline double round(double x) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::round(x);
#else
    return cuda::std::round(x);
#endif
#else
    return std::round(x);
#endif
}

XSF_HOST_DEVICE inline double trunc(double x) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::trunc(x);
#else
    return cuda::std::trunc(x);
#endif
#else
    return std::trunc(x);
#endif
}

XSF_HOST_DEVICE inline double fma(double x, double y, double z) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::fma(x, y, z);
#else
    return cuda::std::fma(x, y, z);
#endif
#else
    return std::fma(x, y, z);
#endif
}

XSF_HOST_DEVICE inline double copysign(double x, double y) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::copysign(x, y);
#else
    return cuda::std::copysign(x, y);
#endif
#else
    return std::copysign(x, y);
#endif
}

XSF_HOST_DEVICE inline double modf(double value, double *iptr) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::modf(value, iptr);
#else
    return cuda::std::modf(value, iptr);
#endif
#else
    return std::modf(value, iptr);
#endif
}

XSF_HOST_DEVICE inline double fmax(double x, double y) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::fmax(x, y);
#else
    return cuda::std::fmax(x, y);
#endif
#else
    return std::fmax(x, y);
#endif
}

XSF_HOST_DEVICE inline double fmin(double x, double y) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::fmin(x, y);
#else
    return cuda::std::fmin(x, y);
#endif
#else
    return std::fmin(x, y);
#endif
}

XSF_HOST_DEVICE inline double log10(double num) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::log10(num);
#else
    return cuda::std::log10(num);
#endif
#else
    return std::log10(num);
#endif
}

XSF_HOST_DEVICE inline double log1p(double num) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::log1p(num);
#else
    return cuda::std::log1p(num);
#endif
#else
    return std::log1p(num);
#endif
}

XSF_HOST_DEVICE inline double frexp(double num, int *exp) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::frexp(num, exp);
#else
    return cuda::std::frexp(num, exp);
#endif
#else
    return std::frexp(num, exp);
#endif
}

XSF_HOST_DEVICE inline double ldexp(double num, int exp) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::ldexp(num, exp);
#else
    return cuda::std::ldexp(num, exp);
#endif
#else
    return std::ldexp(num, exp);
#endif
}

XSF_HOST_DEVICE inline double fmod(double x, double y) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::fmod(x, y);
#else
    return cuda::std::fmod(x, y);
#endif
#else
    return std::fmod(x, y);
#endif
}

XSF_HOST_DEVICE inline double nextafter(double from, double to) {
#if defined(XSF_TARGET_CUDA)
#if defined(__CUDACC_RTC__)
    return ::nextafter(from, to);
#else
    return cuda::std::nextafter(from, to);
#endif
#else
    return std::nextafter(from, to);
#endif
}

template <typename T>
XSF_HOST_DEVICE void swap(T &a, T &b) {
#if defined(XSF_TARGET_CUDA)
    return cuda::std::swap(a, b);
#else
    return std::swap(a, b);
#endif
}

// Reimplement std::min, std::max, std::clamp until they are available in CuPy
template <typename T>
XSF_HOST_DEVICE constexpr const T &min(const T &a, const T &b) {
    return a < b ? a : b;
}

template <typename T>
XSF_HOST_DEVICE constexpr const T &max(const T &a, const T &b) {
    return a < b ? b : a;
}

template <typename T>
XSF_HOST_DEVICE constexpr const T &clamp(const T &v, const T &lo, const T &hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

template <typename T>
XSF_HOST_DEVICE T real(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::real(z);
#else
    return std::real(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE T imag(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::imag(z);
#else
    return std::imag(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE T abs(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::abs(z);
#else
    return std::abs(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> exp(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::exp(z);
#else
    return std::exp(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> log(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::log(z);
#else
    return std::log(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE T norm(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::norm(z);
#else
    return std::norm(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> sqrt(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::sqrt(z);
#else
    return std::sqrt(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> conj(const complex<T> &z) {
#if defined(XSF_TARGET_CUDA)
    return thrust::conj(z);
#else
    return std::conj(z);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> pow(const complex<T> &x, int y) {
#if defined(XSF_TARGET_CUDA)
    return thrust::pow(x, y);
#else
    return std::pow(x, y);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> pow(const complex<T> &x, const complex<T> &y) {
#if defined(XSF_TARGET_CUDA)
    return thrust::pow(x, y);
#else
    return std::pow(x, y);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> pow(const complex<T> &x, const T &y) {
#if defined(XSF_TARGET_CUDA)
    return thrust::pow(x, y);
#else
    return std::pow(x, y);
#endif
}

template <typename T>
XSF_HOST_DEVICE complex<T> pow(const T &x, const complex<T> &y) {
#if defined(XSF_TARGET_CUDA)
    return thrust::pow(x, y);
#else
    return std::pow(x, y);
#endif
}

} // namespace xsf::cxx

#ifdef __CUDACC__

// Fallback to global namespace for functions unsupported on NVRTC Jit
#ifdef _LIBCUDACXX_COMPILER_NVRTC
#include <cuda_runtime.h>
#endif

#endif

#ifdef __CUDACC__
#define XSF_ASSERT(a)
#else
#ifdef DEBUG
#define XSF_ASSERT(a) assert(a)
#else
#define XSF_ASSERT(a)
#endif
#endif

#ifndef __CUDACC__

namespace xsf {

// basic
using std::abs;

// exponential
using std::exp;

// power
using std::sqrt;

// trigonometric
using std::cos;
using std::sin;

// floating-point manipulation
using std::copysign;

// classification and comparison
using std::isfinite;
using std::isinf;
using std::isnan;
using std::signbit;

// complex
using std::imag;
using std::real;

template <typename T>
struct remove_complex {
    using type = T;
};

template <typename T>
struct remove_complex<std::complex<T>> {
    using type = T;
};

template <typename T>
using remove_complex_t = typename remove_complex<T>::type;

template <typename T>
struct complex_type {
    using type = std::complex<T>;
};

template <typename T>
using complex_type_t = typename complex_type<T>::type;

template <typename T>
using complex = complex_type_t<T>;

} // namespace xsf

#endif // __CUDACC__
