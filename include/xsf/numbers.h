#pragma once

#include "config.h"

namespace xsf {
namespace numbers {

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, cxx::complex<T>>::type i_v =
        cxx::complex<T>(0.0, 1.0);

    // These constants match llvm's libcxx <numbers> header as coded in
    // this (open, at the time of writing) pull request: https://github.com/llvm/llvm-project/pull/222830
    // at https://github.com/llvm/llvm-project/pull/222830/changes/14751b686988f11205327f0a12836dc1aae30b8a
    // which in turn states that it matches GCC's <numbers>

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type e_v =
        2.718281828459045235360287471352662498L;
    inline constexpr double e = e_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type log2e_v =
        1.442695040888963407359924681001892137L;
    inline constexpr double log2e = log2e_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type log10e_v =
        0.434294481903251827651128918916605082L;
    inline constexpr double log10e = log10e_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type pi_v =
        3.141592653589793238462643383279502884L;
    inline constexpr double pi = pi_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type inv_pi_v =
        0.318309886183790671537767526745028724L;
    inline constexpr double inv_pi = inv_pi_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type inv_sqrtpi_v =
        0.564189583547756286948079451560772586L;
    inline constexpr double inv_sqrtpi = inv_sqrtpi_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type ln2_v =
        0.693147180559945309417232121458176568L;
    inline constexpr double ln2 = ln2_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type ln10_v =
        2.302585092994045684017991454684364208L;
    inline constexpr double ln10 = ln10_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type sqrt2_v =
        1.414213562373095048801688724209698079L;
    inline constexpr double sqrt2 = sqrt2_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type sqrt3_v =
        1.732050807568877293527446341505872367L;
    inline constexpr double sqrt3 = sqrt3_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type inv_sqrt3_v =
        0.577350269189625764509148780501957456L;
    inline constexpr double inv_sqrt3 = inv_sqrt3_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type egamma_v =
        0.577215664901532860606512090082402431L;
    inline constexpr double egamma = egamma_v<double>;

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, T>::type phi_v =
        1.618033988749894848204586834365638118L;
    inline constexpr double phi = phi_v<double>;

} // namespace numbers
} // namespace xsf
