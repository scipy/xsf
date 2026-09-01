#pragma once

#include "config.h"

namespace xsf {
namespace numbers {

    template <typename T>
    std::complex<T> i_v;

    template <>
    std::complex<float> i_v<float> = std::literals::complex_literals::operator""if(1.0L);

    template <>
    std::complex<double> i_v<double> = std::literals::complex_literals::operator""i(1.0L);

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type e_v =
        T(2.718281828459045235360287471352662l);
    inline constexpr double e = e_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type log2e_v =
        T(1.442695040888963407359924681001892l);
    inline constexpr double log2e = log2e_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type log10e_v =
        T(0.434294481903251827651128918916605l);
    inline constexpr double log10e = log10e_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type pi_v =
        T(3.141592653589793238462643383279502l);
    inline constexpr double pi = pi_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type inv_pi_v =
        T(0.318309886183790671537767526745028l);
    inline constexpr double inv_pi = inv_pi_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type inv_sqrtpi_v =
        T(0.564189583547756286948079451560772l);
    inline constexpr double inv_sqrtpi = inv_sqrtpi_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type ln2_v =
        T(0.693147180559945309417232121458176l);
    inline constexpr double ln2 = ln2_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type ln10_v =
        T(2.302585092994045684017991454684364l);
    inline constexpr double ln10 = ln10_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type sqrt2_v =
        T(1.414213562373095048801688724209698l);
    inline constexpr double sqrt2 = sqrt2_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type sqrt3_v =
        T(1.732050807568877293527446341505872l);
    inline constexpr double sqrt3 = sqrt3_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type inv_sqrt3_v =
        T(0.577350269189625764509148780501957l);
    inline constexpr double inv_sqrt3 = inv_sqrt3_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type egamma_v =
        T(0.577215664901532860606512090082402l);
    inline constexpr double egamma = egamma_v<double>;

    template <typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type phi_v =
        T(1.618033988749894848204586834365638l);
    inline constexpr double phi = phi_v<double>;


} // namespace numbers
} // namespace xsf
