#pragma once

#include "../config.h"
#include "../ellint_carlson/ellint_carlson.hh"
#include "../error.h"

namespace xsf {
namespace cpu {
    namespace detail {

        constexpr double ellint_carlson_rerr = 5e-16;

        template <typename T>
        inline T elliprc(T x, T y) {
            T result;
            const auto status = ellint_carlson::rc(x, y, ellint_carlson_rerr, result);
            set_error("elliprc", static_cast<sf_error_t>(status), nullptr);
            return result;
        }

        template <typename T>
        inline T elliprd(T x, T y, T z) {
            T result;
            const auto status = ellint_carlson::rd(x, y, z, ellint_carlson_rerr, result);
            set_error("elliprd", static_cast<sf_error_t>(status), nullptr);
            return result;
        }

        template <typename T>
        inline T elliprf(T x, T y, T z) {
            T result;
            const auto status = ellint_carlson::rf(x, y, z, ellint_carlson_rerr, result);
            set_error("elliprf", static_cast<sf_error_t>(status), nullptr);
            return result;
        }

        template <typename T>
        inline T elliprg(T x, T y, T z) {
            T result;
            const auto status = ellint_carlson::rg(x, y, z, ellint_carlson_rerr, result);
            set_error("elliprg", static_cast<sf_error_t>(status), nullptr);
            return result;
        }

        template <typename T>
        inline T elliprj(T x, T y, T z, T p) {
            T result;
            const auto status = ellint_carlson::rj(x, y, z, p, ellint_carlson_rerr, result);
            set_error("elliprj", static_cast<sf_error_t>(status), nullptr);
            return result;
        }

    } // namespace detail

    inline double elliprc(double x, double y) { return detail::elliprc(x, y); }

    inline float elliprc(float x, float y) {
        return static_cast<float>(elliprc(static_cast<double>(x), static_cast<double>(y)));
    }

    inline std::complex<double> elliprc(std::complex<double> x, std::complex<double> y) {
        return detail::elliprc(x, y);
    }

    inline std::complex<float> elliprc(std::complex<float> x, std::complex<float> y) {
        return static_cast<std::complex<float>>(
            elliprc(static_cast<std::complex<double>>(x), static_cast<std::complex<double>>(y))
        );
    }

    inline double elliprd(double x, double y, double z) { return detail::elliprd(x, y, z); }

    inline float elliprd(float x, float y, float z) {
        return static_cast<float>(elliprd(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)));
    }

    inline std::complex<double> elliprd(std::complex<double> x, std::complex<double> y, std::complex<double> z) {
        return detail::elliprd(x, y, z);
    }

    inline std::complex<float> elliprd(std::complex<float> x, std::complex<float> y, std::complex<float> z) {
        return static_cast<std::complex<float>>(elliprd(
            static_cast<std::complex<double>>(x), static_cast<std::complex<double>>(y),
            static_cast<std::complex<double>>(z)
        ));
    }

    inline double elliprf(double x, double y, double z) { return detail::elliprf(x, y, z); }

    inline float elliprf(float x, float y, float z) {
        return static_cast<float>(elliprf(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)));
    }

    inline std::complex<double> elliprf(std::complex<double> x, std::complex<double> y, std::complex<double> z) {
        return detail::elliprf(x, y, z);
    }

    inline std::complex<float> elliprf(std::complex<float> x, std::complex<float> y, std::complex<float> z) {
        return static_cast<std::complex<float>>(elliprf(
            static_cast<std::complex<double>>(x), static_cast<std::complex<double>>(y),
            static_cast<std::complex<double>>(z)
        ));
    }

    inline double elliprg(double x, double y, double z) { return detail::elliprg(x, y, z); }

    inline float elliprg(float x, float y, float z) {
        return static_cast<float>(elliprg(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)));
    }

    inline std::complex<double> elliprg(std::complex<double> x, std::complex<double> y, std::complex<double> z) {
        return detail::elliprg(x, y, z);
    }

    inline std::complex<float> elliprg(std::complex<float> x, std::complex<float> y, std::complex<float> z) {
        return static_cast<std::complex<float>>(elliprg(
            static_cast<std::complex<double>>(x), static_cast<std::complex<double>>(y),
            static_cast<std::complex<double>>(z)
        ));
    }

    inline double elliprj(double x, double y, double z, double p) { return detail::elliprj(x, y, z, p); }

    inline float elliprj(float x, float y, float z, float p) {
        return static_cast<float>(
            elliprj(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z), static_cast<double>(p))
        );
    }

    inline std::complex<double>
    elliprj(std::complex<double> x, std::complex<double> y, std::complex<double> z, std::complex<double> p) {
        return detail::elliprj(x, y, z, p);
    }

    inline std::complex<float>
    elliprj(std::complex<float> x, std::complex<float> y, std::complex<float> z, std::complex<float> p) {
        return static_cast<std::complex<float>>(elliprj(
            static_cast<std::complex<double>>(x), static_cast<std::complex<double>>(y),
            static_cast<std::complex<double>>(z), static_cast<std::complex<double>>(p)
        ));
    }

} // namespace cpu
} // namespace xsf
