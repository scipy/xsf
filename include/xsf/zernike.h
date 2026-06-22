#pragma once

#include "config.h"
#include "error.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace xsf {

namespace detail {

enum class zernike_arg_status {
    ok,
    domain,
    nan
};

inline zernike_arg_status validate_zernike_radial_args(
    std::ptrdiff_t n, std::ptrdiff_t m, double rho, std::ptrdiff_t &m_abs
) {
    if (std::isnan(rho)) {
        return zernike_arg_status::nan;
    }

    if (n < 0 || m == std::numeric_limits<std::ptrdiff_t>::min()) {
        return zernike_arg_status::domain;
    }

    m_abs = m < 0 ? -m : m;

    if (m_abs > n || ((n - m_abs) & 1) != 0 || rho < 0.0 || rho > 1.0) {
        return zernike_arg_status::domain;
    }

    return zernike_arg_status::ok;
}

inline double zernike_radial_domain_error(const char *func_name) {
    set_error(func_name, SF_ERROR_DOMAIN, nullptr);
    return std::numeric_limits<double>::quiet_NaN();
}

} // namespace detail

inline double eval_zernike_radial_barmak(std::ptrdiff_t n, std::ptrdiff_t m, double rho) {
    std::ptrdiff_t m_abs = 0;
    const auto status = detail::validate_zernike_radial_args(n, m, rho, m_abs);

    if (status == detail::zernike_arg_status::nan) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (status == detail::zernike_arg_status::domain) {
        return detail::zernike_radial_domain_error("eval_zernike_radial");
    }

    try {
        const std::size_t width = static_cast<std::size_t>(n) + 3;
        std::vector<double> prev2(width, 0.0);
        std::vector<double> prev1(width, 0.0);
        std::vector<double> current(width, 0.0);

        prev1[0] = 1.0;

        if (n == 0) {
            return prev1[static_cast<std::size_t>(m_abs)];
        }

        for (std::ptrdiff_t nn = 1; nn <= n; ++nn) {
            std::fill(current.begin(), current.end(), 0.0);

            for (std::ptrdiff_t mm = 0; mm <= nn; ++mm) {
                if (((nn - mm) & 1) != 0) {
                    continue;
                }

                double left = 0.0;
                double right = 0.0;
                double grand = 0.0;

                const std::ptrdiff_t m_left = mm == 0 ? 1 : mm - 1;

                if (m_left <= nn - 1) {
                    left = prev1[static_cast<std::size_t>(m_left)];
                }
                if (mm + 1 <= nn - 1) {
                    right = prev1[static_cast<std::size_t>(mm + 1)];
                }
                if (nn >= 2) {
                    grand = prev2[static_cast<std::size_t>(mm)];
                }

                current[static_cast<std::size_t>(mm)] = rho * (left + right) - grand;
            }

            prev2.swap(prev1);
            prev1.swap(current);
        }

        return prev1[static_cast<std::size_t>(m_abs)];
    } catch (...) {
        set_error("eval_zernike_radial", SF_ERROR_MEMORY, nullptr);
        return std::numeric_limits<double>::quiet_NaN();
    }
}

inline double eval_zernike_radial_kintner(std::ptrdiff_t n, std::ptrdiff_t m, double rho) {
    std::ptrdiff_t m_abs = 0;
    const auto status = detail::validate_zernike_radial_args(n, m, rho, m_abs);

    if (status == detail::zernike_arg_status::nan) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (status == detail::zernike_arg_status::domain) {
        return detail::zernike_radial_domain_error("eval_zernike_radial_kintner");
    }

    const double rho2 = rho * rho;

    if (n == m_abs) {
        return std::pow(rho, static_cast<double>(m_abs));
    }

    double r_nm4 = std::pow(rho, static_cast<double>(m_abs));
    double r_nm2 = ((static_cast<double>(m_abs) + 2.0) * rho2 - (static_cast<double>(m_abs) + 1.0)) * r_nm4;

    if (n == m_abs + 2) {
        return r_nm2;
    }

    for (std::ptrdiff_t current = m_abs + 4; current <= n; current += 2) {
        const double cd = static_cast<double>(current);
        const double md = static_cast<double>(m_abs);

        const double a =
            2.0 * (cd - 1.0) * (2.0 * cd * (cd - 2.0) * rho2 - md * md - cd * (cd - 2.0));
        const double b = cd * (cd + md - 2.0) * (cd - md - 2.0);
        const double denom = (cd + md) * (cd - md) * (cd - 2.0);

        const double r_n = (a * r_nm2 - b * r_nm4) / denom;

        r_nm4 = r_nm2;
        r_nm2 = r_n;
    }

    return r_nm2;
}

inline double eval_zernike_radial(std::ptrdiff_t n, std::ptrdiff_t m, double rho) {
    return eval_zernike_radial_barmak(n, m, rho);
}

} // namespace xsf
