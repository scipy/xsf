#pragma once

#include "config.h"

namespace xsf {
namespace numbers {

    template <typename T>
    inline constexpr typename cxx::enable_if<cxx::is_floating_point<T>::value, cxx::complex<T>>::type i_v =
        cxx::complex<T>(0.0, 1.0);

} // namespace numbers
} // namespace xsf
