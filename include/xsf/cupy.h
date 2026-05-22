#pragma once

#include "config.h"
#include "third_party/kokkos/mdspan.hpp"

namespace xsf {

template <typename T, int ndim, int is_c_contiguous, int index_32_bits, int core_ndim>
__device__ inline auto as_mdspan(const CArray<T, ndim, is_c_contiguous, index_32_bits, core_ndim>& arr) {
    std::array<std::ptrdiff_t, ndim> exts;
    std::array<std::ptrdiff_t, ndim> strs;

    for (int i = 0; i < ndim; ++i) {
        exts[i] = static_cast<std::ptrdiff_t>(arr.shape_[i]);
        strs[i] = static_cast<std::ptrdiff_t>(arr.strides_[i]) / sizeof(T);
    }

    using Extents = std::dextents<std::ptrdiff_t, ndim>;
    using Mapping = std::layout_stride::mapping<Extents>;

    return std::mdspan<T, Extents, std::layout_stride>(
        arr.data_, Mapping(Extents(exts), strs)
    );
}

} // namespace xsf
