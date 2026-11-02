#pragma once

#include "config.h"
#include "mdspan.h"
#include "relaxed_layout.h"

namespace xsf {

template <typename T, int ndim, bool is_c_contiguous, bool index_32_bits, int core_ndim>
__device__ inline auto as_mdspan(const CArray<T, ndim, is_c_contiguous, index_32_bits, core_ndim> &arr) {
    cxx::array<cxx::ptrdiff_t, ndim> exts;
    cxx::array<cxx::ptrdiff_t, ndim> strs;

    for (int i = 0; i < ndim; ++i) {
        exts[i] = static_cast<cxx::ptrdiff_t>(arr.shape_[i]);
        strs[i] = static_cast<cxx::ptrdiff_t>(arr.strides_[i]) / static_cast<cxx::ptrdiff_t>(sizeof(T));
    }

    using Extents = cxx::dextents<cxx::ptrdiff_t, ndim>;
    using Mapping = relaxed_layout::mapping<Extents>;

    Mapping mapping{Extents(exts), strs};

    auto *data = arr.data_;
    if (mapping.offset() != 0) {
        data -= mapping.offset();
    }

    return cxx::mdspan<T, Extents, relaxed_layout>(data, mapping);
}

} // namespace xsf
