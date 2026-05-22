#pragma once

#include <cuda/std/array>
#include <cuda/std/cstddef>
#include <cuda/std/mdspan>

namespace xsf {

template <typename T, int ndim, bool is_c_contiguous, bool index_32_bits, int core_ndim>
__device__ inline auto as_mdspan(CArray<T, ndim, is_c_contiguous, index_32_bits, core_ndim> arr) {
    cuda::std::array<cuda::std::ptrdiff_t, ndim> exts;
    cuda::std::array<cuda::std::ptrdiff_t, ndim> strs;

    for (int i = 0; i < ndim; ++i) {
        exts[i] = static_cast<cuda::std::ptrdiff_t>(arr.shape_[i]);
        strs[i] = static_cast<cuda::std::ptrdiff_t>(arr.strides_[i]) / sizeof(T);
    }

    using Extents = cuda::std::dextents<cuda::std::ptrdiff_t, ndim>;
    using Mapping = cuda::std::layout_stride::mapping<Extents>;

    return cuda::std::mdspan<T, Extents, cuda::std::layout_stride>(arr.data(), Mapping(Extents(exts), strs));
}

} // namespace xsf
