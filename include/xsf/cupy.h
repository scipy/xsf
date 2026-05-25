#pragma once

#include "config.h"

/* Inject libcu++'s mdspan into std when compiling as CUDA, and use kokkos otherwise.
 * This allows this to work with both CUDA and HIP/ROCm. */
#ifdef __CUDACC__
#include <cuda/std/mdspan>
namespace std {

using cuda::std::layout_left;
using cuda::std::layout_right;
using cuda::std::layout_stride;

template <typename ElementType>
using default_accessor = cuda::std::default_accessor<ElementType>;

template <typename IndexType, size_t Rank>
using dextents = cuda::std::dextents<IndexType, Rank>;

template <typename IndexType, size_t... Extents>
using extents = cuda::std::extents<IndexType, Extents...>;

template <
    typename ElementType, typename Extents, typename LayoutPolicy = cuda::std::layout_right,
    typename AccessorPolicy = cuda::std::default_accessor<ElementType>>
using mdspan = cuda::std::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>;
} // namespace std
#else
#include "third_party/kokkos/mdspan.hpp"
#endif

namespace xsf {

template <typename T, int ndim, bool is_c_contiguous, bool index_32_bits, int core_ndim>
__device__ inline auto as_mdspan(const CArray<T, ndim, is_c_contiguous, index_32_bits, core_ndim> &arr) {
    std::array<std::ptrdiff_t, ndim> exts;
    std::array<std::ptrdiff_t, ndim> strs;

    for (int i = 0; i < ndim; ++i) {
        exts[i] = static_cast<std::ptrdiff_t>(arr.shape_[i]);
        strs[i] = static_cast<std::ptrdiff_t>(arr.strides_[i]) / sizeof(T);
    }

    using Extents = std::dextents<std::ptrdiff_t, ndim>;
    using Mapping = std::layout_stride::mapping<Extents>;

    return std::mdspan<T, Extents, std::layout_stride>(arr.data_, Mapping(Extents(exts), strs));
}

} // namespace xsf
