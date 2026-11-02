#pragma once

#include "config.h"

#ifndef __CUDACC__

#define MDSPAN_IMPL_STANDARD_NAMESPACE xsf::cxx
// Force defining the parenthesis operator even when compiling with a compiler
// defaulting to C++ >= 23.
#define MDSPAN_USE_PAREN_OPERATOR 1
#include "third_party/kokkos/mdspan.hpp"

#else

#include <cuda/std/mdspan>

namespace xsf::cxx {

template <class IndexType, size_t... Extents>
using extents = cuda::std::extents<IndexType, Extents...>;

template <class IndexType, size_t Rank>
using dextents = cuda::std::dextents<IndexType, Rank>;

using layout_left = cuda::std::layout_left;
using layout_right = cuda::std::layout_right;
using layout_stride = cuda::std::layout_stride;

template <class ElementType>
using default_accessor = cuda::std::default_accessor<ElementType>;

template <
    class ElementType, class Extents, class LayoutPolicy = layout_right,
    class AccessorPolicy = default_accessor<ElementType>>
using mdspan = cuda::std::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>;

} // namespace xsf::cxx

#endif
