#pragma once

#include "config.h"

namespace xsf {

template <typename T, typename = void>
constexpr inline bool valid_extent = false;

template <typename T>
constexpr inline bool
    valid_extent<T, cxx::void_t<typename T::index_type, typename T::size_type, typename T::rank_type>> = true;

struct relaxed_layout {
    template <typename TExtents>
    class mapping {
      public:
        static_assert(valid_extent<TExtents>, "TExtents must be a valid extents type.");

        using extents_type = TExtents;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;
        using strides_type = cxx::array<index_type, extents_type::rank()>;
        using layout_type = relaxed_layout;

      private:
        extents_type m_extents;
        strides_type m_strides;

        XSF_HOST_DEVICE constexpr auto wrap(index_type index, index_type extent) const {
            return (index >= index_type{0}) ? index : extent + index;
        }

        template <cxx::size_t Rank, cxx::size_t MaxRank>
        struct rank_counter {};

        template <cxx::size_t Rank, cxx::size_t MaxRank, typename TIndex, typename... TIndices>
        XSF_HOST_DEVICE constexpr auto
        compute_offset(rank_counter<Rank, MaxRank>, const TIndex &index, TIndices... indices) const noexcept {
            const auto shifted_index = wrap(index, m_extents.extent(Rank));
            return shifted_index * m_strides[Rank] + compute_offset(rank_counter<Rank + 1, MaxRank>{}, indices...);
        }

        template <typename TIndex>
        XSF_HOST_DEVICE constexpr auto compute_offset(
            rank_counter<extents_type::rank() - 1, extents_type::rank()>, const TIndex &index
        ) const noexcept {
            return wrap(index, m_extents.extent(extents_type::rank() - 1)) * m_strides[extents_type::rank() - 1];
        }

        XSF_HOST_DEVICE constexpr auto compute_offset(rank_counter<0, 0>) const noexcept { return index_type{0}; }

      public:
        XSF_HOST_DEVICE mapping() = default;
        XSF_HOST_DEVICE mapping(const extents_type &extents) : m_extents{extents}, m_strides{} {
            if constexpr (extents_type::rank() > 0) {
                m_strides[extents_type::rank() - 1] = index_type{1};
                for (auto i = extents_type::rank() - 1; i > 0; --i) {
                    m_strides[i - 1] = m_strides[i] * m_extents.extent(i);
                }
            }
        }
        XSF_HOST_DEVICE mapping(const extents_type &extents, const strides_type &strides)
            : m_extents{extents}, m_strides{strides} {}

        XSF_HOST_DEVICE mapping(const mapping &) noexcept = default;
        XSF_HOST_DEVICE mapping &operator=(const mapping &) = default;
        XSF_HOST_DEVICE mapping(mapping &&) noexcept = default;
        XSF_HOST_DEVICE mapping &operator=(mapping &&) noexcept = default;

        XSF_HOST_DEVICE const auto &extents() const noexcept { return m_extents; }
        XSF_HOST_DEVICE const auto &strides() const noexcept { return m_strides; }
        XSF_HOST_DEVICE constexpr auto stride(rank_type i) const noexcept { return m_strides[i]; }

        XSF_HOST_DEVICE constexpr auto required_span_size() const noexcept {
            auto size = index_type{1};
            for (auto i = rank_type{0}; i < extents_type::rank(); ++i) {
                if (m_extents.extent(i) == index_type{0}) {
                    return index_type{0};
                }
                size += (m_extents.extent(i) - 1) * ((m_strides[i] > 0) ? m_strides[i] : -m_strides[i]);
            }
            return size;
        }

        template <typename... TIndex>
        XSF_HOST_DEVICE constexpr auto operator()(TIndex... indices) const noexcept {
            static_assert(sizeof...(TIndex) == extents_type::rank(), "Number of indices must match rank.");
            return compute_offset(rank_counter<0, extents_type::rank()>{}, indices...);
        }

        XSF_HOST_DEVICE static constexpr auto is_unique() noexcept { return false; }
        XSF_HOST_DEVICE static constexpr auto is_exhaustive() noexcept { return false; }
        XSF_HOST_DEVICE static constexpr auto is_strided() noexcept { return true; }
        XSF_HOST_DEVICE static constexpr auto is_always_unique() noexcept { return false; }
        XSF_HOST_DEVICE static constexpr auto is_always_exhaustive() noexcept { return false; }
        XSF_HOST_DEVICE static constexpr auto is_always_strided() noexcept { return true; }

        friend XSF_HOST_DEVICE constexpr auto operator==(const mapping &lhs, const mapping &rhs) -> bool {
            return lhs.m_extents == rhs.m_extents && lhs.m_strides == rhs.m_strides;
        }
    };
};

} // namespace xsf
