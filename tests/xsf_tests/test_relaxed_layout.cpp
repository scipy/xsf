#include "../testing_utils.h"

#include "xsf/mdspan.h"
#include "xsf/relaxed_layout.h"

using index_t = xsf::cxx::ptrdiff_t;

TEST_CASE("relaxed_layout computes zero offset for positive strides") {
    using Extents = xsf::cxx::dextents<index_t, 2>;
    using Mapping = xsf::relaxed_layout::mapping<Extents>;

    Extents extents{2, 3};
    Mapping::strides_type strides{3, 1};
    Mapping mapping{extents, strides};

    CHECK(mapping.offset() == 0);
    CHECK(mapping.required_span_size() == 6);

    CHECK(mapping(0, 0) == 0);
    CHECK(mapping(0, 1) == 1);
    CHECK(mapping(0, 2) == 2);
    CHECK(mapping(1, 0) == 3);
    CHECK(mapping(1, 1) == 4);
    CHECK(mapping(1, 2) == 5);
}

TEST_CASE("relaxed_layout computes canonical offset for a negative stride") {
    using Extents = xsf::cxx::dextents<index_t, 1>;
    using Mapping = xsf::relaxed_layout::mapping<Extents>;

    Extents extents{5};
    Mapping::strides_type strides{-1};
    Mapping mapping{extents, strides};

    CHECK(mapping.offset() == 4);
    CHECK(mapping.required_span_size() == 5);

    CHECK(mapping(0) == 4);
    CHECK(mapping(1) == 3);
    CHECK(mapping(2) == 2);
    CHECK(mapping(3) == 1);
    CHECK(mapping(4) == 0);
}

TEST_CASE("relaxed_layout handles mixed positive and negative strides") {
    using Extents = xsf::cxx::dextents<index_t, 2>;
    using Mapping = xsf::relaxed_layout::mapping<Extents>;

    Extents extents{2, 3};
    Mapping::strides_type strides{3, -1};
    Mapping mapping{extents, strides};

    CHECK(mapping.offset() == 2);
    CHECK(mapping.required_span_size() == 6);

    CHECK(mapping(0, 0) == 2);
    CHECK(mapping(0, 1) == 1);
    CHECK(mapping(0, 2) == 0);
    CHECK(mapping(1, 0) == 5);
    CHECK(mapping(1, 1) == 4);
    CHECK(mapping(1, 2) == 3);
}

TEST_CASE("relaxed_layout uses zero offset for an empty index space") {
    using Extents = xsf::cxx::dextents<index_t, 3>;
    using Mapping = xsf::relaxed_layout::mapping<Extents>;

    Extents extents{3, 0, 4};
    Mapping::strides_type strides{12, -4, 1};
    Mapping mapping{extents, strides};

    CHECK(mapping.offset() == 0);
    CHECK(mapping.required_span_size() == 0);
}

TEST_CASE("relaxed_layout works with static extents") {
    using Extents = xsf::cxx::extents<index_t, 2, 3>;
    using Mapping = xsf::relaxed_layout::mapping<Extents>;

    Mapping::strides_type strides{3, -1};
    Mapping mapping{Extents{}, strides};

    CHECK(mapping.offset() == 2);
    CHECK(mapping.required_span_size() == 6);

    CHECK(mapping(0, 0) == 2);
    CHECK(mapping(1, 2) == 3);
}

TEST_CASE("relaxed_layout can be used by mdspan with negative strides") {
    using Extents = xsf::cxx::dextents<index_t, 1>;
    using Mapping = xsf::relaxed_layout::mapping<Extents>;
    using Mdspan = xsf::cxx::mdspan<int, Extents, xsf::relaxed_layout>;

    int data[] = {0, 1, 2, 3, 4};

    Extents extents{5};
    Mapping::strides_type strides{-1};
    Mapping mapping{extents, strides};

    Mdspan view{data, mapping};

    CHECK(view(0) == 4);
    CHECK(view(1) == 3);
    CHECK(view(2) == 2);
    CHECK(view(3) == 1);
    CHECK(view(4) == 0);
}
