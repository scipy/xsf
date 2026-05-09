#include "../testing_utils.h"
#include <xsf/stats.h>

/*
# Reference values computed with scipy.stats._kendall_dis
import numpy as np
from scipy.stats._stats import _kendall_dis

rng = np.random.default_rng(123456789)
n = 150
x_no_ties = rng.permutation(n) + 1
y_no_ties = rng.permutation(n) + 1
x_with_ties = rng.integers(1, n + 1, size=n)
y_with_ties = rng.integers(1, n + 1, size=n)

# check "no_ties" arrays have no ties and "with_ties" arrays have ties
assert not len(np.unique(x_no_ties)) < len(x_no_ties)
assert len(np.unique(x_with_ties)) < len(x_with_ties)
assert not len(np.unique(y_no_ties)) < len(y_no_ties)
assert len(np.unique(y_with_ties)) < len(y_with_ties)

for key_x, x in {"no_ties": x_no_ties, "with_ties": x_with_ties}.items():
    for key_y, y in {"no_ties": y_no_ties, "with_ties": y_with_ties}.items():
        print(f"x {key_x} and y {key_y}: disc = {_kendall_dis(x, y)}")
*/
TEST_CASE("kendall_dis test", "[kendall_dis][xsf_tests]") {
    std::vector<intptr_t> x_no_ties = {
        81,  147, 93,  88,  65,  140, 89,  139, 72,  17,  43,  79,  19,  35,  101, 6,   53,  18,  149, 120, 90,  26,
        115, 42,  103, 98,  97,  96,  84,  54,  32,  102, 27,  108, 8,   68,  123, 94,  64,  134, 22,  124, 44,  41,
        50,  117, 116, 69,  138, 58,  91,  85,  21,  37,  49,  33,  80,  130, 75,  132, 83,  15,  51,  82,  127, 24,
        95,  114, 56,  109, 59,  34,  150, 119, 133, 107, 14,  4,   77,  31,  100, 92,  137, 122, 128, 136, 145, 55,
        48,  63,  29,  36,  105, 76,  5,   125, 38,  141, 113, 39,  45,  11,  148, 71,  111, 66,  28,  25,  110, 52,
        112, 118, 62,  74,  16,  3,   144, 67,  87,  57,  7,   70,  104, 30,  131, 86,  73,  9,   23,  142, 106, 60,
        143, 78,  129, 126, 99,  46,  10,  12,  2,   1,   121, 20,  40,  47,  146, 135, 13,  61
    };
    std::vector<intptr_t> y_no_ties = {
        39,  8,   58,  40,  55,  91,  57,  79,  117, 149, 29,  147, 150, 60,  93, 138, 36, 112, 122, 145, 137, 34,
        87,  19,  126, 26,  102, 86,  41,  105, 111, 146, 76,  143, 140, 136, 32, 18,  75, 67,  10,  17,  25,  22,
        119, 134, 43,  77,  14,  84,  101, 80,  103, 144, 46,  6,   124, 73,  38, 1,   68, 49,  74,  133, 114, 12,
        116, 118, 92,  63,  54,  106, 5,   121, 82,  33,  109, 69,  96,  16,  97, 95,  89, 61,  35,  27,  20,  37,
        56,  127, 108, 141, 142, 104, 135, 30,  13,  65,  45,  130, 83,  51,  23, 78,  44, 59,  110, 52,  81,  129,
        128, 2,   15,  53,  4,   132, 21,  64,  123, 148, 48,  139, 115, 100, 72, 71,  24, 90,  120, 9,   88,  98,
        113, 7,   70,  131, 66,  42,  62,  125, 31,  28,  3,   50,  107, 85,  99, 47,  11, 94
    };
    std::vector<intptr_t> x_with_ties = {
        45,  34,  80,  111, 115, 17,  129, 24,  29,  16,  98,  133, 59,  71,  131, 99,  41,  123, 44,  38,  42,  106,
        144, 82,  89,  59,  32,  75,  33,  99,  13,  84,  96,  88,  28,  85,  141, 101, 123, 65,  55,  70,  78,  124,
        4,   8,   83,  140, 24,  86,  40,  19,  120, 44,  111, 77,  105, 97,  94,  34,  46,  74,  119, 38,  41,  148,
        27,  145, 88,  88,  56,  30,  118, 35,  52,  123, 34,  97,  6,   15,  140, 27,  99,  54,  104, 18,  103, 102,
        12,  140, 134, 129, 90,  143, 95,  43,  70,  102, 22,  64,  62,  135, 147, 46,  123, 147, 55,  118, 70,  23,
        43,  41,  65,  134, 10,  29,  124, 105, 110, 119, 104, 102, 34,  26,  3,   15,  16,  52,  45,  24,  79,  49,
        89,  36,  149, 15,  113, 58,  108, 147, 28,  76,  99,  47,  95,  55,  146, 39,  42,  50
    };
    std::vector<intptr_t> y_with_ties = {
        97,  89, 31,  48,  54,  30,  120, 48,  61,  116, 35,  77,  33,  4,   139, 46,  52,  61,  104, 25,  109, 128,
        44,  77, 1,   28,  87,  43,  36,  139, 41,  76,  40,  82,  93,  49,  59,  110, 70,  110, 11,  130, 88,  65,
        27,  1,  60,  70,  98,  126, 103, 53,  109, 109, 133, 136, 25,  117, 36,  45,  20,  7,   89,  52,  129, 39,
        6,   93, 132, 120, 42,  66,  148, 141, 94,  106, 127, 86,  116, 41,  31,  47,  63,  23,  111, 48,  63,  108,
        120, 48, 46,  72,  18,  147, 65,  5,   59,  52,  111, 64,  20,  23,  105, 99,  92,  42,  107, 82,  16,  87,
        65,  37, 87,  31,  111, 15,  5,   76,  125, 91,  127, 30,  67,  100, 9,   136, 106, 68,  14,  74,  143, 97,
        86,  18, 28,  51,  23,  145, 133, 28,  94,  84,  72,  46,  103, 68,  12,  40,  80,  138
    };

    // x_no_ties and y_no_ties
    const int64_t res_no_no = xsf::kendall_dis(x_no_ties, y_no_ties);
    const int64_t expected_no_no = 5943;
    CAPTURE(x_no_ties, y_no_ties, res_no_no, expected_no_no);
    REQUIRE(res_no_no == expected_no_no);

    // x_no_ties and y_with_ties
    const int64_t res_no_with = xsf::kendall_dis(x_no_ties, y_with_ties);
    const int64_t expected_no_with = 5444;
    CAPTURE(x_no_ties, y_with_ties, res_no_with, expected_no_with);
    REQUIRE(res_no_with == expected_no_with);

    // x_with_ties and y_no_ties
    const int64_t res_with_no = xsf::kendall_dis(x_with_ties, y_no_ties);
    const int64_t expected_with_no = 5942;
    CAPTURE(x_with_ties, y_no_ties, res_with_no, expected_with_no);
    REQUIRE(res_with_no == expected_with_no);

    // x_with_ties and y_with_ties
    const int64_t res_with_with = xsf::kendall_dis(x_with_ties, y_with_ties);
    const int64_t expected_with_with = 5443;
    CAPTURE(x_with_ties, y_with_ties, res_with_with, expected_with_with);
    REQUIRE(res_with_with == expected_with_with);
}
