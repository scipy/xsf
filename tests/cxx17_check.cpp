#include <xsf/config.h>

#ifdef _MSVC_LANG
static_assert(_MSVC_LANG == 201703L, "Tests must compile as C++17");
#else
static_assert(__cplusplus == 201703L, "Tests must compile as C++17");
#endif
