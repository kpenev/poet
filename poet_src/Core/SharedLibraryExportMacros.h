#if defined TOOLCHAIN_MSVC

    #ifdef BUILDING_LIBRARY
        #define LIB_PUBLIC __declspec(dllexport)
    #else
        #define LIB_PUBLIC __declspec(dllimport)
    #endif

    #define LIB_LOCAL

#elif defined TOOLCHAIN_GCC

    #define LIB_PUBLIC __attribute__ ((visibility ("default")))
    #define LIB_LOCAL  __attribute__ ((visibility ("hidden")))

#elif defined TOOLCHAIN_CLANG

    #define LIB_PUBLIC __attribute__ ((visibility ("default")))
    #define LIB_LOCAL  __attribute__ ((visibility ("hidden")))

#else

    #warning "No toolchain defined"

#endif
