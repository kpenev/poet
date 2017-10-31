#ifndef __SHARED_LIBRARY_EXPORT_MACROS_H
#define __SHARED_LIBRARY_EXPORT_MACROS_H

    #if defined TOOLCHAIN_MSVC

        #warning "Compiling using MSVC"

        #ifdef BUILDING_DLL
            #define LIB_PUBLIC __declspec(dllexport)
        #else
            #define LIB_PUBLIC __declspec(dllimport)
        #endif

        #define LIB_LOCAL

    #elif defined TOOLCHAIN_GCC

        #warning "Compiling using GCC"

        #define LIB_PUBLIC __attribute__ ((visibility ("default")))
        #define LIB_LOCAL  __attribute__ ((visibility ("hidden")))

    #elif defined TOOLCHAIN_CLANG

        #warning "Compiling using clang"

        #define LIB_PUBLIC __attribute__ ((visibility ("default")))
        #define LIB_LOCAL  __attribute__ ((visibility ("hidden")))

    #else

        #warning "No toolchain defined"

    #endif

#endif
