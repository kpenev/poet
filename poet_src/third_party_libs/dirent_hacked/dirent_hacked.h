
#include "../../Core/SharedLibraryExportMacros.h"

#ifdef TOOLCHAIN_MSVC

#ifndef DIRENT_INCLUDED
#define DIRENT_INCLUDED


/*

    Declaration of POSIX directory browsing functions and types for Win32.

    Author:  Kevlin Henney (kevlin@acm.org, kevlin@curbralan.com)
    History: Created March 1997. Updated June 2003.
    Rights:  See end of file.
    
*/

extern "C"
{

typedef struct LIB_PUBLIC DIR DIR;

struct LIB_PUBLIC dirent
{
    char *d_name;
};

LIB_PUBLIC DIR           *opendir(const char *);
LIB_PUBLIC int           closedir(DIR *);
LIB_PUBLIC struct dirent *readdir(DIR *);
LIB_PUBLIC void          rewinddir(DIR *);

/*

    Copyright Kevlin Henney, 1997, 2003. All rights reserved.

    Permission to use, copy, modify, and distribute this software and its
    documentation for any purpose is hereby granted without fee, provided
    that this copyright and permissions notice appear in all copies and
    derivatives.
    
    This software is supplied "as is" without express or implied warranty.

    But that said, if there are any problems please get in touch.

*/

}

#endif

#else

#include <dirent.h>

#endif
