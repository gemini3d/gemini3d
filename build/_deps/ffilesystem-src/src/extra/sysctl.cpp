#include "ffilesystem.h"

#if defined(FFS_DARWIN)
#include <cerrno>
#include <sys/sysctl.h>
#endif

bool
fs_is_rosetta()
{
#if defined(FFS_DARWIN)
// https://developer.apple.com/documentation/apple-silicon/about-the-rosetta-translation-environment
    int ret = 0;
    std::size_t size = sizeof(ret);

    if (sysctlbyname("sysctl.proc_translated", &ret, &size, nullptr, 0) < 0) {
        if (errno != ENOENT)
            fs_print_error("");
        return false;
    }

    return ret == 1;
#else
    return false;
#endif
}

// Microsoft Prism detection
// https://github.com/scivision/cmake-arm64-prism-rosetta/blob/main/is_prism.cpp
