#include <string>
#include <string_view>

#include <locale>

#include "ffilesystem.h"


std::string fs_get_locale_name() {
// https://en.cppreference.com/w/cpp/locale/locale/name
    return std::locale("").name();
}
