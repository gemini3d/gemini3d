#include "ffilesystem.h"
#include <string>

#include <boost/ut.hpp>


int main() {

using namespace boost::ut;

suite TestExpand = [] {
    "Expanduser"_test = [] {
        std::string const h = fs_get_homedir();
        expect(!h.empty() >> fatal) << "Home directory should not be empty";
        expect(fs_is_dir(h) >> fatal) << "Home directory should be a directory: " << h;

        expect(fs_expanduser("").empty());
        expect(fs_expanduser(".") == ".");
        expect(fs_expanduser("~") == h);
        expect(fs_expanduser("~/") == h);
        expect(fs_expanduser("~/test") == h + "/test");
        expect(fs_expanduser("~test") == "~test");
        expect(fs_expanduser("test~") == "test~");
        expect(fs_expanduser("test~test") == "test~test");

        expect(fs_expanduser("日本語") == "日本語");
    };
};
}
