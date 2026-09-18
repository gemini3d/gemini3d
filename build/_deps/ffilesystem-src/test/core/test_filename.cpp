#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"HasFilename"_test = [] {

    expect(!fs_has_filename("") >> fatal);

    expect(!fs_has_filename("/"));
    expect(fs_has_filename("."));
    expect(!fs_has_filename("./"));
    expect(fs_has_filename(".."));
    expect(!fs_has_filename("../"));
    expect(fs_has_filename("a"));
    expect(!fs_has_filename("a/"));
    expect(fs_has_filename("a/."));
    expect(fs_has_filename("a/.."));
    expect(fs_has_filename("a/b"));
    expect(!fs_has_filename("a/b/"));
    expect(fs_has_filename("a/b/c"));
    expect(fs_has_filename("a/b sdc/some space"));
    expect(fs_has_filename("ab/.parent"));
    expect(fs_has_filename("ab/.parent.txt"));
    expect(fs_has_filename("a/b/../.parent.txt"));
    expect(fs_has_filename("/.fil"));
    expect(fs_has_filename("./日本語"));

if (fs_is_windows()){
    expect(!fs_has_filename("C:/"));
    expect(fs_has_filename(R"(C:\ab\asb)"));
    expect(fs_has_filename(R"(\\server\share\some space here.txt)"));

    if(fs_win32_long_paths_enabled()) {
        expect(!fs_has_filename(R"(\\.\)"));

        expect(!fs_has_filename(R"(\\?\C:\)"));
        expect(!fs_has_filename(R"(\\.\C:\)"));
        expect(fs_has_filename(R"(\\?\UNC\server\share)"));
        expect(fs_has_filename(R"(\\?\UNC\server\share\日本語)"));
        expect(fs_has_filename(R"(\\?\C:\some space here.txt)"));
  }
}

};

"FileName"_test = [] {

    expect(fs_file_name("") == "");
    expect(fs_file_name("/") == "");
    expect(fs_file_name(".") == ".");
    expect(fs_file_name("./") == "");
    expect(fs_file_name("..") == "..");
    expect(fs_file_name("../") == "");
    expect(fs_file_name("a") == "a");
    expect(fs_file_name("a/") == "");
    expect(fs_file_name("a/.") == ".");
    expect(fs_file_name("a/..") == "..");
    expect(fs_file_name("a/b") == "b");
    expect(fs_file_name("a/b/") == "");
    expect(fs_file_name("a/b/c") == "c");
    expect(fs_file_name("a/b sdc/some space") == "some space");
    expect(fs_file_name("ab/.parent") == ".parent");
    expect(fs_file_name("ab/.parent.txt") == ".parent.txt");
    expect(fs_file_name("a/b/../.parent.txt") == ".parent.txt");
    expect(fs_file_name("/.fil") == ".fil");
    expect(fs_file_name("./日本語") == "日本語");

if (fs_is_windows()){
    expect(fs_file_name("C:/") == "");
    expect(fs_file_name(R"(C:\ab\asb)") == "asb");

if(fs_win32_long_paths_enabled()){
    expect(fs_file_name(R"(\\?\)") == "");
    expect(fs_file_name(R"(\\.\)") == "");

    expect(fs_file_name(R"(\\?\C:\)") == "");
    expect(fs_file_name(R"(\\.\C:\)") == "");
    expect(fs_file_name(R"(\\?\UNC\server\share)") == "share");
    expect(fs_file_name(R"(\\?\UNC\server\share\日本語)") == "日本語");
    expect(fs_file_name(R"(\\server\share\some space here.txt)") == "some space here.txt");
    expect(fs_file_name(R"(\\?\C:\some space here.txt)") == "some space here.txt");
}
}

};

}
