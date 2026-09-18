#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"memory"_test = [] {
  expect(gt(fs_get_free_memory(), 0));
};

"total_memory"_test = [] {
  expect(gt(fs_total_sys_memory(), 0));
};

"free_less_than_total"_test = [] {
  expect(le(fs_get_free_memory(), fs_total_sys_memory()));
};

"max_open_files"_test = [] {
  expect(gt(fs_get_max_open_files(), 0));
};
}
