// Audit addition 2026-09-16. Apache-2.0; see ../LICENSE.
#ifndef GEMINI_CLI_UTILS_H
#define GEMINI_CLI_UTILS_H
#include <charconv>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include "gemini3d.h"
namespace gemini_cli {
inline bool positive_int(std::string_view s, int& value) {
  if (s.empty()) return false;
  const auto result = std::from_chars(s.data(), s.data()+s.size(), value);
  return result.ec == std::errc{} && result.ptr == s.data()+s.size() && value > 0;
}
inline bool copy_outdir(std::string_view s, char* target, std::size_t capacity) {
  if (s.size() >= capacity || s.find('\0') != std::string_view::npos) return false;
  std::memcpy(target, s.data(), s.size());
  target[s.size()] = '\0';
  return true;
}
inline bool parse_options(int argc, char** argv, params& p, int& d2, int& d3,
                          bool& help, std::string& error) {
  help = false;
  for (int i=2; i<argc; ++i) {
    const std::string_view a(argv[i]);
    if (a == "-d" || a == "-debug") p.debug = 1;
    else if (a == "-dryrun") p.dryrun = 1;
    else if (a == "-h" || a == "-help") help = true;
    else if (a == "-manual_grid") {
      if (i+2 >= argc || !positive_int(argv[i+1], d2) || !positive_int(argv[i+2], d3)) {
        error = "-manual_grid requires two positive integer counts";
        return false;
      }
      i += 2;
    } else {
      error = "Unsupported C++ option: " + std::string(a) + "; use gemini.bin for additional Fortran CLI options";
      return false;
    }
  }
  return true;
}
inline std::size_t checked_product(std::size_t a, std::size_t b) {
  if (b && a > std::numeric_limits<std::size_t>::max()/b)
    throw std::overflow_error("size_t allocation overflow");
  return a*b;
}
inline std::size_t cell_count(int x1, int x2, int x3) {
  if (x1 <= 0 || x2 <= 0 || x3 <= 0) throw std::invalid_argument("grid sizes must be positive");
  auto cells = std::size_t{1};
  for (const int x : {x1, x2, x3}) cells = checked_product(cells, static_cast<std::size_t>(x)+4);
  return cells;
}
inline std::size_t array_bytes(std::size_t cells, std::size_t fields) {
  return checked_product(checked_product(cells, fields), sizeof(double));
}
}
#endif
