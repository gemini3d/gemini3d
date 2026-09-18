#include <chrono>
#include <string>
#include <algorithm>
#include <iostream>
#include <functional>
#include <cstdlib>
#include <vector>
#include <array>
#include <numeric>

// unconditional includes for benchmark output formatting
#include <iomanip> // for std::setw, std::setprecision
#include <sstream> // for std::ostringstream
#if __has_include(<format>)
#include <format>
#endif

#include <variant>
#include <unordered_map>

#include "ffilesystem.h"

constexpr int kBenchBatches = 9;


void print_cpp(std::chrono::duration<double> avg_t,
               std::chrono::duration<double> med_t,
               std::string_view path,
               std::string_view func,
               std::string_view w,
               bool b)
{
  const double avg_us =
    std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(avg_t).count();
  const double med_us =
    std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(med_t).count();

  std::string unit;
  double avg_v = 0.0;
  double med_v = 0.0;
  int prec = 3;

  if (avg_us < 1.0 && med_us < 1.0) {
    unit = "ns";
    avg_v = avg_us * 1000.0;
    med_v = med_us * 1000.0;
    prec = 1;
  } else {
    unit = "us";
    avg_v = avg_us;
    med_v = med_us;
    prec = 3;
  }

#ifdef __cpp_lib_format
  const auto avg_cell = std::format("{:.{}f} {}", avg_v, prec, unit);
  const auto med_cell = std::format("{:.{}f} {}", med_v, prec, unit);
  const auto result = w.empty() ? std::string(b ? "1" : "0") : std::string(w);

  std::cout << std::format("{:>12}{:>12}  {}({}) = {}\n",
                           avg_cell,
                           med_cell,
                           func,
                           path,
                           result);
#else
  const auto old_flags = std::cout.flags();
  const auto old_prec = std::cout.precision();

  std::ostringstream avg_cell;
  avg_cell << std::fixed << std::setprecision(prec) << avg_v << ' ' << unit;

  std::ostringstream med_cell;
  med_cell << std::fixed << std::setprecision(prec) << med_v << ' ' << unit;

  std::cout << std::right << std::setw(12) << avg_cell.str()
            << std::setw(12) << med_cell.str()
            << "  " << func << "(" << path << ") = ";
  if (w.empty())
    std::cout << b;
  else
    std::cout << w;
  std::cout << "\n";

  std::cout.flags(old_flags);
  std::cout.precision(old_prec);
#endif
}


std::chrono::duration<double> bench_cpp(int n, std::string_view path, std::string_view fname, bool verbose)
{
auto invalid = std::chrono::duration<double>::max();

constexpr bool strict = false;

using fs_function = std::function<std::variant<std::string, bool>(std::string_view)>;

std::unordered_map<std::string_view, fs_function> fs_function_map = {
  {"absolute", [=](std::string_view p) { return fs_absolute(p); }},
  {"canonical", [=](std::string_view p) { return fs_canonical(p, strict); }},
  {"resolve", [=](std::string_view p) { return fs_resolve(p, strict); }},
  {"drop_slash", [](std::string_view p) { return fs_drop_slash(p); }},
  {"parent", [](std::string_view p) { return fs_parent(p); }},
  {"file_name", [](std::string_view p) { return fs_file_name(p); }},
  {"suffix", [](std::string_view p) { return fs_suffix(p); }},
  {"normal", [](std::string_view p) { return fs_normal(p); }},
  {"reserved", [](std::string_view p) { return fs_is_reserved(p); }},
  {"exists", [](std::string_view p) { return fs_exists(p); }},
  {"is_dir", [](std::string_view p) { return fs_is_dir(p); }},
  {"is_char", [](std::string_view p) { return fs_is_char_device(p); }},
  {"is_file", [](std::string_view p) { return fs_is_file(p); }},
  {"is_symlink", [](std::string_view p) { return fs_is_symlink(p); }},
  {"read_symlink", [](std::string_view p) { return fs_read_symlink(p); }},
  {"which", [](std::string_view p) { return fs_which(p); }},
  {"homedir", [](std::string_view) { return fs_get_homedir(); }},
  {"expanduser", [](std::string_view p) { return fs_expanduser(p); }},
  {"cwd", [](std::string_view) { return fs_get_cwd(); }}
};

std::string h;
bool b = false;

auto it = fs_function_map.find(fname);
if (it == fs_function_map.end()) {
  std::cerr << "Error: unknown function " << fname << "\n";
  return invalid;
}

auto result = it->second(path);
if (std::holds_alternative<std::string>(result)) {
  h = std::get<std::string>(result);
  if (h.empty()) {
    std::cerr << "ERROR: " << fname << " " << path << " failed on warmup\n";
    return invalid;
  }
} else {
  b = std::get<bool>(result);
}

std::vector<std::chrono::duration<double>> per_call;
per_call.reserve(kBenchBatches);

for (int batch = 0; batch < kBenchBatches; ++batch) {
  auto t0 = std::chrono::steady_clock::now();

  for (int i = 0; i < n; ++i) {
    auto iter_result = it->second(path);
    if (std::holds_alternative<std::string>(iter_result))
      h = std::get<std::string>(iter_result);
    else
      b = std::get<bool>(iter_result);
  }

  auto t1 = std::chrono::steady_clock::now();
  per_call.push_back((t1 - t0) / static_cast<double>(n));
}

std::chrono::duration<double> sum{0};
for (const auto& d : per_call)
  sum += d;
const auto avg = sum / static_cast<double>(per_call.size());

std::sort(per_call.begin(), per_call.end());
std::chrono::duration<double> med;
if (per_call.size() % 2 == 1) {
  med = per_call[per_call.size() / 2];
} else {
  const auto mid = per_call.size() / 2;
  med = (per_call[mid - 1] + per_call[mid]) / 2.0;
}

if (verbose)
  print_cpp(avg, med, path, fname, h, b);

return avg;
}


int main(int argc, char** argv){

if (!fs_is_optimized())
  std::cerr << "WARNING: ffilesystem might not have been compiled with optimizations\n";

int n = 1000;
if(argc > 1)
    n = std::stoi(argv[1]);

std::string_view path;

std::cout << fs_compiler() << "\n";
std::cout << "Ffilesystem backend: " << fs_backend() << "  C++ standard: " << fs_cpp_lang() << "  C standard: " << fs_c_lang() << "\n";
std::cout << "Benchmark C++ standard " << __cplusplus << "\n";
std::cout << "Optimized: " << fs_is_optimized() << "\n";
std::cout << "C++20 <format> support: " << fs_cpp_format() << "\n";
std::cout << "Benchmark parameters: n=" << n << ", batches=" << kBenchBatches << "\n";

#ifdef __cpp_lib_format
std::cout << std::format("{:>12}{:>12}\n", "avg", "median");
#else
std::cout << "         avg          median\n";
#endif

std::vector<std::string_view> funcs;
if(argc > 3)
  funcs = {argv[3]};
else
  funcs = {"absolute", "canonical", "resolve", "which", "expanduser", "normal", "cwd",
           "homedir", "parent", "file_name", "reserved", "drop_slash",
           "exists", "is_dir", "is_file", "is_symlink"};

// only use reliable calls e.g. not to non-existent paths to avoid avalanche of error messages
// not doing read_symlink because it generates errors without a test symlink

for (std::string_view func : funcs)
  {
  // in sorted ascending order for binary search
  constexpr std::array<std::string_view, 8> tf = {
    "canonical", "drop_slash", "expanduser", "file_name", "normal", "parent", "resolve"
  };

  if (argc > 2)
    path = argv[2];
  else {
    if (std::binary_search(tf.begin(), tf.end(), func))
      path = "./..";
    else if (func == "which")
      path = (fs_is_windows()) ? "cmake.exe" : "sh";
    else
      path = ".";
  }

 bench_cpp(n, path, func, true);
}


return EXIT_SUCCESS;

}
