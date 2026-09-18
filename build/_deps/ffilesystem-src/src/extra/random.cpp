#include <algorithm>            // for generate, generate_n
#include <array>                // for array
#include <functional>           // for ref
#include <iterator>             // for begin, end
#include <random>

#if __has_include(<ranges>)
#include <ranges>
#endif

#include <string>
#include <string_view>

#include "ffilesystem.h"


template <typename T = std::mt19937>
static auto fs_random_generator() -> T {
  auto constexpr seed_bytes = sizeof(typename T::result_type) * T::state_size;
  auto constexpr seed_len = seed_bytes / sizeof(std::seed_seq::result_type);
  auto seed = std::array<std::seed_seq::result_type, seed_len>();
  auto dev = std::random_device();
#if defined(__cpp_lib_ranges)
  std::ranges::generate(seed, std::ref(dev));
#else
  std::generate_n(std::begin(seed), seed_len, std::ref(dev));
#endif
  auto seed_seq = std::seed_seq(std::begin(seed), std::end(seed));
  return T{seed_seq};
}

std::string fs_generate_random_alphanumeric_string(const std::string::size_type len)
{
  constexpr std::string_view chars =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";

  // Keep one engine per thread to avoid repeated expensive initialization.
  thread_local std::mt19937 rng = fs_random_generator<std::mt19937>();
  auto dist = std::uniform_int_distribution<std::string::size_type>(0, chars.size() - 1);

  std::string result;
  result.resize(len);

#if defined(__cpp_lib_ranges)
  std::ranges::generate(result, [&]() { return chars[dist(rng)]; });
#else
  std::generate_n(std::begin(result), len, [&]() { return chars[dist(rng)]; });
#endif

  return result;
}
