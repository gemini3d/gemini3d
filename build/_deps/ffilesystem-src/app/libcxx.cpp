// https://gcc.gnu.org/onlinedocs/libstdc++/manual/abi.html
// https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_macros.html
// https://en.cppreference.com/w/cpp/header/ciso646

#include <cstddef>  // IWYU pragma: keep

#if __has_include(<version>)
#include <version>  // IWYU pragma: keep
#elif __has_include(<ciso646>)
// < C++20 standard
#include <ciso646>
#endif

#include <cstdlib>
#include <iostream>

// #include <iosfwd> more standard to use <version> and <ciso646>
// https://libcxx.llvm.org/DesignDocs/ABIVersioning.html
// _LIBCPP_ABI_VER

// Microsoft STL doesn't have a version number.

unsigned long libcxx_release()
{
#if defined(_LIBCPP_VERSION)  // LLVM libc++
  return _LIBCPP_VERSION;
#elif defined(_GLIBCXX_RELEASE)  // GNU libstdc++
  return _GLIBCXX_RELEASE;
#elif defined(__MSVC_STL_UPDATE)  // Microsoft STL
  return _MSVC_STL_UPDATE;
#else
  return 0;
#endif
}


int main()
{
  std::cout
#ifdef _LIBCPP_VERSION
  << "LLVM libc++"
#elif defined(_GLIBCXX_RELEASE)
  << "GNU GCC libstdc++"
#elif defined(_MSVC_STL_UPDATE)
  << "Microsoft STL"
#endif
  << " " << libcxx_release() << "\n";

  return EXIT_SUCCESS;
}
