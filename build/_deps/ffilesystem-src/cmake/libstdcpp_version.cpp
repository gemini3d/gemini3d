#if __has_include(<version>)
#include <version>
#elif __has_include(<ciso646>)
// < C++20 standard
#include <ciso646>
#endif

#include <cstdlib>
#include <iostream>

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


int main(){

#if defined(_LIBCPP_VERSION)
  std::cout << "LLVM ";
#elif defined(_GLIBCXX_RELEASE)
  std::cout << "GNU ";
#elif defined(_MSVC_STL_UPDATE)
  std::cout << "Microsoft ";
#endif

  std::cout << libcxx_release() << "\n";

  return EXIT_SUCCESS;
}
