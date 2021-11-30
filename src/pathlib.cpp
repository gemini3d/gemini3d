// MSVC 2017 or newer:
// cl /std:c++17 /EHsc

// Clang >= 9  or  g++ >= 8
// clang++ -std=c++17

#include <string>

std::string expanduser(std::string path){
// does not handle ~username/foo

if(path.front() != '~') return path;

std::string s, expanded = path;
#ifdef _WIN32
s = std::getenv("USERPROFILE");
#else
s = std::getenv("HOME");
#endif

if (!s.empty()) {
  expanded.erase(expanded.begin());
  expanded.insert(0, s);
}

return expanded;

}
