#include <iostream>
#include <ostream> // for std::endl
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>
#include <functional>
#include <ctime>
#include <system_error>
#include <variant>
#include <unordered_map>
#include <type_traits>
#include <sstream>
#include <optional>

#include <chrono> // IWYU pragma: keep
// needed to std::format() std::filesystem::file_time_type

#if __has_include(<format>)
#include <format> // IWYU pragma: keep
#endif

#if defined(HAVE_CXX_FILESYSTEM)
#include <filesystem>
namespace Filesystem = std::filesystem;
#endif

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <crtdbg.h>
#endif

#include "ffilesystem.h"


static bool no_arg(std::string_view fun){

  static const std::unordered_map<std::string_view, std::function<bool()>> mbool =
  {
    {"optimized", fs_is_optimized},
    {"is_admin", fs_is_admin},
    {"is_android", fs_is_android},
    {"is_bsd", fs_is_bsd},
    {"is_linux", fs_is_linux},
    {"is_macos", fs_is_macos},
    {"is_rosetta", fs_is_rosetta},
    {"is_unix", fs_is_unix},
    {"is_windows", fs_is_windows},
    {"is_mingw", fs_is_mingw},
    {"is_cygwin", fs_is_cygwin},
    {"has_statx", fs_has_statx},
    {"long_paths", fs_win32_long_paths_enabled}
  };

using fs_function = std::function<std::variant<std::string, bool, int, char, long, unsigned long long>()>;

static const std::unordered_map<std::string_view, fs_function> fs_function_map = {
  {"backend", []() { return fs_backend(); }},
  {"is_wsl", []() { return fs_is_wsl(); }},
  {"pathsep", []() { return fs_pathsep(); }},
  {"cpp_lang", []() { return fs_cpp_lang(); }},
  {"c_lang", []() { return fs_c_lang(); }},
  {"locale", []() { return fs_get_locale_name(); }},
  {"username", []() { return fs_get_username(); }},
  {"hostname", []() { return fs_hostname(); }},
  {"shell", []() { return fs_get_shell(); }},
  {"stdin_tty", []() { return fs_stdin_tty(); }},
  {"arch", []() { return fs_cpu_arch(); }},
  {"profile", []() { return fs_get_profile_dir(); }},
  {"os_version", []() { return fs_os_version(); }},
  {"libpath", []() { return fs_lib_path(); }},
  {"exepath", []() { return fs_exe_path(); }},
  {"compiler", []() { return fs_compiler(); }},
  {"homedir", []() { return fs_get_homedir(); }},
  {"terminal", []() { return fs_get_terminal(); }},
  {"tempdir", []() { return fs_get_tempdir(); }},
  {"max_path", []() { return static_cast<unsigned long long>(fs_get_max_path()); }},
  {"cwd", []() { return fs_get_cwd(); }},
  {"free_ram", []() { return fs_get_free_memory(); }},
  {"total_ram", []() { return fs_total_sys_memory(); }},
  {"max_open_files", []() { return fs_get_max_open_files(); }}
};

  bool ok = true;

  auto it = fs_function_map.find(fun);
  if (it != fs_function_map.end()) {
    auto result = it->second();
    if (std::holds_alternative<std::string>(result))
      std::cout << std::get<std::string>(result);
    else if (std::holds_alternative<int>(result))
      std::cout << std::get<int>(result);
    else if (std::holds_alternative<bool>(result))
      std::cout << std::get<bool>(result);
    else if (std::holds_alternative<char>(result))
      std::cout << std::get<char>(result);
    else if (std::holds_alternative<long>(result))
      std::cout << std::get<long>(result);
    else if (std::holds_alternative<unsigned long long>(result))
      std::cout << std::get<unsigned long long>(result);
    else
      ok = false;

  } else if (const auto bit = mbool.find(fun); bit != mbool.end())
    std::cout << bit->second();
  else if (fun == "pid")
    std::cout << fs_getpid();
  else
    ok = false;

  if(ok)
    std::cout << "\n";

  return ok;
}


static bool one_arg(std::string_view fun, std::string_view a1)
{
  // each possible return type for the function
  using fs_one_arg_function = std::function<std::variant<std::string, bool, int, std::uintmax_t>(std::string_view)>;

  static const std::unordered_map<std::string_view, fs_one_arg_function> fs_one_arg_function_map = {
    {"lexically_normal", [](std::string_view a1) { return fs_lexically_normal(a1); }},
    {"make_preferred", [](std::string_view a1) { return fs_make_preferred(a1); }},
    {"parent", [](std::string_view a1) { return fs_parent(a1); }},
    {"suffix", [](std::string_view a1) { return fs_suffix(a1); }},
    {"canonical", [](std::string_view a1) { return fs_canonical(a1, true); }},
    {"weakly_canonical", [](std::string_view a1) { return fs_canonical(a1, false); }},
    {"resolve", [](std::string_view a1) { return fs_resolve(a1, true); }},
    {"weakly_resolve", [](std::string_view a1) { return fs_resolve(a1, false); }},
    {"normal", [](std::string_view a1) { return fs_normal(a1); }},
    {"type", [](std::string_view a1) { return fs_filesystem_type(a1); }},
    {"is_reserved", [](std::string_view a1) { return fs_is_reserved(a1); }},
    {"is_safe", [](std::string_view a1) { return fs_is_safe_name(a1); }},
    {"is_appexec", [](std::string_view a1) { return fs_is_appexec_alias(a1); }},
    {"long", [](std::string_view a1) { return fs_longname(a1); }},
    {"short", [](std::string_view a1) { return fs_shortname(a1); }},
    {"touch", [](std::string_view a1) { return "touch " + std::string(a1) + " " + std::to_string(fs_touch(a1)); }},
    {"set_modtime", [](std::string_view a1) { return fs_set_modtime(a1); }},
    {"getenv", [](std::string_view a1) { auto e = fs_getenv(a1);
      if (e){ return e.value(); };
      std::cerr << a1 << " environment variable not set\n";
      return std::string();
    }},
    {"setenv", [](std::string_view a1) { std::cout << "unset " << a1 << "\n"; return fs_setenv(a1, ""); }},
    {"realpath", [](std::string_view a1) { return fs_realpath(a1); }},
    {"as_posix", [](std::string_view a1) { std::string s(a1); fs_as_posix(s); return s; }},
    {"as_windows", [](std::string_view a1) { std::string s(a1); fs_as_windows(s); return s; }},
    {"max_component", [](std::string_view a1) { return static_cast<std::uintmax_t>(fs_max_component(a1)); }},
    {"mkdir", [](std::string_view a1) { return fs_mkdir(a1); }},
    {"owner", [](std::string_view a1) { return fs_get_owner_name(a1) + "\n" + fs_get_owner_group(a1); }},
    {"expanduser", [](std::string_view a1) { return fs_expanduser(a1); }},
    {"final_path", [](std::string_view a1) { return fs_win32_final_path(a1); }},
    {"root", [](std::string_view a1) { return fs_root(a1); }},
    {"drop_slash", [](std::string_view a1) { return fs_drop_slash(a1); }},
    {"root_name", [](std::string_view a1) { return fs_root_name(a1); }},
    {"filename", [](std::string_view a1) { return fs_file_name(a1); }},
    {"has_filename", [](std::string_view a1) { return fs_has_filename(a1); }},
    {"is_absolute", [](std::string_view a1) { return fs_is_absolute(a1); }},
    {"is_exe", [](std::string_view a1) { return fs_is_exe(a1); }},
    {"is_exe_bin", [](std::string_view a1) { return fs_is_executable_binary(a1); }},
    {"is_case_sensitive", [](std::string_view a1) { return fs_is_case_sensitive(a1); }},
    {"is_dir", [](std::string_view a1) { return fs_is_dir(a1); }},
    {"is_char", [](std::string_view a1) { return fs_is_char_device(a1); }},
    {"is_file", [](std::string_view a1) { return fs_is_file(a1); }},
    {"is_fifo", [](std::string_view a1) { return fs_is_fifo(a1); }},
    {"is_symlink", [](std::string_view a1) { return fs_is_symlink(a1); }},
    {"is_removable", [](std::string_view a1) { return fs_is_removable(a1); }},
    {"is_readable", [](std::string_view a1) { return fs_is_readable(a1); }},
    {"is_writable", [](std::string_view a1) { return fs_is_writable(a1); }},
    {"device", [](std::string_view a1) { return static_cast<std::uintmax_t>(fs_st_dev(a1)); }},
    {"mode", [](std::string_view a1) { return static_cast<std::uintmax_t>(fs_st_mode(a1)); }},
    {"inode", [](std::string_view a1) { return static_cast<std::uintmax_t>(fs_inode(a1)); }},
    {"perm", [](std::string_view a1) { return fs_get_permissions(a1); }},
    {"read_symlink", [](std::string_view a1) { return fs_read_symlink(a1); }},
    {"stem", [](std::string_view a1) { return fs_stem(a1); }},
    {"exists", [](std::string_view a1) { return fs_exists(a1); }},
    {"blksize", [](std::string_view a1) { return static_cast<std::uintmax_t>(fs_get_blksize(a1)); }},
    {"absolute", [](std::string_view a1) { return fs_absolute(a1); }},
    {"is_empty", [](std::string_view a1) { return fs_is_empty(a1); }},
    {"remove", [](std::string_view a1) { return fs_remove(a1); }},
    {"which", [](std::string_view a1) { return fs_which(a1, "", false); }},
    {"which_all", [](std::string_view a1) { return fs_which(a1, "", true); }},
  };

  static const std::unordered_map<std::string_view, std::function<std::uintmax_t(std::string_view)>> m1uintm =
  {
    {"file_size", [](std::string_view a1) { return fs_file_size(a1); }},
    {"space_available", [](std::string_view a1) { return fs_space_available(a1); }},
    {"space_capacity", [](std::string_view a1) { return fs_space_capacity(a1); }},
    {"hard_link_count", [](std::string_view a1) { return fs_hard_link_count(a1); }}
  };

  bool ok = true;

  auto it = fs_one_arg_function_map.find(fun);
  // if-else each possible return type for the function
  if (it != fs_one_arg_function_map.end()) {
    auto r = it->second(a1);
    if (std::holds_alternative<std::string>(r))
      std::cout << std::get<std::string>(r);
    else if (std::holds_alternative<bool>(r))
      std::cout << std::get<bool>(r);
    else if (std::holds_alternative<int>(r))
      std::cout << std::get<int>(r);
    else if (std::holds_alternative<std::uintmax_t>(r))
      std::cout << std::get<std::uintmax_t>(r);
    else
      ok = false;

  } else if (const auto mit = m1uintm.find(fun); mit != m1uintm.end()) {
    if(auto r = mit->second(a1); r != fs_unknown_size)
      std::cout << r;
  } else if (fun == "modtime"){

#if defined(HAVE_CXX_FILESYSTEM) && defined(__cpp_lib_format) // C++20
    const auto t = fs_get_modtime_fs(a1);
    if(t)
      std::cout << std::format("{}\n", t.value());
#else
    const std::time_t t = fs_get_modtime(a1);
    #if defined(_MSC_VER)
      std::string buf;
      buf.resize(26);
      ctime_s(buf.data(), buf.size(), &t);
      std::cout << buf;
    #else
      std::cout << std::ctime(&t); // NOSONAR
    #endif
#endif
  } else if (fun == "fs_modtime")
    std::cout << fs_get_modtime(a1);
  else if (fun == "random")
    std::cout << fs_generate_random_alphanumeric_string(std::strtoul(a1.data(), nullptr, 10));
  else if (fun == "split"){
    std::vector<std::string> v = fs_split(a1);
    for (const auto &s : v)
      std::cout << s << "\n";
  } else if (fun == "normal_vector"){
    std::vector<std::string> v = fs_normal_vector(a1);
    for (const auto &s : v)
      std::cout << s << "\n";
  } else if (fun == "chdir" || fun == "set_cwd") {
    auto cwd = fs_get_cwd();
    if(!cwd.empty()){
      std::cout << "cwd: " << cwd;
      if(fs_set_cwd(a1)){
        cwd = fs_get_cwd();
        if(!cwd.empty())
          std::cout << "new cwd: " << cwd;
        else
          std::cerr << "ERROR get_cwd() after chdir\n";
      }
    } else {
      std::cerr << "ERROR get_cwd() before chdir\n";
    }
  } else if (fun == "ls") {
#if defined(HAVE_CXX_FILESYSTEM)
    std::error_code ec;
    for (auto const& dir_entry :
         Filesystem::directory_iterator{a1,
          Filesystem::directory_options::skip_permission_denied, ec}){

      if(ec)
        std::cerr << "ERROR: " << ec.message() << " " << ec.value() << "\n";

      Filesystem::path p = dir_entry.path();
      std::cout << p;

      if (Filesystem::is_regular_file(p, ec)){
        if (const auto &s = Filesystem::file_size(p, ec); s && !ec)
          std::cout << " " << s;
      }

      std::cout << " " << fs_get_permissions(p.generic_string());
      std::cout << std::endl;
    }
#else
      std::cerr << "ERROR: ls requires C++17 filesystem\n";
#endif
  } else
    ok = false;

  if (ok)
    std::cout << "\n";

  return ok;
}


static bool two_arg(std::string_view fun, std::string_view a1, std::string_view a2)
{
  using fs_two_arg_function = std::function<std::variant<std::string, bool>(std::string_view, std::string_view)>;

  static const std::unordered_map<std::string_view, fs_two_arg_function> fs_two_arg_function_map = {
    {"which", [](std::string_view a1, std::string_view a2) { return fs_which(a1, a2, false); }},
    {"which_all", [](std::string_view a1, std::string_view a2) { return fs_which(a1, a2, true); }},
    {"same", [](std::string_view a1, std::string_view a2) { return fs_equivalent(a1, a2); }},
    {"join", [](std::string_view a1, std::string_view a2) { return fs_join(a1, a2); }},
    {"with_suffix", [](std::string_view a1, std::string_view a2) { return fs_with_suffix(a1, a2); }},
    {"is_prefix", [](std::string_view a1, std::string_view a2) { return fs_is_prefix(a1, a2); }},
    {"is_subdir", [](std::string_view a1, std::string_view a2) { return fs_is_subdir(a1, a2); }},
    {"relative", [](std::string_view a1, std::string_view a2) { return fs_relative_to(a1, a2); }},
    {"proximate", [](std::string_view a1, std::string_view a2) { return fs_proximate_to(a1, a2); }},
    {"setenv", [](std::string_view a1, std::string_view a2) { return fs_setenv(a1, a2); }},
    {"copy", [](std::string_view a1, std::string_view a2) { return fs_copy_file(a1, a2, false); }},
    {"create_symlink", [](std::string_view a1, std::string_view a2) { return fs_create_symlink(a1, a2); }},
    {"absolute", [](std::string_view a1, std::string_view a2) { return fs_absolute(a1, a2); }},
    {"rename", [](std::string_view a1, std::string_view a2) { return fs_rename(a1, a2); }}
  };

  bool ok = true;

  auto it = fs_two_arg_function_map.find(fun);
  if (it != fs_two_arg_function_map.end()) {
    auto r = it->second(a1, a2);
    if (std::holds_alternative<std::string>(r))
      std::cout << std::get<std::string>(r);
    else if (std::holds_alternative<bool>(r))
      std::cout << std::get<bool>(r);
    else
      ok = false;

  } else if (fun == "set_perm"){
    int r = 0;
    int w = 0;
    int x = 0;

#ifdef __cpp_lib_string_contains
    if(a2.contains("+r"))
      r = 1;
    if(a2.contains("+w"))
      w = 1;
    if(a2.contains("+x"))
      x = 1;
    if(a2.contains("-r"))
      r = -1;
    if(a2.contains("-w"))
      w = -1;
    if(a2.contains("-x"))
      x = -1;
#else
    if(a2.find("+r") != std::string_view::npos)
      r = 1;
    if(a2.find("+w") != std::string_view::npos)
      w = 1;
    if(a2.find("+x") != std::string_view::npos)
      x = 1;
    if(a2.find("-r") != std::string_view::npos)
      r = -1;
    if(a2.find("-w") != std::string_view::npos)
      w = -1;
    if(a2.find("-x") != std::string_view::npos)
      x = -1;
#endif

    auto p = fs_get_permissions(a1);
    if(p.empty()){
      std::cerr << "ERROR get_permissions(" << a1 << ") before chmod\n";
      return false;
    }

    std::cout << "before chmod " << a1 << " " << p << "\n";
    std::cout << "setting r,w,x " << r << " " << w << " " << x << "\n";
    if(!fs_set_permissions(a1, r, w, x))
      std::cerr << "ERROR set_permissions(" << a1 << ")\n";

    p = fs_get_permissions(a1);
    if(p.empty())
      std::cerr << "ERROR get_permissions(" << a1 << ") after chmod\n";
    else
      std::cout << "after chmod " << a1 << " " << p << "\n";
  } else
    ok = false;

  if(ok)
    std::cout << "\n";

  return ok;
}


int main()
{
#ifdef _MSC_VER
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
  _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
  _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
#endif

const std::string cwd = fs_get_cwd();

std::cout << "Ffilesystem. Backend: " << fs_backend() << "\n";
std::cout << fs_compiler() << "\n";
std::cout << "Ffilesystem C++ standard " << fs_cpp_lang() << "  CLI C++ standard " << __cplusplus << "\n";
std::cout << "Ffilesystem C++ <format> support: " << fs_cpp_format() << "\n";
std::cout << "Ffilesystem C++ <ranges> support: " << fs_cpp_ranges() << "\n";
std::cout << "Ffilesystem C standard " << fs_c_lang() << "\n";
std::cout << "libcpp: " << fs_libcxx() << "\n";
std::cout << "libc: " << fs_libc() << "\n";
std::cout << "maximum path length: " << fs_get_max_path() << "\n";
std::cout << "current working directory (CWD): " << cwd << "\n";
std::cout << "Compiler: " << fs_compiler() << "\n";
std::cout << "Optimized: " << fs_is_optimized() << " Trace: " << fs_trace << "\n";

std::cout << "Username: " << fs_get_username() << "\n";
std::cout << "Homedir: " << fs_get_homedir() << "\n";

if (fs_is_windows())
  std::cout << "Windows long paths enabled: " << fs_win32_long_paths_enabled() << "\n";

// doesn't work usefully on Cygwin
#if defined(ffilesystem_extra)

std::cout << "maximum path component length: " << fs_max_component(cwd) << "\n";
std::cout << "CPU arch: " << fs_cpu_arch() << "\n";
std::cout << "Shell: " << fs_get_shell() << "\n";
std::cout << "Hostname: " << fs_hostname() << "\n";

if(!fs_is_cygwin())
  std::cout << "CWD filesystem Type: " << fs_filesystem_type(cwd) << "\n";
#endif

// commented out  because some systems (Cygwin) silently exit with error code 0 on std::locale("") call
//std::cout << "Locale: " << fs_get_locale_name() << "\n";

#if defined(ffilesystem_extra)
if(fs_is_admin())
  std::cout << "NOTE: running as admin / sudo\n";
#endif

std::cout << std::endl;
// flush for CI etc.

while (true){

  std::string inp;

  std::cout << "Ffs> ";

  std::getline(std::cin, inp);

  // "\x04" is Ctrl-D on Windows.
  // EOF for non-Windows
  if (std::cin.eof() || inp ==
#ifdef __cpp_named_character_escapes
  "\x{04}"
#else
  "\x04"
#endif
  || inp == "q" || inp == "quit" || inp == "exit")
    break;

  // split variable inp on space-delimiters
  std::vector<std::string> args;
  {
    std::istringstream iss(inp);
    std::string token;
    while (iss >> token) {
      args.push_back(token);
    }
  }

  if(args.empty() || args.at(0).empty())
    continue;

  const std::vector<std::string>::size_type argc = args.size();

  switch (argc){
    case 1:
      if(no_arg(args.at(0)))
        continue;
      break;
    case 2:
      if(one_arg(args.at(0), args.at(1)))
        continue;
      break;
    case 3:
      if(two_arg(args.at(0), args.at(1), args.at(2)))
        continue;
      break;
    default:
      std::cerr << "ERROR: too many arguments\n";
      break;
  }

  std::cerr << "Unknown command or missing arguments for: " << inp << "\n";
}

return EXIT_SUCCESS;
}
