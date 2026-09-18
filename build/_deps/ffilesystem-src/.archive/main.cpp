#include <iostream>
#include <cstdlib>
#include <string>

#ifdef _MSC_VER
#include <crtdbg.h>
#endif

#include "ffilesystem.h"


int main(int argc, char* argv[]){
#ifdef _MSC_VER
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
  _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
  _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
#endif

  if (argc == 1) {
    std::cerr << "fs_cli <function_name> [<arg1> ...]\n";
    return EXIT_FAILURE;
  }

  std::string arg1 = argv[1];
  std::string arg2, arg3;
  // default initializer for std::string is empty string
  if(argc > 2)
    arg2 = argv[2];
  if(argc > 3)
    arg3 = argv[3];

  if (arg1 == "pathsep")
    std::cout << fs_pathsep() << "\n";
  else if (arg1 == "touch"){
    std::cout << "touch " << arg2 << "\n";
    Ffs::touch(arg2);
  }
  else if (arg1 == "remove")
    std::cout << "remove " << arg2 << " " << Ffs::remove(arg2) << "\n";
  else if (arg1 == "expanduser")
    std::cout << Ffs::expanduser(arg2) << "\n";
  else if (arg1 == "canonical")
    std::cout << Ffs::canonical(arg2, false, false).value() << "\n";
  else if (arg1 == "resolve")
    std::cout << Ffs::resolve(arg2, false, false).value() << "\n";
  else if (arg1 == "compiler")
    std::cout << fs_compiler() << "\n";
  else if (arg1 == "chdir" || arg1 == "set_cwd"){
    std::cout << "cwd: " << Ffs::get_cwd() << "\n";
    Ffs::chdir(arg2);
    std::cout << "new cwd: " << Ffs::get_cwd() << "\n";
  }
  else if (arg1 == "perm")
    std::cout << Ffs::get_permissions(arg2) << "\n";
  else if (arg1 == "set_perm" && argc == 6){
    if(fs_is_windows())
      std::cerr << "chmod is not supported on Windows\n";

    int r = std::stoi(arg3);
    int w = std::stoi(argv[4]);
    int x = std::stoi(argv[5]);
    std::cout << "chmod " << Ffs::get_permissions(arg2) << " " << arg2 << " => ";
    Ffs::set_permissions(arg2, r, w, x);
    std::cout << Ffs::get_permissions(arg2) << " " << arg2 << "\n";
  }
  else if (arg1 == "cpp")
    std::cout << fs_cpp() << "\n";
  else if (arg1 == "lang")
    std::cout << fs_lang() << "\n";
  else if (arg1 == "homedir")
    std::cout << Ffs::get_homedir() << "\n";
  else if (arg1 == "tempdir")
    std::cout << Ffs::get_tempdir() << "\n";
  else if (arg1 == "lib_path")
    std::cout << Ffs::lib_path() << "\n";
  else if (arg1 == "exe_path")
    std::cout << Ffs::exe_path() << "\n";
  else if (arg1 == "is_admin")
    std::cout << fs_is_admin() << "\n";
  else if (arg1 == "is_bsd")
    std::cout << fs_is_bsd() << "\n";
  else if (arg1 == "is_linux")
    std::cout << fs_is_linux() << "\n";
  else if (arg1 == "is_macos")
    std::cout << fs_is_macos() << "\n";
  else if (arg1 == "is_unix")
    std::cout << fs_is_unix() << "\n";
  else if (arg1 == "is_windows")
    std::cout << fs_is_windows() << "\n";
  else if (arg1 == "is_wsl")
    std::cout << fs_is_wsl() << "\n";
  else if (arg1 == "is_mingw")
    std::cout << fs_is_mingw() << "\n";
  else if (arg1 == "is_cygwin")
    std::cout << fs_is_cygwin() << "\n";
  else if (arg1 == "parent")
    std::cout << Ffs::parent(arg2) << "\n";
  else if (arg1 == "root")
    std::cout << Ffs::root(arg2) << "\n";
  else if (arg1 == "file_size")
    std::cout << Ffs::file_size(arg2) << "\n";
  else if (arg1 == "exists")
    std::cout << Ffs::exists(arg2) << "\n";
  else if (arg1 == "is_dir")
    std::cout << Ffs::is_dir(arg2) << "\n";
  else if (arg1 == "is_subdir")
    std::cout << Ffs::is_subdir(arg2, arg3) << "\n";
  else if (arg1 == "is_exe")
    std::cout << Ffs::is_exe(arg2) << "\n";
  else if (arg1 == "which")
    std::cout << Ffs::which(arg2) << "\n";
  else if (arg1 == "shortname")
    std::cout << Ffs::shortname(arg2) << "\n";
  else if (arg1 == "longname")
    std::cout << Ffs::longname(arg2) << "\n";
  else if (arg1 == "is_char")
    std::cout << Ffs::is_char_device(arg2) << "\n";
  else if (arg1 == "is_file")
    std::cout << Ffs::is_file(arg2) << "\n";
  else if (arg1 == "same")
    std::cout << Ffs::equivalent(arg2, arg3) << "\n";
  else if (arg1 == "create_symlink"){
    std::cout << "create_symlink " << arg2 << " <= " << arg3 << "\n";
    Ffs::create_symlink(arg2, arg3);
  } else if (arg1 == "is_symlink")
    std::cout << Ffs::is_symlink(arg2) << "\n";

  else if (arg1 == "copy_file"){
    std::cout << "copy_file " << arg2 << " => " << arg3 << "\n";
    Ffs::copy_file(arg2, arg3, false);
  } else if (arg1 == "mkdir"){
    std::cout << "mkdir " << arg2 << "\n";
    Ffs::mkdir(arg2);
  }
  else if (arg1 == "relative_to")
    std::cout << Ffs::relative_to(arg2, arg3) << "\n";
  else if (arg1 == "normal")
    std::cout << Ffs::normal(arg2) << "\n";
  else{
    std::cerr << "fs_cli <function_name> [<arg1> ...]\n";
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;

}
