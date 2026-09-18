#include <string>
#include <string_view>

#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
  using namespace boost::ut;

  "username"_test = [] {
    std::string user = fs_get_username();
    expect(!user.empty());
    std::cout << "Username " << user << "\n";
  };

  "environment"_test = [] {
    expect(!fs_get_cwd().empty());

    std::string pdir = fs_get_profile_dir();
    expect(!pdir.empty());
    std::cout << "Profile directory " << pdir << "\n";

    std::string p = fs_get_homedir();
    expect(!p.empty());
    std::cout << "Home directory " << p << "\n";
    expect(fs_is_dir(p));

    // NOTE: profiledir does not always (but may) equal homedir, for example when root user.

    expect(eq(fs_expanduser("~"), p));

    // --- tempdir
    auto t = fs_get_tempdir();
    expect(!t.empty());
    std::cout << "Temp directory " << t << "\n";
    expect(fs_is_dir(t));
  };

  "getenv"_test = [] {
    const std::string name = fs_is_windows() ? "USERPROFILE" : "HOME";
    std::string in_name = name + "-walkoff-end-of-buffer";
    std::string_view nonnull_name(in_name.data(), name.size());
    expect(nonnull_name.back() != '\0') << "Environment variable name should not be null-terminated in test";

    auto h = fs_getenv(fs_is_windows() ? "USERPROFILE" : "HOME");
    expect(h.has_value() >> fatal) << "Environment variable HOME or USERPROFILE should be set";

    h = fs_getenv(nonnull_name);
    expect(h.has_value()) << "Environment variable " << nonnull_name << " should be set in test";
  };

  "setenv"_test = [] {
    const std::string k = "FORTtest";
    const std::string v = "FORTvalue";
    fs_setenv(k, v);

    auto e = fs_getenv(k);
    expect(e.has_value() >> fatal) << "Environment variable " << k << " not set";
    expect(eq(e.value(), v));

    fs_setenv(k, ""); // unset
    e = fs_getenv(k);
    expect(!e.has_value()) << "Environment variable " << k << " should be unset";

    std::string in_k = k + "-walkoff-end-of-buffer";
    std::string_view nonnull_k(in_k.data(), k.size());
    expect(nonnull_k.back() != '\0' >> fatal) << "Environment variable name should not be null-terminated in test";

    expect(fs_setenv(nonnull_k, v) >> fatal) << "Failed to set environment variable with non-null-terminated name";
    e = fs_getenv(k);
    expect(e.has_value() >> fatal) << "Environment variable " << k << " not set with non-null-terminated name";
    expect(eq(e.value(), v)) << "Environment variable " << k << " has wrong value with non-null-terminated name";
  };

  "no_environment"_test = [] {
    expect(fs_setenv("PATH", "") >> fatal);
    expect(fs_setenv("HOME", "") >> fatal);
    expect(fs_setenv("XDG_CONFIG_HOME", "") >> fatal);
    expect(fs_setenv("TMPDIR", "") >> fatal);
    if (fs_is_windows()) {
      expect(fs_setenv("USERPROFILE", "") >> fatal);
    }

    std::string pdir = fs_get_profile_dir();
    expect(!pdir.empty());
    std::cout << "Profile directory " << pdir << "\n";

    auto t = fs_get_tempdir();
    expect(!t.empty());
    std::cout << "Temp directory " << t << "\n";
    expect(fs_is_dir(t));
  };

  "no_home"_test = [] {
    expect(fs_setenv("HOME", "") >> fatal);
    expect(fs_setenv("XDG_CONFIG_HOME", "") >> fatal);
    if (fs_is_windows()) {
      expect(fs_setenv("USERPROFILE", "") >> fatal);
    }

    auto h = fs_getenv(fs_is_windows() ? "USERPROFILE" : "HOME");
    expect(!h.has_value()) << "Environment variable HOME or USERPROFILE should not be set in test";

    std::string p = fs_get_homedir();
    expect(!p.empty());
    std::cout << "Home directory " << p << "\n";
    expect(fs_is_dir(p));

    // NOTE: profiledir does not always (but may) equal homedir, for example when root user.

    expect(eq(fs_expanduser("~"), p));
  };
}
