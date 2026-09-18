#include <boost/filesystem.hpp>
#include <filesystem>

#include <string>
#include <iostream>
#include <cstdlib>


int main()
{

#ifndef _WIN32
std::cerr << "This example is for Windows.\n";
#endif

std::string s1 = R"(C:\a\ball.text)";
std::string s2 = R"(\\?\C:\a\ball.text)";
std::string s3 = R"(\\.\C:\a\ball.text)";

for (const auto& s : {s1, s2, s3}) {
    boost::filesystem::path bp(s);
    std::filesystem::path p(s);

    std::cout << "Boost: " << bp << "    std::filesystem: " << p << "\n";

    std::cout << "  root_name: " << bp.root_name() << "    " << p.root_name() << "\n";
    std::cout << "  root_path: " << bp.root_path() << "    " << p.root_path() << "\n";
    std::cout << "  root_directory: " << bp.root_directory() << "    " << p.root_directory() << "\n";
    std::cout << "---------------------" << "\n";
}

return EXIT_SUCCESS;
}
