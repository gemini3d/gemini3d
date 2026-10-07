// Audit addition 2026-09-16. Apache-2.0.
#include "cli_utils.h"
#include <cassert>
#include <climits>
#include <vector>
#include <iostream>
#include <type_traits>
// Remains active in Release builds, unlike assert().
void require(bool condition) { if (!condition) throw std::runtime_error("CLI regression"); }
bool parse(std::vector<std::string> args, int& x, int& y) {
 std::vector<char*> argv;
 for (auto& s : args) argv.push_back(s.data());
 params p{}; bool help=false; std::string error;
 return gemini_cli::parse_options(static_cast<int>(argv.size()),argv.data(),p,x,y,help,error);
}
int main() {
 static_assert(std::is_same_v<decltype(&get_config_vars_C),void(*)(void**,int*,int*,double*,double*)>);
 int x=-1,y=-1;
 require(parse({"gemini",".","-manual_grid","2","3","-dryrun"},x,y));
 require(x==2 && y==3);
 for (auto bad : {"0","-1","2x","2147483648","","1.5"})
   require(!parse({"gemini",".","-manual_grid",bad,"2"},x,y));
 require(!parse({"gemini",".","-manual_grid"},x,y));
 require(!parse({"gemini",".","-manual_grid","2"},x,y));
 require(!parse({"gemini",".","-unknown"},x,y));
 char out[1000]; require(gemini_cli::copy_outdir(std::string(999,'x'),out,sizeof(out)));
 require(out[999]=='\0'); require(!gemini_cli::copy_outdir(std::string(1000,'x'),out,sizeof(out)));
 const auto cells=gemini_cli::cell_count(500,500,500);
 require(gemini_cli::array_bytes(cells,35)==static_cast<std::size_t>(504)*504*504*35*8);
 bool caught=false;
 try { (void)gemini_cli::cell_count(INT_MAX,INT_MAX,INT_MAX); } catch(const std::overflow_error&) {caught=true;}
 require(caught); caught=false;
 try { (void)gemini_cli::cell_count(0,2,2); } catch(const std::invalid_argument&) {caught=true;}
 require(caught);
 std::cout << "PASS CLI, ABI declaration and allocation boundaries\n";
}
