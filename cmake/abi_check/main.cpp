#include <iostream>
#include <cstdlib>

extern "C" int addone(int);


int main() {

  if(addone(2) != 3) {
    std::cerr << "2 + 1 != " << addone(2) << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
