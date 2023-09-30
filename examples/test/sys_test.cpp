#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"

int main(int argc, char *argv[])
{ 
  std::string myname {};

  SYS_T::get_option_string("-name", myname, argc, argv);

  std::cout << myname << std::endl;

  return EXIT_SUCCESS;
}

// EOF
