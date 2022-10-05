#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"


int main(int argc, char *argv[])
{
  std::vector<std::string> a;
  a.push_back("hello a\n");
  a.push_back("cmame \t");
  a.push_back("july 31th \n");
  a.push_back("aug 31th \n");
  double jj = 0;
  // for(auto &out : a) out = "hello\n";
  for (int ii = 0; ii < 4; ++ii)
  {
    if (ii > 2)
    {
      jj = ii + 0.5;
      a.push_back("\n std::to_string: ");
      a.push_back(std::to_string(jj));
      a.push_back("\n   std::to_string: ");
      a.push_back(std::to_string(jj));
    }
    else
    {
      a.push_back("\n std::to_string: ");
      a.push_back(std::to_string(ii));
      a.push_back("\n   std::to_string: ");
      a.push_back(std::to_string(ii));
    }
  }

  for (auto out : a)
    std::cout << out;

  printf("\n");
  SymmMatrix_3x3 A;

  A.gen_id();

  // A.print_in_row();

  return EXIT_SUCCESS;
}

// EOF
