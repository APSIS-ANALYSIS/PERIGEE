#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include "FlowRateFactory.hpp"

int main(int argc, char *argv[])
{
  #if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  auto fr = FlowRateFactory::createFlowRate("test.txt");

  fr->print_info();

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
