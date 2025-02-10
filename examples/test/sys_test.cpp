#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>

#include "GenBCFactory.hpp"

int main()
{
  auto gbc = GenBCFactory::createGenBC("test.txt", 0.3, 0.14, 1, 200);

  gbc -> print_info();

  return EXIT_SUCCESS;
}

// EOF
