#include "Sys_Tools.hpp"

double SYS_T::gen_randomD_closed(const double &min, const double &max)
{
  return ( rand() % 1000001 ) * 1.0e-6 * (max - min) + min;
}

double SYS_T::gen_randomD_open(const double &min, const double &max)
{
  return ( rand() % 999998 + 1 ) * 1.0e-6 * (max - min) + min;
}

int SYS_T::gen_randomI_closed(const int &min, const int &max)
{
  return ( rand() % (max - min + 1)) + min;
}

// EOF
