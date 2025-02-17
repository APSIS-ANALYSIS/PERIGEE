#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>

#include "AGlobal_Mesh_Info.hpp"
#include "Vis_Tools.hpp"

int main()
{
  const std::string epart_file = "epart.h5";
  const std::string part_file="postpart";
  const int rank = 0;
  auto GMIptr = new AGlobal_Mesh_Info(part_file,rank);
  
  const int nElem = GMIptr->get_nElem();

  std::vector<int> epart_map_1;
  VIS_T::read_epart( epart_file, nElem, epart_map_1 );

  auto epart_map_2 = VIS_T::read_epart(epart_file, nElem);

  VEC_T::print(epart_map_1);
  VEC_T::print(epart_map_2);

  if(VEC_T::is_equal(epart_map_1, epart_map_2, 0)) std::cout<<"good!\n";
  else std::cout<<"bad\n";

  delete GMIptr;
  return EXIT_SUCCESS;
}

// EOF
