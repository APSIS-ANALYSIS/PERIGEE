#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"
#include "Mesh_Tet.hpp"
#include "IEN_FEM.hpp"
#include "Matrix_double_3by3_Array.hpp"
#include "Matrix_double_6by6_Array.hpp"
#include "VTK_Tools.hpp"

int main(int argc, char *argv[])
{
  const std::string vtufilename = "whole_vol.vtu";
  const std::string vtpfilename = "outflow_vol_000.vtp";

  const std::string filename = vtufilename;

  std::cout << "Filename: " << filename << std::endl;

  // old version
  int numpts_old, numcels_old;
  std::vector<double> pts_old;
  std::vector<int> ien_old;

  VTK_T::read_vtu_grid(filename, numpts_old, numcels_old, pts_old, ien_old);

  std::cout << "Old reading succeeded." << std::endl; 

  // new version
  int numpts_new, numcels_new;
  std::vector<double> pts_new;
  std::vector<int> ien_new;

  VTK_T::read_grid(filename, numpts_new, numcels_new, pts_new, ien_new);

  std::cout << "New reading succeeded." << std::endl; 

  int numpts_newnew = VTK_T::read_num_pt(filename);
  int numcels_newnew = VTK_T::read_num_cl(filename);

  std::cout << "New new reading succeeded." << std::endl; 

  // Comparing
  std::cout << "Compare: " << std::endl;
  std::cout << "Old numpts: " << numpts_old << std::endl;
  std::cout << "New numpts: " << numpts_new << std::endl;
  std::cout << "New new numpts: " << numpts_newnew << std::endl;
  std::cout << std::endl;

  std::cout << "Old numcels: " << numcels_old << std::endl;
  std::cout << "New numcels: " << numcels_new << std::endl;
  std::cout << "New new numcels: " << numcels_newnew << std::endl;
  std::cout << std::endl;

  std::cout << "Old pts: (first 100)" << std::endl;
  for (int it{0}; it != 100; ++it) {
    std::cout << pts_old[it] << " ";
  }
  std::cout << std::endl;

  std::cout << "New pts: (first 100)" << std::endl;
  for (int it{0}; it != 100; ++it) {
    std::cout << pts_new[it] << " ";
  }
  std::cout << std::endl;

  std::cout << "Old ien: (first 100)" << std::endl;
  for (int it{0}; it != 100; ++it) {
    std::cout << ien_old[it] << " ";
  }
  std::cout << std::endl;

  std::cout << "New ien: (first 100)" << std::endl;
  for (int it{0}; it != 100; ++it) {
    std::cout << ien_new[it] << " ";
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

// EOF
