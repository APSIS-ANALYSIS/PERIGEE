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
#include "NodalBC.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_vtu.hpp"

int main(int argc, char *argv[])
{
  const int nFunc = 7500;

  std::vector<std::string> fname { "wall_vol.vtu" };

  INodalBC * a = new NodalBC_3D_vtu( fname, nFunc );

  INodalBC * b = new NodalBC( fname, nFunc );

  for(int ii=0; ii<nFunc; ++ii)
  {
    if( a-> get_ID(ii) != b-> get_ID(ii) )
      std::cout<<"Error "<<a-> get_ID(ii)<<'\t'<<b-> get_ID(ii)<<'\n';
  }

  for(unsigned int ii=0; ii< a->get_num_dir_nodes(); ++ii)
  {
    if( a->get_dir_nodes(ii) != b->get_dir_nodes(ii) )
      std::cout<<"Error! \n";
  }

  for( unsigned int ii=0; ii< a->get_num_per_nodes(); ++ii )
  {
    if( a->get_per_slave_nodes(ii) != b->get_per_slave_nodes(ii) )
      std::cout<<"ERROR! \n";
    
    if( a->get_per_master_nodes(ii) != b->get_per_master_nodes(ii) )
      std::cout<<"ERROR! \n";
  }
  
  delete a; delete b;

  return EXIT_SUCCESS;
}

// EOF
