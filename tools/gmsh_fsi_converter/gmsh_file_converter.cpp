// ============================================================================
// to be added
//
// Author: Jiayi Huang and Ju Liu
// Date: Sept. 20 2022
// ============================================================================
#include "Tet_Tools.hpp"
#include "ElemBC_3D_tet_wall.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  // const int elemType = 502;
  // const int num_outlet = 0;
  std::string geo_file("./mesh/whole_vol.vtu");
  // std::string sur_file_in("./mesh/inflow_vol.vtu");
  // std::string sur_file_wall("./mesh/wall_vol.vtu");
  // std::string sur_file_out_base("./mesh/outflow_vol_");

  // std::vector< std::string > sur_file_out;
  // sur_file_out.resize( num_outlet );

  // for(int ii=0; ii<num_outlet; ++ii)
  // {
  //   if(elemType == 501 )
  //     sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtp" );
  //   else
  //     sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtu" );

  //   SYS_T::file_check(sur_file_out[ii]);
  //   cout<<sur_file_out[ii]<<" found. \n";
  // }

  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);
  printf("There are %d points in the whole_vol.vtu\n",nFunc);
  printf("There are %d cells in the whole_vol.vtu\n",nElem);

  PetscFinalize();
  return 0;
}

//EOF
