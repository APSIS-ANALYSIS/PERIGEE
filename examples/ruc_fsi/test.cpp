// Test the new TET_T functions

#include "ElemBC_3D_tet_wall.hpp"
#include "NodalBC_3D_wall.hpp"
#include "NodalBC_3D_CMM.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const int elemType = 502;
  const int num_outlet = 1;
  std::string geo_file("./whole_vol.vtu");
  std::string sur_file_in("./inflow_vol.vtu");
  std::string sur_file_wall("./wall_vol.vtu");
  std::string sur_file_out_base("./outflow_vol_");

  std::vector< std::string > sur_file_out;
  sur_file_out.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    if(elemType == 501 )
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtp" );
    else
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtu" );

    SYS_T::file_check(sur_file_out[ii]);
    cout<<sur_file_out[ii]<<" found. \n";
  }

  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);


  INodalBC * wallBC = new NodalBC_3D_wall( sur_file_in,
      sur_file_wall, sur_file_out, nFunc, elemType );

  wallBC -> print_info();

  INodalBC * nbc = new NodalBC_3D_CMM( wallBC, nFunc );

  nbc -> print_info();

  delete wallBC;
  delete nbc;
  PetscFinalize();
  return 0;
}

//EOF
