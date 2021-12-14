// ============================================================================
// preprocess_fsi_p1p0.cpp
// ----------------------------------------------------------------------------
// This is the preprocess code for handling FSI 3D meshes that uses tet for the
// geometry and discontinuous pressure space.
//
// Author: Ju Liu
// Date: Dec. 13 2021
// ============================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Global_Part_Reload.hpp"
#include "Part_Tet_FSI.hpp"
#include "NodalBC_3D_vtu.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_tet_outflow.hpp"
#include "NBC_Partition.hpp"
#include "NBC_Partition_inflow.hpp"
#include "EBC_Partition_outflow.hpp"

int main( int argc, char * argv[] )
{
  // Remove previously existing hdf5 files
  SYS_T::execute("rm -rf apart");
  SYS_T::execute("mkdir apart");

  // Define basic settings
  const int dofNum = 7;     // degree-of-freedom for the physical problem
  const int dofMat = 4;     // degree-of-freedom in the matrix problem
  const int elemType = 501; // first order simplicial element

  // Input files
  std::string geo_file("./whole_vol.vtu");

  std::string geo_f_file("./lumen_vol.vtu");
  std::string geo_s_file("./tissue_vol.vtu");

  std::string sur_f_file_wall("./lumen_wall_vol.vtp");
  std::string sur_f_file_in_base( "./lumen_inlet_vol_" );
  std::string sur_f_file_out_base("./lumen_outlet_vol_");

  std::string sur_s_file_interior_wall("./tissue_interior_wall_vol.vtp");

  std::string sur_s_file_wall("./tissue_wall_vol.vtp");
  std::string sur_s_file_in_base( "./tissue_inlet_vol_" );
  std::string sur_s_file_out_base("./tissue_outlet_vol_");

  int num_outlet = 1, num_inlet = 1;

  const std::string part_file("./apart/part");

  // fsiBC_type : 0 deformable wall, 1 rigid wall
  int fsiBC_type = 0;

  // ringBC_type : 0 fully clamped, 1 in-plane motion allowed
  int ringBC_type = 0;

  // Mesh partition setting
  int cpu_size = 1;
  int in_ncommon = 2;
  const bool isDualGraph = true;

  bool isReload = false;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor needs to be run in serial.\n");





  PetscFinalize();
  cout<<"===> Preprocessing completes successfully!\n";
  return EXIT_SUCCESS;
}



// EOF
