#include "Math_Tools.hpp"
#include "Mesh_Tet.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Global_Part_Reload.hpp"
#include "Part_FEM_FSI.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_tet_outflow.hpp"
#include "NBC_Partition_MF.hpp"
#include "NBC_Partition_inflow_MF.hpp"
#include "EBC_Partition_outflow_MF.hpp"
#include "Vector_3.hpp"

int main( int argc, char * argv[] )
{

  // Define basic settings
  const int elemType = 601; // first order simplicial element

  // Input files
  std::string geo_file("./whole_vol.vtu");

  std::string geo_f_file("./fluid.vtu");
  std::string geo_s_file("./solid.vtu");

  std::string sur_s_file_interior_wall("./finter_solid.vtu");

  std::string sur_f_file_fro("./ffro_fluid.vtu");

  std::string sur_s_file_fro("./sfro_solid.vtu");  

  std::string sur_f_file_bac("./fbac_fluid.vtu");

  std::string sur_s_file_bac("./sbac_solid.vtu");

  std::string sur_f_file_lef("./flef_fluid.vtu");

  std::string sur_s_file_lef("./slef_solid.vtu");

  std::string sur_f_file_rig("./frig_fluid.vtu");

  std::string sur_s_file_rig("./srig_solid.vtu");

  std::string sur_f_file_top("./ftop_fluid.vtu");

  std::string sur_s_file_top("./stop_solid.vtu");

  std::cout<<"===== Command Line Arguments ====="<<std::endl;
  std::cout<<" -geo_file: "           <<geo_file           <<std::endl;
  std::cout<<" -geo_f_file: "         <<geo_f_file         <<std::endl;
  std::cout<<" -geo_s_file: "         <<geo_s_file         <<std::endl;
  std::cout<<" -sur_f_file_fro: "     <<sur_f_file_fro     <<std::endl;
  std::cout<<" -sur_s_file_fro: "     <<sur_s_file_fro     <<std::endl;
  std::cout<<" -sur_f_file_bac: "     <<sur_f_file_bac     <<std::endl;
  std::cout<<" -sur_s_file_bac: "     <<sur_s_file_bac     <<std::endl;
  std::cout<<" -sur_f_file_lef: "     <<sur_f_file_lef     <<std::endl;
  std::cout<<" -sur_s_file_lef: "     <<sur_s_file_lef     <<std::endl;
  std::cout<<" -sur_f_file_rig: "     <<sur_f_file_rig     <<std::endl;
  std::cout<<" -sur_s_file_rig: "     <<sur_s_file_rig     <<std::endl;
  std::cout<<" -sur_f_file_top: "     <<sur_f_file_top     <<std::endl;
  std::cout<<" -sur_s_file_top: "     <<sur_s_file_top     <<std::endl;
  std::cout<<"----------------------------------\n";
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<"===== Command Line Arguments ====="<<std::endl;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  // Check if the geometrical file exist on disk
  SYS_T::file_check(geo_file); std::cout<<geo_file<<" found. \n";

  SYS_T::file_check(geo_f_file); std::cout<<geo_f_file<<" found. \n";

  SYS_T::file_check(geo_s_file); std::cout<<geo_s_file<<" found. \n";

  SYS_T::file_check(sur_f_file_fro); std::cout<<sur_f_file_fro<<" found. \n";

  SYS_T::file_check(sur_s_file_fro); std::cout<<sur_s_file_fro<<" found. \n";

  SYS_T::file_check(sur_f_file_bac); std::cout<<sur_f_file_bac<<" found. \n";

  SYS_T::file_check(sur_s_file_bac); std::cout<<sur_s_file_bac<<" found. \n";

  SYS_T::file_check(sur_f_file_lef); std::cout<<sur_f_file_lef<<" found. \n";

  SYS_T::file_check(sur_s_file_lef); std::cout<<sur_s_file_lef<<" found. \n";

  SYS_T::file_check(sur_f_file_rig); std::cout<<sur_f_file_rig<<" found. \n";

  SYS_T::file_check(sur_s_file_rig); std::cout<<sur_s_file_rig<<" found. \n";

  SYS_T::file_check(sur_f_file_top); std::cout<<sur_f_file_top<<" found. \n";

  SYS_T::file_check(sur_s_file_top); std::cout<<sur_s_file_top<<" found. \n";

  SYS_T::file_check(sur_s_file_interior_wall); std::cout<<sur_s_file_interior_wall<<" found. \n";


  // Read the geometry file for the whole FSI domain for the velocity /
  // displacement field
  int nFunc_v, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid( geo_file, nFunc_v, nElem, ctrlPts, vecIEN );

  // Generate IEN
  IIEN * IEN_v = new IEN_FEM( nElem, vecIEN );

  // ----------------------------------------------------------------
  // Setup boundary conditions.
  
  // Physical ElemBC
  std::cout<<"4. Elem boundary for the implicit solver: \n";

  std::vector<std::string> ebclist {sur_f_file_fro, sur_f_file_bac, sur_f_file_lef, sur_f_file_rig, sur_f_file_top, 
                                    sur_s_file_fro, sur_s_file_bac, sur_s_file_lef, sur_s_file_rig, sur_s_file_top};

  //std::vector<std::string> ebclist {sur_s_file_fro};

  ElemBC * ebc = new ElemBC_3D( ebclist, 602 );

  ebc -> resetSurIEN_outwardnormal( IEN_v ); // assign outward orientation for surface

  //int nFunc_v_sur, nElem_sur;
  //std::vector<int> vecIEN_sur;
  //std::vector<double> ctrlPts_sur;

  //VTK_T::read_vtp_grid( sur_s_file_out, nFunc_v_sur, nElem_sur, ctrlPts_sur, vecIEN_sur );

  for (int ii=0; ii<static_cast<int>(ebclist.size()); ++ii)
  {  
    std::cout<<ebclist[ii]<<" the num of surface eLement:"<<ebc->get_num_cell(ii)<<std::endl;

    for (int ee=0; ee<ebc->get_num_cell(ii); ++ee)
    {
      std::cout<<"surface eLement:"<<ee<<std::endl;

      const double qua_r = 0.5;
      const double qua_s = 0.5;

      const double Nr[3] = { (2.0 * qua_r - 1.0) * (qua_r - 1.0),
          - 4.0 * qua_r * (qua_r - 1.0), qua_r * (2.0 * qua_r - 1.0) };
      const double Ns[3] = { (2.0 * qua_s - 1.0) * (qua_s - 1.0),
          - 4.0 * qua_s * (qua_s - 1.0), qua_s * (2.0 * qua_s - 1.0) };

      const double dNr[3] = { 4.0 * qua_r - 3.0, 
          - 8.0 * qua_r + 4.0, 4.0 * qua_r - 1.0 };
      const double dNs[3] = { 4.0 * qua_s - 3.0, 
          - 8.0 * qua_s + 4.0, 4.0 * qua_s - 1.0 };

      const double Rr[9] { 
      dNr[0] * Ns[0], dNr[2] * Ns[0], dNr[2] * Ns[2],
      dNr[0] * Ns[2], dNr[1] * Ns[0], dNr[2] * Ns[1],
      dNr[1] * Ns[2], dNr[0] * Ns[1], dNr[1] * Ns[1] };
      const double Rs[9] { 
      Nr[0] * dNs[0], Nr[2] * dNs[0], Nr[2] * dNs[2],
      Nr[0] * dNs[2], Nr[1] * dNs[0], Nr[2] * dNs[1],
      Nr[1] * dNs[2], Nr[0] * dNs[1], Nr[1] * dNs[1] };

      double ctrl_x[9] {};
      double ctrl_y[9] {};
      double ctrl_z[9] {};

      for (int jj=0; jj<9; ++jj)
      {
        ctrl_x[jj] = ebc->get_pt_xyz(ii, ebc->get_ien(ii, ee, jj), 0);
        ctrl_y[jj] = ebc->get_pt_xyz(ii, ebc->get_ien(ii, ee, jj), 1);
        ctrl_z[jj] = ebc->get_pt_xyz(ii, ebc->get_ien(ii, ee, jj), 2);

        std::cout<<jj<<": "<<ebc->get_ien(ii, ee, jj)<<"\t"<<ebc->get_global_node(ii, ebc->get_ien(ii, ee, jj))<<"\t("<<ctrl_x[jj]<<","<<ctrl_y[jj]<<","<<ctrl_z[jj]<<")"<< std::endl;
      }

      Vector_3 dx_dr( 0.0, 0.0, 0.0 );
      Vector_3 dx_ds( 0.0, 0.0, 0.0 );
      Vector_3 un( 0.0, 0.0, 0.0 );

      for( int nn=0; nn<9; ++nn )
      {
        dx_dr += Vector_3( ctrl_x[nn] * Rr[nn], ctrl_y[nn] * Rr[nn], ctrl_z[nn] * Rr[nn] );
        dx_ds += Vector_3( ctrl_x[nn] * Rs[nn], ctrl_y[nn] * Rs[nn], ctrl_z[nn] * Rs[nn] );
      }

      un = VEC3_T::cross_product( dx_dr, dx_ds );

      un.normalize();

      std::cout<<"surface element "<<ee<<" out_normal: ("<<un.x()<<","<<un.y()<<","<<un.z()<<")"<<std::endl;
    }

  }


  PetscFinalize();

  return EXIT_SUCCESS;
}

// EOF
