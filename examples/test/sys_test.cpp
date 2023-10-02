#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "IEN_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"

int main(int argc, char *argv[])
{
  const int elemType = 601;

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

  std::cout<<"===== Command Line Arguments ====="<<std::endl;
  std::cout<<" -num_inlet: "          <<num_inlet          <<std::endl;
  std::cout<<" -num_outlet: "         <<num_outlet         <<std::endl;
  std::cout<<" -geo_file: "           <<geo_file           <<std::endl;
  std::cout<<" -geo_f_file: "         <<geo_f_file         <<std::endl;
  std::cout<<" -geo_s_file: "         <<geo_s_file         <<std::endl;
  std::cout<<" -sur_f_file_wall: "    <<sur_f_file_wall    <<std::endl;
  std::cout<<" -sur_s_file_wall: "    <<sur_s_file_wall    <<std::endl;
  std::cout<<" -sur_s_file_int_wall: "<<sur_s_file_interior_wall <<std::endl;
  std::cout<<" -sur_f_file_in_base: " <<sur_f_file_in_base <<std::endl;
  std::cout<<" -sur_f_file_out_base: "<<sur_f_file_out_base<<std::endl;
  std::cout<<" -sur_s_file_in_base: " <<sur_s_file_in_base <<std::endl;
  std::cout<<" -sur_s_file_out_base: "<<sur_s_file_out_base<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  // Check if the geometrical file exist on disk
  SYS_T::file_check(geo_file); std::cout<<geo_file<<" found. \n";

  SYS_T::file_check(geo_f_file); std::cout<<geo_f_file<<" found. \n";

  SYS_T::file_check(geo_s_file); std::cout<<geo_s_file<<" found. \n";

  SYS_T::file_check(sur_f_file_wall); std::cout<<sur_f_file_wall<<" found. \n";

  SYS_T::file_check(sur_s_file_wall); std::cout<<sur_s_file_wall<<" found. \n";

  SYS_T::file_check(sur_s_file_interior_wall); std::cout<<sur_s_file_interior_wall<<" found. \n";

  std::vector< std::string > sur_f_file_in(  num_inlet ) , sur_s_file_in(  num_inlet );
  std::vector< std::string > sur_f_file_out( num_outlet ), sur_s_file_out( num_outlet );

  for(int ii=0; ii<num_inlet; ++ii)
  {
    sur_f_file_in[ii] = SYS_T::gen_capfile_name( sur_f_file_in_base, ii, ".vtp" );
    sur_s_file_in[ii] = SYS_T::gen_capfile_name( sur_s_file_in_base, ii, ".vtp" );

    SYS_T::file_check( sur_f_file_in[ii] );
    std::cout<<sur_f_file_in[ii]<<" found. \n";
    SYS_T::file_check( sur_s_file_in[ii] );
    std::cout<<sur_s_file_in[ii]<<" found. \n";
  }

  for(int ii=0; ii<num_outlet; ++ii)
  {
    sur_f_file_out[ii] = SYS_T::gen_capfile_name( sur_f_file_out_base, ii, ".vtp" );
    sur_s_file_out[ii] = SYS_T::gen_capfile_name( sur_s_file_out_base, ii, ".vtp" );

    SYS_T::file_check( sur_f_file_out[ii] );
    std::cout<<sur_f_file_out[ii]<<" found. \n";
    SYS_T::file_check( sur_s_file_out[ii] );
    std::cout<<sur_s_file_out[ii]<<" found. \n";
  }

  // Read the geometry file for the whole FSI domain for the velocity /
  // displacement field
  int nFunc_v, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid( geo_file, nFunc_v, nElem, ctrlPts, vecIEN );

  // Generate IEN
  IIEN * IEN_v = new IEN_FEM( nElem, vecIEN );

  // InflowBC info
  std::cout<<"1. Inflow cap surfaces: \n";
  std::vector<Vector_3> inlet_outvec( num_inlet );
  for(int ii=0; ii<num_inlet; ++ii)
    inlet_outvec[ii] = HEX_T::get_out_normal( sur_f_file_in[ii], ctrlPts, IEN_v );

  INodalBC * InFBC = new NodalBC_3D_inflow( sur_f_file_in, sur_f_file_wall, nFunc_v, inlet_outvec, elemType);

  InFBC -> resetSurIEN_outwardnormal( IEN_v ); // assign outward orientation for inlet surface mesh

  for (int ii=0; ii<static_cast<int>(sur_f_file_in.size()); ++ii)
  {  
    std::cout<<sur_f_file_in[ii]<<" the num of lumen inlet eLement:"<<InFBC->get_num_cell(ii)<<std::endl;

    const std::vector<double> intBasis = InFBC->get_intNA(ii);

    for (int ee=0; ee<InFBC->get_num_cell(ii); ++ee)
    {
      std::cout<<"the lumen inlet eLement :"<<ee<<std::endl;

      const double qua_r = 0.5;
      const double qua_s = 0.5;

      const double Rr[4] { qua_s - 1.0, 1.0 - qua_s, qua_s, -qua_s };
      const double Rs[4] { qua_r - 1.0, -qua_r, qua_r, 1.0 - qua_r };

      double ctrl_x[4] {};
      double ctrl_y[4] {};
      double ctrl_z[4] {};

      for (int jj=0; jj<4; ++jj)
      {
        ctrl_x[jj] = InFBC->get_pt_xyz(ii, InFBC->get_ien(ii, ee, jj), 0);
        ctrl_y[jj] = InFBC->get_pt_xyz(ii, InFBC->get_ien(ii, ee, jj), 1);
        ctrl_z[jj] = InFBC->get_pt_xyz(ii, InFBC->get_ien(ii, ee, jj), 2);

        std::cout<<jj<<": "<<InFBC->get_ien(ii, ee, jj)<<"\t"<<InFBC->get_global_node(ii, InFBC->get_ien(ii, ee, jj))
                 <<"\t("<<ctrl_x[jj]<<","<<ctrl_y[jj]<<","<<ctrl_z[jj]<<")"
                 <<"\t"<< intBasis[InFBC->get_ien(ii, ee, jj)]<<std::endl;
      }

      Vector_3 dx_dr( 0.0, 0.0, 0.0 );
      Vector_3 dx_ds( 0.0, 0.0, 0.0 );
      Vector_3 un( 0.0, 0.0, 0.0 );

      for( int nn=0; nn<4; ++nn )
      {
        dx_dr += Vector_3( ctrl_x[nn] * Rr[nn], ctrl_y[nn] * Rr[nn], ctrl_z[nn] * Rr[nn] );
        dx_ds += Vector_3( ctrl_x[nn] * Rs[nn], ctrl_y[nn] * Rs[nn], ctrl_z[nn] * Rs[nn] );
      }

      un = Vec3::cross_product( dx_dr, dx_ds );

      un.normalize();

      std::cout<<"the lumen inlet eLement "<<ee<<" out_normal: ("<<un.x()<<","<<un.y()<<","<<un.z()<<")"<<std::endl;
      std::cout<<"-------------------------------------------------"<<std::endl;
    }

    //for (auto &ite: intBasis) std::cout<<ite<<std::endl;

  }

  // Physical ElemBC
  cout<<"2. Elem boundary for the implicit solver: \n";
  std::vector< Vector_3 > outlet_outvec( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
    outlet_outvec[ii] = HEX_T::get_out_normal( sur_f_file_out[ii], ctrlPts, IEN_v );  

  ElemBC * ebc = new ElemBC_3D_outflow( sur_f_file_out, outlet_outvec, elemType);

  ebc -> resetSurIEN_outwardnormal( IEN_v ); // assign outward orientation for outlet surface mesh 
  
  for (int ii=0; ii<static_cast<int>(sur_f_file_out.size()); ++ii)
  {  
    std::cout<<sur_f_file_out[ii]<<" the num of lumen outlet eLement:"<<ebc->get_num_cell(ii)<<std::endl;

    const std::vector<double> intBasis = ebc->get_intNA(ii);

    for (int ee=0; ee<ebc->get_num_cell(ii); ++ee)
    {
      std::cout<<"the lumen outlet eLement:"<<ee<<std::endl;

      const double qua_r = 0.31;
      const double qua_s = 0.65;

      const double Rr[4] { qua_s - 1.0, 1.0 - qua_s, qua_s, -qua_s };
      const double Rs[4] { qua_r - 1.0, -qua_r, qua_r, 1.0 - qua_r };

      double ctrl_x[4] {};
      double ctrl_y[4] {};
      double ctrl_z[4] {};

      for (int jj=0; jj<4; ++jj)
      {
        ctrl_x[jj] = ebc->get_pt_xyz(ii, ebc->get_ien(ii, ee, jj), 0);
        ctrl_y[jj] = ebc->get_pt_xyz(ii, ebc->get_ien(ii, ee, jj), 1);
        ctrl_z[jj] = ebc->get_pt_xyz(ii, ebc->get_ien(ii, ee, jj), 2);

        std::cout<<jj<<": "<<ebc->get_ien(ii, ee, jj)<<"\t"<<ebc->get_global_node(ii, ebc->get_ien(ii, ee, jj))
                 <<"\t("<<ctrl_x[jj]<<","<<ctrl_y[jj]<<","<<ctrl_z[jj]<<")"
                 <<"\t"<< intBasis[ebc->get_ien(ii, ee, jj)]<<std::endl;
      }

      Vector_3 dx_dr( 0.0, 0.0, 0.0 );
      Vector_3 dx_ds( 0.0, 0.0, 0.0 );
      Vector_3 un( 0.0, 0.0, 0.0 );

      for( int nn=0; nn<4; ++nn )
      {
        dx_dr += Vector_3( ctrl_x[nn] * Rr[nn], ctrl_y[nn] * Rr[nn], ctrl_z[nn] * Rr[nn] );
        dx_ds += Vector_3( ctrl_x[nn] * Rs[nn], ctrl_y[nn] * Rs[nn], ctrl_z[nn] * Rs[nn] );
      }

      un = Vec3::cross_product( dx_dr, dx_ds );

      un.normalize();

      std::cout<<"the lumen outlet eLement "<<ee<<" out_normal: ("<<un.x()<<","<<un.y()<<","<<un.z()<<")"<<std::endl;
      std::cout<<"-------------------------------------------------"<<std::endl;
    }

  }

  delete IEN_v; delete InFBC; delete ebc;

  PetscFinalize();

  return EXIT_SUCCESS;
}

// EOF
