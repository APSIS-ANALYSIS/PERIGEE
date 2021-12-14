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
  const int num_fields = 2; // Two fields : pressure + velocity/displacement

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

  SYS_T::GetOptionInt(   "-cpu_size",            cpu_size);
  SYS_T::GetOptionInt(   "-in_ncommon",          in_ncommon);
  SYS_T::GetOptionInt(   "-fsiBC_type",          fsiBC_type);
  SYS_T::GetOptionInt(   "-ringBC_type",         ringBC_type);
  SYS_T::GetOptionInt(   "-num_outlet",          num_outlet);
  SYS_T::GetOptionInt(   "-num_inlet",           num_inlet);
  SYS_T::GetOptionString("-geo_file",            geo_file);
  SYS_T::GetOptionString("-geo_f_file",          geo_f_file);
  SYS_T::GetOptionString("-geo_s_file",          geo_s_file);
  SYS_T::GetOptionString("-sur_f_file_wall",     sur_f_file_wall);
  SYS_T::GetOptionString("-sur_s_file_wall",     sur_s_file_wall);
  SYS_T::GetOptionString("-sur_f_file_in_base",  sur_f_file_in_base);
  SYS_T::GetOptionString("-sur_f_file_out_base", sur_f_file_out_base);
  SYS_T::GetOptionString("-sur_s_file_in_base",  sur_s_file_in_base);
  SYS_T::GetOptionString("-sur_s_file_out_base", sur_s_file_out_base);
  SYS_T::GetOptionBool("-isReload",              isReload);

  SYS_T::print_fatal_if( fsiBC_type != 0 && fsiBC_type != 1 && fsiBC_type != 2, "Error: fsiBC_type should be 0, 1, or 2.\n" );
  SYS_T::print_fatal_if( ringBC_type != 0 && ringBC_type != 1, "Error: ringBC_type should be 0 or 1.\n" );

  std::cout<<"===== Command Line Arguments ====="<<std::endl;
  std::cout<<" -fsiBC_type: "         <<fsiBC_type         <<std::endl;
  std::cout<<" -ringBC_type: "        <<ringBC_type        <<std::endl;
  std::cout<<" -num_inlet: "          <<num_inlet          <<std::endl;
  std::cout<<" -num_outlet: "         <<num_outlet         <<std::endl;
  std::cout<<" -geo_file: "           <<geo_file           <<std::endl;
  std::cout<<" -geo_f_file: "         <<geo_f_file         <<std::endl;
  std::cout<<" -geo_s_file: "         <<geo_s_file         <<std::endl;
  std::cout<<" -sur_f_file_wall: "    <<sur_f_file_wall    <<std::endl;
  std::cout<<" -sur_s_file_wall: "    <<sur_s_file_wall    <<std::endl;
  std::cout<<" -sur_f_file_in_base: " <<sur_f_file_in_base <<std::endl;
  std::cout<<" -sur_f_file_out_base: "<<sur_f_file_out_base<<std::endl;
  std::cout<<" -sur_s_file_in_base: " <<sur_s_file_in_base <<std::endl;
  std::cout<<" -sur_s_file_out_base: "<<sur_s_file_out_base<<std::endl;
  std::cout<<" -part_file: "          <<part_file          <<std::endl;
  std::cout<<" -cpu_size: "           <<cpu_size           <<std::endl;
  std::cout<<" -in_ncommon: "         <<in_ncommon         <<std::endl;
  std::cout<<" -isDualGraph: true \n";
  if(isReload) std::cout<<" -isReload : true \n";
  else std::cout<<" -isReload : false \n";
  std::cout<<"----------------------------------\n";
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<"===== Command Line Arguments ====="<<std::endl;

  // Check if the geometrical file exist on disk
  SYS_T::file_check(geo_file); std::cout<<geo_file<<" found. \n";

  SYS_T::file_check(geo_f_file); std::cout<<geo_f_file<<" found. \n";

  SYS_T::file_check(geo_s_file); std::cout<<geo_s_file<<" found. \n";

  SYS_T::file_check(sur_f_file_wall); std::cout<<sur_f_file_wall<<" found. \n";

  SYS_T::file_check(sur_s_file_wall); std::cout<<sur_s_file_wall<<" found. \n";

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

  // If we can still detect additional files on disk, throw an warning
  if( SYS_T::file_exist(SYS_T::gen_capfile_name(sur_f_file_in_base, num_inlet, ".vtp")) ||
      SYS_T::file_exist(SYS_T::gen_capfile_name(sur_s_file_in_base, num_inlet, ".vtp")) )
    cout<<endl<<"Warning: there are additional inlet surface files on disk. Check num_inlet please.\n\n";

  if( SYS_T::file_exist(SYS_T::gen_capfile_name(sur_f_file_out_base, num_outlet, ".vtp")) ||
      SYS_T::file_exist(SYS_T::gen_capfile_name(sur_s_file_out_base, num_outlet, ".vtp")) )
    cout<<endl<<"Warning: there are additional outlet surface files on disk. Check num_outlet please.\n\n";

  // ----- Write the input argument into a HDF5 file
  SYS_T::execute("rm -rf preprocessor_cmd.h5");
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("num_outlet",       num_outlet);
  cmdh5w->write_intScalar("num_inlet",        num_inlet);
  cmdh5w->write_intScalar("cpu_size",         cpu_size);
  cmdh5w->write_intScalar("in_ncommon",       in_ncommon);
  cmdh5w->write_intScalar("dofNum",           dofNum);
  cmdh5w->write_intScalar("dofMat",           dofMat);
  cmdh5w->write_intScalar("elemType",         elemType);
  cmdh5w->write_intScalar("fsiBC_type",       fsiBC_type);
  cmdh5w->write_intScalar("ringBC_type",      ringBC_type);
  cmdh5w->write_string("geo_file",            geo_file);
  cmdh5w->write_string("geo_f_file",          geo_f_file);
  cmdh5w->write_string("geo_s_file",          geo_s_file);
  cmdh5w->write_string("sur_f_file_in_base",  sur_f_file_in_base);
  cmdh5w->write_string("sur_f_file_out_base", sur_f_file_out_base);
  cmdh5w->write_string("sur_f_file_wall",     sur_f_file_wall);
  cmdh5w->write_string("sur_s_file_in_base",  sur_s_file_in_base);
  cmdh5w->write_string("sur_s_file_out_base", sur_s_file_out_base);
  cmdh5w->write_string("sur_s_file_wall",     sur_s_file_wall);
  cmdh5w->write_string("part_file",           part_file);
  cmdh5w->write_string("date",                SYS_T::get_date() );
  cmdh5w->write_string("time",                SYS_T::get_time() );

  delete cmdh5w; H5Fclose(cmd_file_id);
  // ----- Finish writing

  // Read the geometry file for the whole FSI domain
  int nFunc_v, nElem;
  std::vector<int> vecIEN, phy_tag;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid( geo_file, nFunc_v, nElem, ctrlPts, vecIEN, phy_tag );

  for(unsigned int ii=0; ii<phy_tag.size(); ++ii)
  {
    if(phy_tag[ii] != 0 && phy_tag[ii] != 1) SYS_T::print_fatal("Error: FSI problem, the physical tag for element should be 0 (fluid domain) or 1 (solid domain).\n");
  }

  // Generate IEN
  IIEN * IEN_v = new IEN_Tetra_P1( nElem, vecIEN );

  // --------------------------------------------------------------------------
  // The fluid-solid interface file will be read and the nodal index will be
  // mapped to a new value by the following rule. The ii-th node in the
  // interface wall node will be assgiend nFunc + ii. 
  // Read the F-S interface vtp file
  const std::vector<int> wall_node_id = TET_T::read_int_PointData( sur_s_file_interior_wall, "GlobalNodeID" );

  VEC_T::print(wall_node_id);
  
  const int nFunc_interface = static_cast<int>( wall_node_id.size() );
  const int nFunc_p = nFunc_v + nFunc_interface;

  // We will generate a new IEN array for the pressure variable by updating the
  // IEN for the solid element. If the solid element has node on the fluid-solid
  // interface, it will be mapped to the new index, that is nFunc + ii.
  std::vector<int> vecIEN_p ( vecIEN );
  
  for(int ee=0; ee<nElem; ++ee)
  {
    if( phy_tag[ee] == 1 )
    {
      //std::cout<<"ee = "<<ee<<'\t'<<vecIEN_p[ee*4]<<'\t'<<vecIEN_p[ee*4+1]<<'\t'<<vecIEN_p[ee*4+2]<<'\t'<<vecIEN_p[ee*4+3]<<'\n';

      // In solid element, loop over its IEN and correct if the node is on the
      // interface
      for(int ii=0; ii<4; ++ii)
      {
        const int pos = VEC_T::get_pos( wall_node_id, vecIEN_p[ee*4 +ii] );
        if( pos >=0 ) vecIEN_p[ee*4+ii] = nFunc_v + pos;     
      }
      
      //std::cout<<"ee = "<<ee<<'\t'<<vecIEN_p[ee*4]<<'\t'<<vecIEN_p[ee*4+1]<<'\t'<<vecIEN_p[ee*4+2]<<'\t'<<vecIEN_p[ee*4+3]<<'\n';
    }
  }
  
  IIEN * IEN_p = new IEN_Tetra_P1( nElem, vecIEN_p );

  VEC_T::clean( vecIEN ); VEC_T::clean( vecIEN_p );
  // --------------------------------------------------------------------------

  // Generate the list of nodes for fluid and solid
  std::vector<int> node_f, node_s; node_f.clear(); node_s.clear();

  for(int ee=0; ee<nElem; ++ee)
  {
    if( phy_tag[ee] == 0 )
    {
      for(int ii=0; ii<4; ++ii) node_f.push_back( IEN_v->get_IEN(ee, ii) );
    }
    else
    {
      for(int ii=0; ii<4; ++ii) node_s.push_back( IEN_v->get_IEN(ee, ii) );
    }
  }

  VEC_T::sort_unique_resize( node_f );
  VEC_T::sort_unique_resize( node_s );

  // Check the mesh
  const double critical_val_aspect_ratio = 3.5;
  TET_T::tetmesh_check( ctrlPts, IEN_v, nElem, critical_val_aspect_ratio );

  // Generate the mesh
  IMesh * mesh_v = new Mesh_Tet4(nFunc_v, nElem);

  // Generate the mesh for pressure
  IMesh * mesh_p = new Mesh_Tet4(nFunc_p, nElem);
  
  std::vector<IMesh const *> mlist;
  mlist.push_back(mesh_p); mlist.push_back(mesh_v);

  mlist[0]->print_info();
  mlist[1]->print_info();

  std::cout<<"Fluid domain: "<<node_f.size()<<" nodes.\n";
  std::cout<<"Solid domain: "<<node_s.size()<<" nodes.\n";
  std::cout<<"Fluid-Solid interface: "<<nFunc_interface<<" nodes.\n";
  
  VEC_T::clean( node_f ); VEC_T::clean( node_s );

  std::vector<IIEN const *> ienlist;
  ienlist.push_back(IEN_p); ienlist.push_back(IEN_v);

  // Partition the mesh
  IGlobal_Part * global_part = nullptr;
  if( isReload ) global_part = new Global_Part_Reload( cpu_size, in_ncommon, isDualGraph );
  else
  {
    if(cpu_size > 1)
    {
      global_part = new Global_Part_METIS( num_fields, cpu_size, in_ncommon, isDualGraph, 
          mlist, ienlist );
    }
    else if(cpu_size == 1)
      global_part = new Global_Part_Serial( num_fields, mlist );
    else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);
  }

  // Re-ordering nodal indices
  const int nFunc_total = nFunc_p + nFunc_v;

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, nFunc_total);

  mnindex->write_hdf5("node_mapping");

  Map_Node_Index * mnindex_p = new Map_Node_Index(global_part, cpu_size, 0, nFunc_p);
  Map_Node_Index * mnindex_v = new Map_Node_Index(global_part, cpu_size, nFunc_p, nFunc_total);

  mnindex_p -> write_hdf5("node_mapping_p");
  mnindex_v -> write_hdf5("node_mapping_v");

  // Generate a list of local node number
  std::vector<int> list_nn_v(cpu_size), list_nn_p(cpu_size);
  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    int num_node_v = 0, num_node_p = 0;
    for(int nn=0; nn<mesh_p -> get_nFunc(); ++nn)
    {
      if(global_part->get_npart(nn) == proc_rank) num_node_p += 1;
    }

    for(int nn=0; nn<mesh_v -> get_nFunc(); ++nn)
    {
      if(global_part->get_npart(nn+nFunc_p) == proc_rank) num_node_v += 1;
    }

    // list stores the number of velo/pres nodes in each cpu
    list_nn_v[proc_rank] = num_node_v;
    list_nn_p[proc_rank] = num_node_p;
  }

  // Now generate the mappings from the gird pt idx to the matrix row idx
  // This is needed because will have a matrix that has a special structure
  // due to the use of mix fem.
  std::vector<int> start_idx_v(cpu_size), start_idx_p(cpu_size);
  start_idx_v[0] = 0;
  start_idx_p[0] = 3 * list_nn_v[0];
  for(int ii = 1; ii < cpu_size; ++ii )
  {
    start_idx_v[ii] = start_idx_v[ii-1] + list_nn_v[ii-1]*3 + list_nn_p[ii-1];
    start_idx_p[ii] = start_idx_v[ii  ] + list_nn_v[ii  ]*3;
  }

  // mapper maps from the new grid point index to the matrix row index
  std::vector< std::vector<int> > mapper_v;
  std::vector<int> mapper_p;
  mapper_v.resize(3);

  for(int ii=0; ii<cpu_size; ++ii)
  {
    for(int jj=0; jj<list_nn_v[ii]; ++jj)
    {
      mapper_v[0].push_back( start_idx_v[ii] + jj * 3 );
      mapper_v[1].push_back( start_idx_v[ii] + jj * 3 + 1 );
      mapper_v[2].push_back( start_idx_v[ii] + jj * 3 + 2 );
    }

    for(int jj=0; jj<list_nn_p[ii]; ++jj)
      mapper_p.push_back( start_idx_p[ii] + jj );
  }

  // Partition the mesh












  PetscFinalize();
  cout<<"===> Preprocessing completes successfully!\n";
  return EXIT_SUCCESS;
}

// EOF
