// ============================================================================
// preprocess_fsi.cpp
// ----------------------------------------------------------------------------
// This is the preprocess code for handling FSI 3D meshes that uses tet for the
// geometry and discontinuous pressure space.
//
// Author: Ju Liu
// Date: Dec. 13 2021
// ============================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet.hpp"
#include "Mesh_FEM.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Global_Part_Reload.hpp"
#include "Part_FEM_FSI_Tissue.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_FSI.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "NBC_Partition_MF.hpp"
#include "NBC_Partition_inflow_MF.hpp"
#include "EBC_Partition_outflow_MF.hpp"
#include "yaml-cpp/yaml.h"

int main( int argc, char * argv[] )
{
  // Set number of threads and  print info of OpenMP
  SYS_T::print_omp_info();
  SYS_T::set_omp_num_threads();
  
  // Remove previously existing hdf5 files
  if( SYS_T::directory_exist("apart") )
  {
    std::cout<<"Clean the folder apart.\n";
    SYS_T::execute("rm -rf apart");
  }
  SYS_T::execute("mkdir apart");

  // Define basic settings
  constexpr int num_fields = 2; // Two fields : pressure + velocity/displacement
  const std::vector<int> dof_fields {1, 3}; // pressure 1 ; velocity/displacement 3

  // Yaml options
  const std::string yaml_file("fsi_preprocess.yml");

  // Check if the yaml file exist on disk
  SYS_T::file_check(yaml_file);

  YAML::Node paras = YAML::LoadFile( yaml_file );

  const int elemType                          = paras["elem_type"].as<int>();
  const int fsiBC_type                        = paras["fsiBC_type"].as<int>();
  const int ringBC_type                       = paras["ringBC_type"].as<int>();
  const int num_inlet                         = paras["num_inlet"].as<int>();
  const int num_outlet                        = paras["num_outlet"].as<int>();
  const std::string geo_file                  = paras["geo_file"].as<std::string>();
  const std::string geo_f_file                = paras["geo_f_file"].as<std::string>();
  const std::string geo_s_file                = paras["geo_s_file"].as<std::string>();

  const std::string sur_f_file_wall           = paras["sur_f_file_wall"].as<std::string>();
  const std::string sur_f_file_in_base        = paras["sur_f_file_in_base"].as<std::string>();
  const std::string sur_f_file_out_base       = paras["sur_f_file_out_base"].as<std::string>();

  const std::string sur_s_file_interior_wall  = paras["sur_s_file_interior_wall"].as<std::string>();
  const std::string sur_s_file_wall           = paras["sur_s_file_wall"].as<std::string>();
  const std::string sur_s_file_in_base        = paras["sur_s_file_in_base"].as<std::string>();
  const std::string sur_s_file_out_base       = paras["sur_s_file_out_base"].as<std::string>();

  const std::string part_file_p               = paras["part_file_p"].as<std::string>();
  const std::string part_file_v               = paras["part_file_v"].as<std::string>();
  const int cpu_size                          = paras["cpu_size"].as<int>();
  const int in_ncommon                        = paras["in_ncommon"].as<int>();
  const bool isDualGraph                      = paras["is_dualgraph"].as<bool>();
  const bool isReload                         = paras["is_reload"].as<bool>();

  const bool isPrintMeshQual                  = paras["is_printmeshqual"].as<bool>();
  const double critical_val_aspect_ratio      = paras["critical_val_aspect_ratio"].as<double>();

  SYS_T::print_fatal_if( elemType != 501 && elemType != 601 , "ERROR: unknown element type %d.\n", elemType );
  SYS_T::print_fatal_if( fsiBC_type != 0 && fsiBC_type != 1 && fsiBC_type != 2, "Error: fsiBC_type should be 0, 1, or 2.\n" );
  SYS_T::print_fatal_if( ringBC_type != 0, "Error: ringBC_type should be 0.\n" );

  std::cout<<"===== Command Line Arguments ====="<<std::endl;
  std::cout<<" -elem_type: "          <<elemType           <<std::endl;
  std::cout<<" -fsiBC_type: "         <<fsiBC_type         <<std::endl;
  std::cout<<" -ringBC_type: "        <<ringBC_type        <<std::endl;
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
  std::cout<<" -cpu_size: "           <<cpu_size           <<std::endl;
  std::cout<<" -in_ncommon: "         <<in_ncommon         <<std::endl;
  if(isDualGraph) std::cout<<" -isDualGraph: true \n";
  else std::cout<<" -isDualGraph: false \n";
  if(isReload) std::cout<<" -isReload : true \n";
  else std::cout<<" -isReload : false \n";
  if(isPrintMeshQual) std::cout<<" -isPrintMeshQual: true \n";
  else std::cout<<" -isPrintMeshQual: false \n";
  std::cout<<" -critical_val_aspect_ratio: "<<critical_val_aspect_ratio<<std::endl;
  std::cout<<"----------------------------------\n";
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<" part_file_p: "<<part_file_p<<std::endl;
  std::cout<<" part_file_v: "<<part_file_v<<std::endl;
  std::cout<<"===== Command Line Arguments ====="<<std::endl;

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
  cmdh5w->write_string("sur_s_file_interior_wall", sur_s_file_interior_wall);
  cmdh5w->write_string("part_file_p",         part_file_p);
  cmdh5w->write_string("part_file_v",         part_file_v);
  cmdh5w->write_string("date",                SYS_T::get_date() );
  cmdh5w->write_string("time",                SYS_T::get_time() );

  delete cmdh5w; H5Fclose(cmd_file_id);
  // ----- Finish writing

  // Read the geometry file for the whole FSI domain for the velocity /
  // displacement field
  int nFunc_v, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid( geo_file, nFunc_v, nElem, ctrlPts, vecIEN );

  const std::vector<int> phy_tag = VTK_T::read_int_CellData(geo_file, "Physics_tag");

  for(unsigned int ii=0; ii<phy_tag.size(); ++ii)
  {
    if(phy_tag[ii] != 0 && phy_tag[ii] != 1) SYS_T::print_fatal("Error: FSI problem, the physical tag for element should be 0 (fluid domain) or 1 (solid domain).\n");
  }

  // Generate IEN
  IIEN * IEN_v = new IEN_FEM( nElem, vecIEN );

  // --------------------------------------------------------------------------
  // The fluid-solid interface file will be read and the nodal index will be
  // mapped to a new value by the following rule. The ii-th node in the
  // interface wall node will be assgiend to nFunc_v + ii. 
  // Read the F-S interface vtp file
  const std::vector<int> wall_node_id = VTK_T::read_int_PointData( sur_s_file_interior_wall, "GlobalNodeID" );

  const int nFunc_interface = static_cast<int>( wall_node_id.size() );
  const int nFunc_p = nFunc_v + nFunc_interface;

  // We will generate a new IEN array for the pressure variable by updating the
  // IEN for the solid element. If the solid element has node on the fluid-solid
  // interface, it will be mapped to the new index, that is nFunc + ii.
  std::vector<int> vecIEN_p ( vecIEN );
  PERIGEE_OMP_PARALLEL_FOR
  for(int ee=0; ee<nElem; ++ee)
  {
    if(phy_tag[ee] == 1)
    {
      // In solid element, loop over its IEN and correct if the node is on the
      // interface
      if(elemType == 501)
      {
        for(int ii=0; ii<4; ++ii)
        {
          const int pos = VEC_T::get_pos( wall_node_id, vecIEN_p[ee*4+ii] );
          if( pos >=0 ) vecIEN_p[ee*4+ii] = nFunc_v + pos;     
        }
      }
      else if(elemType == 601)
      {
        for(int ii=0; ii<8; ++ii)
        {
          const int pos = VEC_T::get_pos( wall_node_id, vecIEN_p[ee*8+ii] );
          if( pos >=0 ) vecIEN_p[ee*8+ii] = nFunc_v + pos;     
        }
      }
      else
        SYS_T::print_fatal("Error: elemType %d is not supported when generating a new IEN array for the pressure variable. \n", elemType);
    }
  }

  IIEN * IEN_p = new IEN_FEM( nElem, vecIEN_p );

  VEC_T::clean( vecIEN ); VEC_T::clean( vecIEN_p );
  // --------------------------------------------------------------------------

  // Generate the list of nodes for fluid and solid
  std::vector<int> v_node_f, v_node_s; v_node_f.clear(); v_node_s.clear();
  PERIGEE_OMP_PARALLEL
  {
    std::vector<int> temp_v_node_f {};
    std::vector<int> temp_v_node_s {};
    PERIGEE_OMP_FOR
    for(int ee=0; ee<nElem; ++ee)
    {
      if( phy_tag[ee] == 0 )
      {
        if(elemType == 501)
        {
          for(int ii=0; ii<4; ++ii) temp_v_node_f.push_back( IEN_v->get_IEN(ee, ii) );
        }
        else if(elemType == 601)
        {
          for(int ii=0; ii<8; ++ii) temp_v_node_f.push_back( IEN_v->get_IEN(ee, ii) );
        }
        else
          SYS_T::print_fatal("Error: elemType %d is not supported when generating the list of velocity nodes for fluid during the preprocessing. \n", elemType);
      }
      else
      {
        if(elemType == 501)
        {
          for(int ii=0; ii<4; ++ii) temp_v_node_s.push_back( IEN_v->get_IEN(ee, ii) );
        }
        else if(elemType == 601)
        {
          for(int ii=0; ii<8; ++ii) temp_v_node_s.push_back( IEN_v->get_IEN(ee, ii) );
        }
        else
          SYS_T::print_fatal("Error: elemType %d is not supported when generating the list of velocity nodes for solid during the preprocessing. \n", elemType);
      }
    }
    PERIGEE_OMP_CRITICAL
    {
      VEC_T::insert_end(v_node_f, temp_v_node_f);
      VEC_T::insert_end(v_node_s, temp_v_node_s);
    }
  }

  VEC_T::sort_unique_resize( v_node_f ); VEC_T::sort_unique_resize( v_node_s );

  std::vector<int> p_node_f, p_node_s; p_node_f.clear(); p_node_s.clear();
  PERIGEE_OMP_PARALLEL
  {
    std::vector<int> temp_p_node_f {};
    std::vector<int> temp_p_node_s {};
    PERIGEE_OMP_FOR
    for(int ee=0; ee<nElem; ++ee)
    {
      if( phy_tag[ee] == 0 )
      {
        if(elemType == 501)
        {
          for(int ii=0; ii<4; ++ii) temp_p_node_f.push_back( IEN_p->get_IEN(ee, ii) );
        }
        else if(elemType == 601)
        {
          for(int ii=0; ii<8; ++ii) temp_p_node_f.push_back( IEN_p->get_IEN(ee, ii) );
        }
        else
          SYS_T::print_fatal("Error: elemType %d is not supported when generating the list of pressure nodes for fluid during the preprocessing. \n", elemType);
      }
      else
      {
        if(elemType == 501)
        {
          for(int ii=0; ii<4; ++ii) temp_p_node_s.push_back( IEN_p->get_IEN(ee, ii) );
        }
        else if(elemType == 601)
        {
          for(int ii=0; ii<8; ++ii) temp_p_node_s.push_back( IEN_p->get_IEN(ee, ii) );
        }
        else
          SYS_T::print_fatal("Error: elemType %d is not supported when generating the list of pressure nodes for solid during the preprocessing. \n", elemType);
      }
    }
    PERIGEE_OMP_CRITICAL
    {
      VEC_T::insert_end(p_node_f, temp_p_node_f);
      VEC_T::insert_end(p_node_s, temp_p_node_s);
    }
  }

  VEC_T::sort_unique_resize( p_node_f ); VEC_T::sort_unique_resize( p_node_s );

  // Check the mesh of kinematics
  if( isPrintMeshQual )
  {
    std::cout<<"Check the mesh quality... \n";
    if(elemType == 501)
    {
      TET_T::tetmesh_check( ctrlPts, IEN_v, nElem, critical_val_aspect_ratio );
    }
    else if(elemType == 601)
    {
      HEX_T::hexmesh_check( ctrlPts, IEN_v, nElem, critical_val_aspect_ratio );
    }
    else
      SYS_T::print_fatal("Error: elemType %d is not supported when checking the mesh of kinematics. \n", elemType);
  }

  // Generate the mesh for kinematics
  IMesh * mesh_v = nullptr;

  // Generate the mesh for pressure (discontinuous over interface)
  IMesh * mesh_p = nullptr;

  switch( elemType )
  {
    case 501:
      mesh_v = new Mesh_Tet(nFunc_v, nElem, 1);
      mesh_p = new Mesh_Tet(nFunc_p, nElem, 1);
      break;
    case 601:
      mesh_v = new Mesh_FEM(nFunc_v, nElem, 8, 1);
      mesh_p = new Mesh_FEM(nFunc_p, nElem, 8, 1);
      break;   
    default:
      SYS_T::print_fatal("Error: elemType %d is not supported when generating the mesh.\n", elemType);
      break;
  }

  SYS_T::print_fatal_if( IEN_v->get_nLocBas() != mesh_v->get_nLocBas(), "Error: the nLocBas from the mesh_v %d and the IEN_v %d classes do not match. \n", mesh_v->get_nLocBas(), IEN_v->get_nLocBas() );
  SYS_T::print_fatal_if( IEN_p->get_nLocBas() != mesh_p->get_nLocBas(), "Error: the nLocBas from the mesh_p %d and the IEN_p %d classes do not match. \n", mesh_p->get_nLocBas(), IEN_p->get_nLocBas() );

  std::vector<IMesh const *> mlist;
  mlist.push_back(mesh_p); mlist.push_back(mesh_v);

  mlist[0]->print_info();
  mlist[1]->print_info();

  std::cout<<"Fluid domain: "<<v_node_f.size()<<" nodes.\n";
  std::cout<<"Solid domain: "<<v_node_s.size()<<" nodes.\n";
  std::cout<<"Fluid-Solid interface: "<<nFunc_interface<<" nodes.\n";

  // --------------------------------------------------------------------------
  // Read the geometry file for the solid domain, generate the list of direction
  // basis vectors of the nodes. The list includes radial, longitudinal, and
  // circumferential basis, denoting by r, l, and c, respectively.
  const std::vector<int> solid_node_id = VTK_T::read_int_PointData(geo_s_file, "GlobalNodeID");
  const std::vector<Vector_3> basis_r  = VTK_T::read_Vector_3_PointData(geo_s_file, "radial_basis");
  const std::vector<Vector_3> basis_l  = VTK_T::read_Vector_3_PointData(geo_s_file, "longitudinal_basis");
  const std::vector<Vector_3> basis_c  = VTK_T::read_Vector_3_PointData(geo_s_file, "circumferential_basis");

  SYS_T::print_fatal_if(v_node_s != solid_node_id, "ERROR: GlobalNodeID for solid geometry file is not equal to the whole FSI domain.");
  SYS_T::print_fatal_if(solid_node_id.size() != basis_r.size(), "ERROR: radial_basis is not matched.");
  SYS_T::print_fatal_if(solid_node_id.size() != basis_l.size(), "ERROR: longitudinal_basis is not matched.");
  SYS_T::print_fatal_if(solid_node_id.size() != basis_c.size(), "ERROR: circumferential_basis is not matched.");
  std::cout<<"=== Direction basis vectors generated.\n";
  // --------------------------------------------------------------------------

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
  Map_Node_Index * mnindex_p = new Map_Node_Index(global_part, cpu_size, nFunc_p, 0);
  Map_Node_Index * mnindex_v = new Map_Node_Index(global_part, cpu_size, nFunc_v, 1);

  mnindex_p -> write_hdf5("node_mapping_p");
  mnindex_v -> write_hdf5("node_mapping_v");

  // Generate a list of local node number
  std::vector<int> list_nn_v(cpu_size), list_nn_p(cpu_size);
  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    // list stores the number of velo/pres nodes in each cpu
    list_nn_p[proc_rank] = 0; list_nn_v[proc_rank] = 0;
    for(int nn=0; nn<mesh_p -> get_nFunc(); ++nn)
    {
      if(global_part->get_npart(nn,0) == proc_rank) list_nn_p[proc_rank] += 1;
    }

    for(int nn=0; nn<mesh_v -> get_nFunc(); ++nn)
    {
      if(global_part->get_npart(nn,1) == proc_rank) list_nn_v[proc_rank] += 1;
    }
  }

  // Now generate the mappings from the gird pt idx to the matrix row idx
  // This is needed because we will have a matrix that has a special structure
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
  std::vector< std::vector<int> > mapper_p, mapper_v;
  mapper_p.resize(1); mapper_v.resize(3);

  for(int ii=0; ii<cpu_size; ++ii)
  {
    for(int jj=0; jj<list_nn_v[ii]; ++jj)
    {
      mapper_v[0].push_back( start_idx_v[ii] + jj * 3     );
      mapper_v[1].push_back( start_idx_v[ii] + jj * 3 + 1 );
      mapper_v[2].push_back( start_idx_v[ii] + jj * 3 + 2 );
    }

    for(int jj=0; jj<list_nn_p[ii]; ++jj)
      mapper_p[0].push_back( start_idx_p[ii] + jj );
  }

  // ----------------------------------------------------------------
  // Setup boundary conditions. Nodal BC is specified by the original nodal
  // indices
  // Physical NodalBC
  std::cout<<"===== Boundary Conditions =====\n";
  std::cout<<"1. Nodal boundary condition for the implicit solver: \n";
  std::vector<INodalBC *> NBC_list_p( 1, nullptr );
  std::vector<INodalBC *> NBC_list_v( 3, nullptr );

  // Here we assumed that the pressure mesh fluid nodal indices are identical to
  // that in the velocity mesh.
  NBC_list_p[0] = new NodalBC_3D_FSI( geo_f_file, nFunc_p, fsiBC_type );

  for( int ii=0; ii<3; ++ii )
    NBC_list_v[ii] = new NodalBC_3D_FSI( geo_f_file, geo_s_file, sur_f_file_wall, 
        sur_s_file_wall, sur_f_file_in, sur_f_file_out, sur_s_file_in, sur_s_file_out, 
        nFunc_v, ii, ringBC_type, fsiBC_type );

  // Mesh solver NodalBC
  std::cout<<"2. Nodal boundary condition for the mesh motion: \n";
  std::vector<INodalBC *> meshBC_list( 3, nullptr );

  std::vector<std::string> meshdir_file_list { geo_s_file };
  VEC_T::insert_end( meshdir_file_list, sur_f_file_in );
  VEC_T::insert_end( meshdir_file_list, sur_f_file_out );

  meshBC_list[0] = new NodalBC( meshdir_file_list, nFunc_v );
  meshBC_list[1] = new NodalBC( meshdir_file_list, nFunc_v );
  meshBC_list[2] = new NodalBC( meshdir_file_list, nFunc_v );

  // InflowBC info
  std::cout<<"3. Inflow cap surfaces: \n";
  std::vector<Vector_3> inlet_outvec( num_inlet );
  if(elemType == 501)
  {
    for(int ii=0; ii<num_inlet; ++ii)
      inlet_outvec[ii] = TET_T::get_out_normal( sur_f_file_in[ii], ctrlPts, IEN_v );
  }
  else if(elemType == 601)
  {
    for(int ii=0; ii<num_inlet; ++ii)
      inlet_outvec[ii] = HEX_T::get_out_normal( sur_f_file_in[ii], ctrlPts, IEN_v );
  }
  else
    SYS_T::print_fatal("Error: elemType %d is not supported when obtaining the outward normal vector for the inflow boundary condition. \n", elemType);

  INodalBC * InFBC = new NodalBC_3D_inflow( sur_f_file_in, sur_f_file_wall, nFunc_v, inlet_outvec, elemType );

  InFBC -> resetSurIEN_outwardnormal( IEN_v ); // assign outward orientation for triangles
  
  // Physical ElemBC
  cout<<"4. Elem boundary for the implicit solver: \n";
  std::vector< Vector_3 > outlet_outvec( num_outlet );

  if(elemType == 501)
  {
    for(int ii=0; ii<num_outlet; ++ii)
      outlet_outvec[ii] = TET_T::get_out_normal( sur_f_file_out[ii], ctrlPts, IEN_v );
  }
  else if(elemType == 601)
  {
    for(int ii=0; ii<num_outlet; ++ii)
      outlet_outvec[ii] = HEX_T::get_out_normal( sur_f_file_out[ii], ctrlPts, IEN_v );    
  }
  else
    SYS_T::print_fatal("Error: elemType %d is not supported when obtaining the outward normal vector for the elemental boundary conditions. \n", elemType);

  std::vector< std::string > ebclist {sur_f_file_wall};

  ElemBC * ebc = nullptr;
  if( fsiBC_type == 0 || fsiBC_type == 1 )
    ebc = new ElemBC_3D_outflow( sur_f_file_out, outlet_outvec, elemType );
  else if( fsiBC_type == 2 )
    ebc = new ElemBC_3D( ebclist, elemType );
  else SYS_T::print_fatal("ERROR: uncognized fsiBC type. \n");

  ebc -> resetSurIEN_outwardnormal( IEN_v ); // assign outward orientation for triangles

  // Mesh solver ElemBC
  cout<<"5. Elem boundary for the mesh solver: \n";
  std::vector<std::string> mesh_ebclist;
  mesh_ebclist.clear();
  ElemBC * mesh_ebc = new ElemBC_3D( mesh_ebclist, elemType );
  std::cout<<"=================================\n";
  // ----------------------------------------------------------------

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();

    IPart * part_p = new Part_FEM_FSI( mesh_p, global_part, mnindex_p, IEN_p,
        ctrlPts, phy_tag, p_node_f, p_node_s, proc_rank, cpu_size, elemType, 
        start_idx_p[proc_rank], {0, dof_fields[0], false, "pressure"} );

    part_p -> print_part_loadbalance_edgecut();
    
    part_p -> write( part_file_p );

    IPart * part_v = new Part_FEM_FSI_Tissue( mesh_v, global_part, mnindex_v, IEN_v,
        ctrlPts, phy_tag, v_node_f, v_node_s, basis_r, basis_l, basis_c,
	proc_rank, cpu_size, elemType, start_idx_v[proc_rank], { 1, dof_fields[1], true, "velocity"} );

    part_v -> print_part_loadbalance_edgecut();

    part_v -> write( part_file_v );

    mytimer -> Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    NBC_Partition * nbcpart_p = new NBC_Partition_MF(part_p, mnindex_p, NBC_list_p, mapper_p);
    nbcpart_p -> write_hdf5( part_file_p );

    NBC_Partition * nbcpart_v = new NBC_Partition_MF(part_v, mnindex_v, NBC_list_v, mapper_v);
    nbcpart_v -> write_hdf5( part_file_v );

    NBC_Partition * mbcpart = new NBC_Partition_MF(part_v, mnindex_v, meshBC_list);
    mbcpart -> write_hdf5( part_file_v, "/mesh_nbc" );

    NBC_Partition_inflow * infpart = new NBC_Partition_inflow_MF(part_v, mnindex_v, InFBC, mapper_v);
    infpart->write_hdf5( part_file_v );

    if( fsiBC_type == 0 || fsiBC_type == 1 )
    {
      EBC_Partition_outflow_MF * ebcpart = new EBC_Partition_outflow_MF(part_v, mnindex_v, ebc, NBC_list_v, mapper_v);
      ebcpart -> write_hdf5( part_file_v );
      delete ebcpart;
    }
    else if( fsiBC_type == 2 )
    {
      EBC_Partition * ebcpart = new EBC_Partition( part_v, mnindex_v, ebc );
      ebcpart -> write_hdf5( part_file_v );
      delete ebcpart;
    }
    else SYS_T::print_fatal("ERROR: unrecognized fsiBC_type. \n");

    EBC_Partition * ebcpart_p = new EBC_Partition( part_p, mnindex_p, ebc );
    ebcpart_p -> write_hdf5( part_file_p );

    EBC_Partition * mebcpart = new EBC_Partition(part_v, mnindex_v, mesh_ebc);

    mebcpart-> write_hdf5( part_file_v, "/mesh_ebc" );

    delete part_p; delete part_v; delete nbcpart_p; delete nbcpart_v;
    delete mbcpart; delete infpart; delete mebcpart; 
  }

  // Clean up the memory
  for(auto &it_nbc : NBC_list_v) delete it_nbc;
  
  for(auto &it_nbc : NBC_list_p) delete it_nbc;

  for(auto &it_nbc : meshBC_list) delete it_nbc;

  delete ebc; delete InFBC; delete mesh_ebc; 
  delete mnindex_p; delete mnindex_v; delete mesh_p; delete mesh_v; 
  delete IEN_p; delete IEN_v; delete mytimer; delete global_part; 

  cout<<"===> Preprocessing completes successfully!\n";
  return EXIT_SUCCESS;
}

// EOF
