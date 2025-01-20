// ==================================================================
// preprocess_ns.cpp
//
// This is a preprocessor code for handling Navier-Stokes equations 
// discretized by tetradedral elements.
//
// Date Created: Jan 01 2020
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet.hpp"
#include "Mesh_FEM.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM_Rotated.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "NodalBC_3D_rotated.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "ElemBC_3D_turbulence_wall_model.hpp"
#include "Interface_pair.hpp"
#include "NBC_Partition.hpp"
#include "NBC_Partition_inflow.hpp"
#include "NBC_Partition_rotated.hpp"
#include "EBC_Partition_outflow.hpp"
#include "EBC_Partition_turbulence_wall_model.hpp"
#include "Interface_Partition.hpp"
#include "yaml-cpp/yaml.h"

int main( int argc, char * argv[] )
{
  // Set number of threads and  print info of OpenMP
  SYS_T::print_omp_info();
  SYS_T::set_omp_num_threads();

  // Clean the potentially pre-existing hdf5 files in the job folder
  SYS_T::execute("rm -rf part_p*.h5");
  SYS_T::execute("rm -rf preprocessor_cmd.h5");
  SYS_T::execute("rm -rf *_itf.h5");

  // Define basic problem settins
  constexpr int dofNum = 4; // degree-of-freedom for the physical problem
  constexpr int dofMat = 4; // degree-of-freedom in the matrix problem

  // Yaml options
  const std::string yaml_file("ns_preprocess.yml");

  // Check if the yaml file exist on disk
  SYS_T::file_check(yaml_file);

  YAML::Node paras = YAML::LoadFile( yaml_file );

  const std::string elemType_str      = paras["elem_type"].as<std::string>();
  const int num_inlet                 = paras["num_inlet"].as<int>();
  const int num_outlet                = paras["num_outlet"].as<int>();
  const std::string fixed_geo_file          = paras["fixed_geo_file"].as<std::string>();
  const std::string sur_file_in_base  = paras["sur_file_in_base"].as<std::string>();
  const std::string sur_file_inner_wall     = paras["sur_file_inner_wall"].as<std::string>();
  const std::string sur_file_outer_wall     = paras["sur_file_outer_wall"].as<std::string>();

  const std::string sur_file_out_base = paras["sur_file_out_base"].as<std::string>();

  const int num_interface_pair        = paras["num_interface_pair"].as<int>();
  const std::string rotated_geo_file  = paras["rotated_geo_file"].as<std::string>();
  const std::string rotated_sur_file  = paras["rotated_sur_file"].as<std::string>();
  const std::string fixed_interface_base   = paras["fixed_interface_base"].as<std::string>();
  const std::string rotated_interface_base = paras["rotated_interface_base"].as<std::string>();

  const std::string part_file         = paras["part_file"].as<std::string>();
  const int cpu_size                  = paras["cpu_size"].as<int>();
  const int in_ncommon                = paras["in_ncommon"].as<int>();
  const bool isDualGraph              = paras["is_dualgraph"].as<bool>();
  const FEType elemType               = FE_T::to_FEType(elemType_str);

  // Optional:
  const int wall_model_type           = paras["wall_model_type"].as<int>();
  // wall_model_type: 0 no weakly enforced Dirichlet bc;
  //                  1 weakly enforced Dirichlet bc in all direction;
  //                  2 strongly enforced in wall-normal direction,
  //                   and weakly enforced in wall-tangent direction

  // Rotated paras:
  const std::vector<double> vec_point_rotated     = paras["point_rotated"].as<std::vector<double>>();
  const std::vector<double> vec_angular_direction = paras["angular_direction"].as<std::vector<double>>();

  SYS_T::print_fatal_if(VEC_T::get_size(vec_point_rotated) != 3, "Error: the size of the input point_rotated vector is not equal to 3. \n");
  SYS_T::print_fatal_if(VEC_T::get_size(vec_angular_direction) != 3, "Error: the size of the input angular_direction vector is not equal to 3. \n");

  // Info of rotation axis
  const Vector_3 point_rotated (vec_point_rotated[0], vec_point_rotated[1], vec_point_rotated[2]);
  const Vector_3 angular_direction = Vec3::normalize(Vector_3(vec_angular_direction[0], vec_angular_direction[1], vec_angular_direction[2]));

  SYS_T::print_fatal_if(std::isnan(angular_direction.x()) || std::isnan(angular_direction.y()) || std::isnan(angular_direction.z()), "Error: the direction vector of rotation axis cannot be zero vector. \n" );

  if( elemType != FEType::Tet4 && elemType != FEType::Tet10 && elemType != FEType::Hex8 && elemType != FEType::Hex27 ) SYS_T::print_fatal("ERROR: unknown element type %s.\n", elemType_str.c_str());

  // Print the command line arguments
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -elem_type: "<<elemType_str<<endl;
  cout<<" -wall_model_type: "<<wall_model_type<<endl;
  cout<<" -num_outlet: "<<num_outlet<<endl;
  cout<<" -fixed_geo_file: "<<fixed_geo_file<<endl;
  cout<<" -rotated_geo_file: "<<fixed_geo_file<<endl;
  cout<<" -sur_file_in_base: "<<sur_file_in_base<<endl;
  cout<<" -sur_file_inner_wall: "<<sur_file_inner_wall<<endl;
  cout<<" -sur_file_outer_wall: "<<sur_file_outer_wall<<endl;
  cout<<" -sur_file_out_base: "<<sur_file_out_base<<endl;
  cout<<" -fixed_interface_base: "<<fixed_interface_base<<endl;
  cout<<" -rotated_interface_base: "<<rotated_interface_base<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -isDualGraph: true \n";
  else cout<<" -isDualGraph: false \n";
  cout<<"---- Problem definition ----\n";
  cout<<" dofNum: "<<dofNum<<endl;
  cout<<" dofMat: "<<dofMat<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // Check if the vtu geometry files exist on disk
  SYS_T::file_check(fixed_geo_file); cout<<fixed_geo_file<<" found. \n";

  SYS_T::file_check(rotated_geo_file); cout<<rotated_geo_file<<" found. \n";

  SYS_T::file_check(sur_file_inner_wall); cout<<sur_file_outer_wall<<" found. \n";

  SYS_T::file_check(sur_file_outer_wall); cout<<sur_file_inner_wall<<" found. \n";

  // Generate the inlet file names and check existance
  std::vector< std::string > sur_file_in;
  sur_file_in.resize( num_inlet );

  for(int ii=0; ii<num_inlet; ++ii)
  {  
    if(elemType == FEType::Tet4 || elemType == FEType::Hex8)
      sur_file_in[ii] = SYS_T::gen_capfile_name( sur_file_in_base, ii, ".vtp" );   
    else if(elemType == FEType::Tet10 || elemType == FEType::Hex27)
      sur_file_in[ii] = SYS_T::gen_capfile_name( sur_file_in_base, ii, ".vtu" );
    else
      SYS_T::print_fatal("Error: unknown element type occurs when generating the inlet file names. \n"); 
  
    SYS_T::file_check(sur_file_in[ii]);
    cout<<sur_file_in[ii]<<" found. \n";
  }

  // Generate the outlet file names and check existance
  std::vector< std::string > sur_file_out;
  sur_file_out.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    if(elemType == FEType::Tet4 || elemType == FEType::Hex8)
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtp" ); 
    else if(elemType == FEType::Tet10 || elemType == FEType::Hex27)
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtu" ); 
    else
      SYS_T::print_fatal("Error: unknown element type occurs when generating the outlet file names. \n");

    SYS_T::file_check(sur_file_out[ii]);
    cout<<sur_file_out[ii]<<" found. \n";
  }

  std::vector< std::string > fixed_interface_file(num_interface_pair);
  std::vector< std::string > rotated_interface_file(num_interface_pair);
  for(int ii=0; ii<num_interface_pair; ++ii)
  {
    if(elemType == FEType::Tet4 || elemType == FEType::Hex8)
    {
      fixed_interface_file[ii] = SYS_T::gen_capfile_name( fixed_interface_base, ii, ".vtp" );
      rotated_interface_file[ii] = SYS_T::gen_capfile_name( rotated_interface_base, ii, ".vtp" );
    } 
    else if(elemType == FEType::Tet10 || elemType == FEType::Hex27)
    {
      fixed_interface_file[ii] = SYS_T::gen_capfile_name( fixed_interface_base, ii, ".vtu" );
      rotated_interface_file[ii] = SYS_T::gen_capfile_name( rotated_interface_base, ii, ".vtu" );
    }  
    else
      SYS_T::print_fatal("Error: unknown element type occurs when generating the outlet file names. \n");

    SYS_T::file_check(fixed_interface_file[ii]);
    cout<<fixed_interface_file[ii]<<" found. \n";

    SYS_T::file_check(rotated_interface_file[ii]);
    cout<<rotated_interface_file[ii]<<" found. \n";
  }

  // Record the problem setting into a HDF5 file: preprocessor_cmd.h5
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("num_inlet", num_inlet);
  cmdh5w->write_intScalar("num_outlet", num_outlet);
  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_string("elemType", elemType_str);
  cmdh5w->write_string("fixed_geo_file", fixed_geo_file);
  cmdh5w->write_string("rotated_geo_file", rotated_geo_file);
  cmdh5w->write_string("sur_file_in_base", sur_file_in_base);
  cmdh5w->write_string("sur_file_out_base", sur_file_out_base);
  cmdh5w->write_string("sur_file_inner_wall", sur_file_inner_wall);
  cmdh5w->write_string("sur_file_outer_wall", sur_file_outer_wall);
  cmdh5w->write_string("fixed_interface_base", fixed_interface_base);
  cmdh5w->write_string("rotated_interface_base", rotated_interface_base);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);

  // Read the volumetric mesh file from the vtu file: fixed_geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;
  
  VTK_T::read_vtu_grid(fixed_geo_file, nFunc, nElem, ctrlPts, vecIEN);
  const int fixed_nFunc = nFunc, fixed_nElem = nElem;

  int rotated_nFunc, rotated_nElem;
  std::vector<int> rotated_vecIEN;
  std::vector<double> rotated_ctrlPts;

  VTK_T::read_vtu_grid(rotated_geo_file, rotated_nFunc, rotated_nElem, rotated_ctrlPts, rotated_vecIEN);
  nFunc += rotated_nFunc;
  nElem += rotated_nElem;

  // fixed_geo: tag = 0
  // rotated_geo: tag = 1
  std::vector<int> rotated_tag (nElem, 1);
  
  for (int ee=0; ee < fixed_nElem; ++ee)
    rotated_tag [ee] = 0;

  for (int &nodeid : rotated_vecIEN)
    nodeid += fixed_nFunc;

  VEC_T::insert_end(vecIEN, rotated_vecIEN);
  VEC_T::insert_end(ctrlPts, rotated_ctrlPts);

  IIEN * IEN = new IEN_FEM(nElem, vecIEN);
  VEC_T::clean( vecIEN ); // clean the vector
  VEC_T::clean( rotated_vecIEN );
  VEC_T::clean( rotated_ctrlPts );

  // Generate the list of fixed and rotated nodes
  std::vector<int> node_f = VTK_T::read_int_PointData( fixed_geo_file, "GlobalNodeID" );

  std::vector<int> node_r = VTK_T::read_int_PointData( rotated_geo_file, "GlobalNodeID" );
  
  for (int &nodeid : node_r)
    nodeid += fixed_nFunc;

  VEC_T::sort_unique_resize( node_f ); VEC_T::sort_unique_resize( node_r );

  IMesh * mesh = nullptr;

  switch( elemType )
  {
    case FEType::Tet4:
      mesh = new Mesh_Tet(nFunc, nElem, 1);
      break;
    case FEType::Tet10:
      mesh = new Mesh_Tet(nFunc, nElem, 2);
      break;
    case FEType::Hex8:
      mesh = new Mesh_FEM(nFunc, nElem, 8, 1);
      break;
    case FEType::Hex27:
      mesh = new Mesh_FEM(nFunc, nElem, 27, 2);
      break;      
    default:
      SYS_T::print_fatal("Error: elemType %s is not supported.\n", elemType_str.c_str());
      break;
  }

  SYS_T::print_fatal_if( IEN->get_nLocBas() != mesh->get_nLocBas(), "Error: the nLocBas from the Mesh %d and the IEN %d classes do not match. \n", mesh->get_nLocBas(), IEN->get_nLocBas() );

  mesh -> print_info();
  
  // Call METIS to partition the mesh 
  IGlobal_Part * global_part = nullptr;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  // Generate the new nodal numbering
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");

  // Partition the interfaces
  for(int ii=0; ii < num_interface_pair; ++ii)
  {
    int sur_fixed_nFunc, sur_fixed_nElem, sur_rotated_nFunc, sur_rotated_nElem;
    std::vector<int> sur_fixed_vecIEN, sur_rotated_vecIEN;
    std::vector<double> sur_fixed_ctrlPts, sur_rotated_ctrlPts;

    VTK_T::read_grid(fixed_interface_file[ii], sur_fixed_nFunc, sur_fixed_nElem, sur_fixed_ctrlPts, sur_fixed_vecIEN);
    VTK_T::read_grid(rotated_interface_file[ii], sur_rotated_nFunc, sur_rotated_nElem, sur_rotated_ctrlPts, sur_rotated_vecIEN);

    IIEN * sur_fixed_IEN = new IEN_FEM(sur_fixed_nElem, sur_fixed_vecIEN);
    VEC_T::clean(sur_fixed_vecIEN);

    IIEN * sur_rotated_IEN = new IEN_FEM(sur_rotated_nElem, sur_rotated_vecIEN);
    VEC_T::clean(sur_rotated_vecIEN);

    IMesh * sur_fixed_mesh = nullptr;
    IMesh * sur_rotated_mesh = nullptr;

    switch( elemType )
    {
      case FEType::Tet4:
        sur_fixed_mesh = new Mesh_FEM(sur_fixed_nFunc, sur_fixed_nElem, 3, 1);
        sur_rotated_mesh = new Mesh_FEM(sur_rotated_nFunc, sur_rotated_nElem, 3, 1);
        break;
      case FEType::Tet10:
        sur_fixed_mesh = new Mesh_FEM(sur_fixed_nFunc, sur_fixed_nElem, 6, 2);
        sur_rotated_mesh = new Mesh_FEM(sur_rotated_nFunc, sur_rotated_nElem, 6, 2);
        break;
      case FEType::Hex8:
        sur_fixed_mesh = new Mesh_FEM(sur_fixed_nFunc, sur_fixed_nElem, 4, 1);
        sur_rotated_mesh = new Mesh_FEM(sur_rotated_nFunc, sur_rotated_nElem, 4, 1);
        break;
      case FEType::Hex27:
        sur_fixed_mesh = new Mesh_FEM(sur_fixed_nFunc, sur_fixed_nElem, 9, 2);
        sur_rotated_mesh = new Mesh_FEM(sur_rotated_nFunc, sur_rotated_nElem, 9, 2);
        break;      
      default:
        SYS_T::print_fatal("Error: elemType %s is not supported.\n", elemType_str.c_str());
        break;
    }
    
    std::string epart_base = "epart_", npart_base = "npart_";
    std::string fixed_epart = SYS_T::gen_capfile_name(epart_base, ii, "_fixed_itf");
    std::string fixed_npart = SYS_T::gen_capfile_name(npart_base, ii, "_fixed_itf");
    std::string rotated_epart = SYS_T::gen_capfile_name(epart_base, ii, "_rotated_itf");
    std::string rotated_npart = SYS_T::gen_capfile_name(npart_base, ii, "_rotated_itf");

    IGlobal_Part * global_part_fixed_itf = nullptr;
    IGlobal_Part * global_part_rotated_itf = nullptr;
    if(cpu_size > 1)
    {
      global_part_fixed_itf = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, sur_fixed_mesh, sur_fixed_IEN, fixed_epart, fixed_npart );

      global_part_rotated_itf = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, sur_rotated_mesh, sur_rotated_IEN, rotated_epart, rotated_npart );
    }
    else if(cpu_size == 1)
    {
      global_part_fixed_itf = new Global_Part_Serial( sur_fixed_mesh, fixed_epart, fixed_npart );

      global_part_rotated_itf = new Global_Part_Serial( sur_rotated_mesh, rotated_epart, rotated_npart );
    }
    else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

    delete global_part_fixed_itf; delete global_part_rotated_itf;
    delete sur_fixed_mesh; delete sur_rotated_mesh; delete sur_fixed_IEN; delete sur_rotated_IEN;
  }

  // Setup Nodal i.e. Dirichlet type Boundary Conditions
  std::vector<INodalBC *> NBC_list( dofMat, nullptr );

  std::vector<std::string> dir_list {};
  std::vector<std::string> weak_list {};

  for(int ii=0; ii<num_inlet; ++ii)
    dir_list.push_back( sur_file_in[ii] );
  
  if (wall_model_type == 0)
  {
    dir_list.push_back( sur_file_outer_wall );
  }
  else if (wall_model_type == 1 || wall_model_type == 2)
  {
    weak_list.push_back( sur_file_outer_wall );
  }
  else
    SYS_T::print_fatal("Unknown wall model type.");

  NBC_list[0] = new NodalBC( nFunc );
  NBC_list[1] = new NodalBC( dir_list, rotated_sur_file, sur_file_inner_wall, fixed_geo_file, nFunc );
  NBC_list[2] = new NodalBC( dir_list, rotated_sur_file, sur_file_inner_wall, fixed_geo_file, nFunc );
  NBC_list[3] = new NodalBC( dir_list, rotated_sur_file, sur_file_inner_wall, fixed_geo_file, nFunc );

  // Rotated BC info
  INodalBC * RotBC = new NodalBC_3D_rotated( rotated_sur_file, fixed_geo_file,
      nFunc, elemType );  

  // Inflow BC info
  std::vector< Vector_3 > inlet_outvec( sur_file_in.size() );

  if(elemType == FEType::Tet4 || elemType == FEType::Tet10)
  {
    for(unsigned int ii=0; ii<sur_file_in.size(); ++ii)
      inlet_outvec[ii] = TET_T::get_out_normal( sur_file_in[ii], ctrlPts, IEN );  
  }
  else if(elemType == FEType::Hex8 || elemType == FEType::Hex27)
  {
    for(unsigned int ii=0; ii<sur_file_in.size(); ++ii)
      inlet_outvec[ii] = HEX_T::get_out_normal( sur_file_in[ii], ctrlPts, IEN );  
  }
  else
    SYS_T::print_fatal("Error: unknown element type occurs when obtaining the outward normal vector for the inflow boundary condition. \n");

  INodalBC * InFBC = new NodalBC_3D_inflow( sur_file_in, sur_file_outer_wall,
      nFunc, inlet_outvec, elemType );

  InFBC -> resetSurIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations

  // Setup Elemental Boundary Conditions
  // Obtain the outward normal vector
  std::vector< Vector_3 > outlet_outvec( sur_file_out.size() );
  
  if(elemType == FEType::Tet4 || elemType == FEType::Tet10)
  {
    for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
      outlet_outvec[ii] = TET_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN );
  }
  else if(elemType == FEType::Hex8 || elemType == FEType::Hex27)
  {
    for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
      outlet_outvec[ii] = HEX_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN );
  }
  else
    SYS_T::print_fatal("Error: unknown element type occurs when obtaining the outward normal vector for the elemental boundary conditions. \n");

  ElemBC * ebc = new ElemBC_3D_outflow( sur_file_out, outlet_outvec, elemType );

  ebc -> resetSurIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations

  // Setup weakly enforced Dirichlet BC on wall if wall_model_type > 0
  ElemBC * wbc = new ElemBC_3D_turbulence_wall_model( weak_list, wall_model_type, IEN, elemType );

  // Set up interface info
  std::vector<double> intervals_0 {0.0, 6.0};

  Interface_pair itf_0(fixed_interface_file[0], rotated_interface_file[0], "epart_000_fixed_itf.h5", "epart_000_rotated_itf.h5",
    fixed_nElem, fixed_nFunc, ctrlPts, IEN, elemType, intervals_0, Vector_3(18.5, 0.0, 0.0));

  std::vector<double> intervals_1 {-4.5, 4.5};

  Interface_pair itf_1(fixed_interface_file[1], rotated_interface_file[1], "epart_001_fixed_itf.h5", "epart_001_rotated_itf.h5",
    fixed_nElem, fixed_nFunc, ctrlPts, IEN, elemType, intervals_1, 0);

  std::vector<Interface_pair> interfaces {itf_0, itf_1};
 
  // Start partition the mesh for each cpu_rank 
  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0; // total number of ghost nodes

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  // Shared data for interfaces
  std::vector<std::vector<std::vector<int>>> distributed_fixed_node_vol_part_tag;
  distributed_fixed_node_vol_part_tag.resize(cpu_size);

  std::vector<std::vector<std::vector<int>>> distributed_fixed_node_loc_pos;
  distributed_fixed_node_loc_pos.resize(cpu_size);

  std::vector<std::vector<std::vector<int>>> distributed_rotated_node_vol_part_tag;
  distributed_rotated_node_vol_part_tag.resize(cpu_size);

  std::vector<std::vector<std::vector<int>>> distributed_rotated_node_loc_pos;
  distributed_rotated_node_loc_pos.resize(cpu_size);

  std::vector<int> max_fixed_nlocalele (num_interface_pair, 0);

  std::vector<int> max_rotated_nlocalele(num_interface_pair, 0);

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();

    // IPart * part = new Part_FEM( mesh, global_part, mnindex, IEN,
    //     ctrlPts, rotated_tag, proc_rank, cpu_size, elemType, {0, dofNum, true, "NS"} );

    IPart * part = new Part_FEM_Rotated( mesh, global_part, mnindex, IEN,
        ctrlPts, rotated_tag, node_f, node_r, proc_rank, cpu_size, elemType, 
        {0, dofNum, true, "ROTATED_NS"} );
    
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file );

    part -> print_part_loadbalance_edgecut();
    
    // Partition Nodal BC and write to h5 file
    NBC_Partition * nbcpart = new NBC_Partition(part, mnindex, NBC_list);
    
    nbcpart -> write_hdf5( part_file );

    //Partition Nodal Rotated BC and write to h5 file
    NBC_Partition_rotated * rotpart = new NBC_Partition_rotated(part, mnindex, RotBC);

    rotpart -> write_hdf5( part_file );

    // Partition Nodal Inflow BC and write to h5 file
    NBC_Partition_inflow * infpart = new NBC_Partition_inflow(part, mnindex, InFBC);
    
    infpart -> write_hdf5( part_file );
    
    // Partition Elemental BC and write to h5 file
    EBC_Partition * ebcpart = new EBC_Partition_outflow(part, mnindex, ebc, NBC_list);

    ebcpart -> write_hdf5( part_file );

    // Partition Weak BC and write to h5 file
    EBC_Partition * wbcpart = new EBC_Partition_turbulence_wall_model(part, mnindex, wbc);

    wbcpart -> write_hdf5( part_file );

    // Writed the info of rotation axis into h5 file
    const std::string fName = SYS_T::gen_partfile_name( part_file, part->get_cpu_rank() );
    hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t g_id = H5Gcreate(file_id, "/rotation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer( file_id );
    h5w -> write_Vector_3( g_id, "point_rotated", point_rotated.to_std_array() );
    h5w -> write_Vector_3( g_id, "angular_direction", angular_direction.to_std_array() );

    delete h5w; H5Gclose( g_id ); H5Fclose( file_id );

    // Partition sliding interface and write to h5 file
    Interface_Partition * itfpart = new Interface_Partition(part, mnindex, interfaces, NBC_list);

    distributed_fixed_node_vol_part_tag[proc_rank] = itfpart -> get_fixed_node_vol_part_tag();
    distributed_fixed_node_loc_pos[proc_rank] = itfpart -> get_fixed_node_loc_pos();

    distributed_rotated_node_vol_part_tag[proc_rank] = itfpart -> get_rotated_node_vol_part_tag();
    distributed_rotated_node_loc_pos[proc_rank] = itfpart -> get_rotated_node_loc_pos();

    for(int ii = 0; ii < VEC_T::get_size(interfaces); ++ii)
    {
      if(max_fixed_nlocalele[ii] < itfpart -> get_fixed_nlocalele(ii))
        max_fixed_nlocalele[ii] = itfpart -> get_fixed_nlocalele(ii);

      if(max_rotated_nlocalele[ii] < itfpart -> get_rotated_nlocalele(ii))
        max_rotated_nlocalele[ii] = itfpart ->get_rotated_nlocalele(ii);
    }

    itfpart -> write_hdf5( part_file );

    // Collect partition statistics
    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete rotpart; delete infpart; delete ebcpart; delete wbcpart; delete itfpart;
  }

  // Combine the fixed/rotated_node_vol_part_tag and rotated_node_loc_pos
  std::vector<std::vector<int>> fixed_node_vol_part_tag, fixed_node_loc_pos;
  fixed_node_vol_part_tag.resize(VEC_T::get_size(interfaces));
  fixed_node_loc_pos.resize(VEC_T::get_size(interfaces));

  std::vector<std::vector<int>> rotated_node_vol_part_tag, rotated_node_loc_pos;
  rotated_node_vol_part_tag.resize(VEC_T::get_size(interfaces));
  rotated_node_loc_pos.resize(VEC_T::get_size(interfaces));

  for(int ii = 0; ii < VEC_T::get_size(interfaces); ++ii)
  { 
    // just a initialization
    fixed_node_vol_part_tag[ii] = distributed_fixed_node_vol_part_tag[0][ii];
    fixed_node_loc_pos[ii] = distributed_fixed_node_loc_pos[0][ii];

    rotated_node_vol_part_tag[ii] = distributed_rotated_node_vol_part_tag[0][ii];
    rotated_node_loc_pos[ii] = distributed_rotated_node_loc_pos[0][ii];

    for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
    {
      PERIGEE_OMP_PARALLEL_FOR
      for(int jj = 0; jj < VEC_T::get_size(fixed_node_vol_part_tag[ii]); ++jj)
      {
        if(distributed_fixed_node_vol_part_tag[proc_rank][ii][jj] != -1)
        {
          fixed_node_vol_part_tag[ii][jj] = distributed_fixed_node_vol_part_tag[proc_rank][ii][jj];
          fixed_node_loc_pos[ii][jj] = distributed_fixed_node_loc_pos[proc_rank][ii][jj];
        }
      }

      PERIGEE_OMP_PARALLEL_FOR
      for(int jj = 0; jj < VEC_T::get_size(rotated_node_vol_part_tag[ii]); ++jj)
      {
        if(distributed_rotated_node_vol_part_tag[proc_rank][ii][jj] != -1)
        {
          rotated_node_vol_part_tag[ii][jj] = distributed_rotated_node_vol_part_tag[proc_rank][ii][jj];
          rotated_node_loc_pos[ii][jj] = distributed_rotated_node_loc_pos[proc_rank][ii][jj];
        }
      }
    }
  }

  // Write the .h5 file
  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    const std::string fName = SYS_T::gen_partfile_name( part_file, proc_rank );

    const std::string GroupName = "/sliding";

    hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    hid_t g_id = H5Gopen( file_id, GroupName.c_str(), H5P_DEFAULT );

    HDF5_Writer * h5w = new HDF5_Writer( file_id );

    h5w -> write_intVector( g_id, "max_num_local_fixed_cell", max_fixed_nlocalele );

    h5w -> write_intVector( g_id, "max_num_local_rotated_cell", max_rotated_nlocalele );

    const std::string groupbase("interfaceid_");

    for(int ii = 0; ii < VEC_T::get_size(interfaces); ++ii)
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );

      hid_t group_id = H5Gopen(g_id, subgroup_name.c_str(), H5P_DEFAULT);

      h5w -> write_intVector( group_id, "fixed_node_part_tag", fixed_node_vol_part_tag[ii] );

      h5w -> write_intVector( group_id, "fixed_node_loc_pos", fixed_node_loc_pos[ii] );

      h5w -> write_intVector( group_id, "rotated_node_part_tag", rotated_node_vol_part_tag[ii] );

      h5w -> write_intVector( group_id, "rotated_node_loc_pos", rotated_node_loc_pos[ii] );

      H5Gclose( group_id );
    }

    delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
  }

  cout<<"\n===> Mesh Partition Quality: "<<endl;
  cout<<"The largest ghost / local node ratio is: "<<VEC_T::max(list_ratio_g2l)<<endl;
  cout<<"The smallest ghost / local node ratio is: "<<VEC_T::min(list_ratio_g2l)<<endl;
  cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<endl;
  cout<<"The maximum badnode number is: "<<VEC_T::max(list_nbadnode)<<endl;

  const int maxpart_nlocalnode = VEC_T::max(list_nlocalnode); 
  const int minpart_nlocalnode = VEC_T::min(list_nlocalnode);

  cout<<"The maximum and minimum local node numbers are ";
  cout<<maxpart_nlocalnode<<"\t"<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of local node is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;

  // Finalize the code and exit
  for(auto &it_nbc : NBC_list) delete it_nbc;

  delete InFBC; delete RotBC; delete ebc; delete wbc; delete mytimer;
  delete mnindex; delete global_part; delete mesh; delete IEN;

  return EXIT_SUCCESS;
}

// EOF
