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
#include "Part_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "ElemBC_3D_turbulence_wall_model.hpp"
#include "Interface_pair.hpp"
#include "NBC_Partition.hpp"
#include "NBC_Partition_inflow.hpp"
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

  const int elemType                  = paras["elem_type"].as<int>();
  const int num_inlet                 = paras["num_inlet"].as<int>();
  const int num_outlet                = paras["num_outlet"].as<int>();
  const std::string geo_file          = paras["geo_file"].as<std::string>();
  const std::string sur_file_in_base  = paras["sur_file_in_base"].as<std::string>();
  const std::string sur_file_wall     = paras["sur_file_wall"].as<std::string>();
  const std::string sur_file_out_base = paras["sur_file_out_base"].as<std::string>();

  const int num_interface_pair        = paras["num_interface_pair"].as<int>();
  const std::string rotated_geo_file  = paras["rotated_geo_file"].as<std::string>();
  const std::string fixed_interface_base   = paras["fixed_interface_base"].as<std::string>();
  const std::string rotated_interface_base = paras["rotated_interface_base"].as<std::string>();

  const std::string part_file         = paras["part_file"].as<std::string>();
  const int cpu_size                  = paras["cpu_size"].as<int>();
  const int in_ncommon                = paras["in_ncommon"].as<int>();
  const bool isDualGraph              = paras["is_dualgraph"].as<bool>();

  // Optional:
  const int wall_model_type           = paras["wall_model_type"].as<int>();
  // wall_model_type: 0 no weakly enforced Dirichlet bc;
  //                  1 weakly enforced Dirichlet bc in all direction;
  //                  2 strongly enforced in wall-normal direction,
  //                   and weakly enforced in wall-tangent direction

  if( elemType != 501 && elemType != 502 && elemType != 601 && elemType != 602 ) SYS_T::print_fatal("ERROR: unknown element type %d.\n", elemType);

  // Print the command line arguments
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -elem_type: "<<elemType<<endl;
  cout<<" -wall_model_type: "<<wall_model_type<<endl;
  cout<<" -num_outlet: "<<num_outlet<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -sur_file_in_base: "<<sur_file_in_base<<endl;
  cout<<" -sur_file_wall: "<<sur_file_wall<<endl;
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
  SYS_T::file_check(geo_file); cout<<geo_file<<" found. \n";

  SYS_T::file_check(sur_file_wall); cout<<sur_file_wall<<" found. \n";

  // Generate the inlet file names and check existance
  std::vector< std::string > sur_file_in;
  sur_file_in.resize( num_inlet );

  for(int ii=0; ii<num_inlet; ++ii)
  {  
    if(elemType == 501 || elemType == 601)
      sur_file_in[ii] = SYS_T::gen_capfile_name( sur_file_in_base, ii, ".vtp" );   
    else if(elemType == 502 || elemType == 602)
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
    if(elemType == 501 || elemType == 601)
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtp" ); 
    else if(elemType == 502 || elemType == 602)
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtu" ); 
    else
      SYS_T::print_fatal("Error: unknown element type occurs when generating the outlet file names. \n");

    SYS_T::file_check(sur_file_out[ii]);
    cout<<sur_file_out[ii]<<" found. \n";
  }

  // Check if the vtu geometry files exist on disk
  SYS_T::file_check(rotated_geo_file); cout<<rotated_geo_file<<" found. \n";

  std::vector< std::string > fixed_interface_file(num_interface_pair);
  std::vector< std::string > rotated_interface_file(num_interface_pair);
  for(int ii=0; ii<num_interface_pair; ++ii)
  {
    if(elemType == 501 || elemType == 601)
    {
      fixed_interface_file[ii] = SYS_T::gen_capfile_name( fixed_interface_base, ii, ".vtp" );
      rotated_interface_file[ii] = SYS_T::gen_capfile_name( rotated_interface_base, ii, ".vtp" );
    } 
    else if(elemType == 502 || elemType == 602)
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
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("sur_file_in_base", sur_file_in_base);
  cmdh5w->write_string("sur_file_out_base", sur_file_out_base);
  cmdh5w->write_string("sur_file_wall", sur_file_wall);
  cmdh5w->write_string("fixed_interface_base", fixed_interface_base);
  cmdh5w->write_string("rotated_interface_base", rotated_interface_base);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);

  // Read the volumetric mesh file from the vtu file: geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;
  
  VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);
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
  
  IMesh * mesh = nullptr;

  switch( elemType )
  {
    case 501:
      mesh = new Mesh_Tet(nFunc, nElem, 1);
      break;
    case 502:
      mesh = new Mesh_Tet(nFunc, nElem, 2);
      break;
    case 601:
      mesh = new Mesh_FEM(nFunc, nElem, 8, 1);
      break;
    case 602:
      mesh = new Mesh_FEM(nFunc, nElem, 27, 2);
      break;      
    default:
      SYS_T::print_fatal("Error: elemType %d is not supported.\n", elemType);
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
    int s_nFunc, s_nElem;
    std::vector<int> s_vecIEN;
    std::vector<double> s_ctrlPts;

    VTK_T::read_grid(fixed_interface_file[ii], s_nFunc, s_nElem, s_ctrlPts, s_vecIEN);

    IIEN * sIEN = new IEN_FEM(s_nElem, s_vecIEN);
    VEC_T::clean(s_vecIEN);

    IMesh * s_mesh = nullptr;

    switch( elemType )
    {
      case 501:
        s_mesh = new Mesh_FEM(s_nFunc, s_nElem, 3, 1);
        break;
      case 502:
        s_mesh = new Mesh_FEM(s_nFunc, s_nElem, 6, 2);
        break;
      case 601:
        s_mesh = new Mesh_FEM(s_nFunc, s_nElem, 4, 1);
        break;
      case 602:
        s_mesh = new Mesh_FEM(s_nFunc, s_nElem, 9, 2);
        break;      
      default:
        SYS_T::print_fatal("Error: elemType %d is not supported.\n", elemType);
        break;
    }
    
    std::string epart_base = "epart_", npart_base = "npart_";
    std::string epart = SYS_T::gen_capfile_name(epart_base, ii, "_itf");
    std::string npart = SYS_T::gen_capfile_name(npart_base, ii, "_itf");

    IGlobal_Part * global_part_itf = nullptr;
    if(cpu_size > 1)
      global_part_itf = new Global_Part_METIS( cpu_size, in_ncommon,
          isDualGraph, s_mesh, sIEN, epart, npart );
    else if(cpu_size == 1)
      global_part_itf = new Global_Part_Serial( s_mesh, epart, npart );
    else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

    delete global_part_itf; delete s_mesh; delete sIEN;
  }

  // Setup Nodal i.e. Dirichlet type Boundary Conditions
  std::vector<INodalBC *> NBC_list( dofMat, nullptr );

  std::vector<std::string> dir_list {};
  std::vector<std::string> weak_list {};

  for(int ii=0; ii<num_inlet; ++ii)
    dir_list.push_back( sur_file_in[ii] );
  
  if (wall_model_type == 0)
    dir_list.push_back( sur_file_wall );
  else if (wall_model_type == 1 || wall_model_type == 2)
    weak_list.push_back( sur_file_wall );
  else
    SYS_T::print_fatal("Unknown wall model type.");

  NBC_list[0] = new NodalBC( nFunc );
  NBC_list[1] = new NodalBC( dir_list, nFunc );
  NBC_list[2] = new NodalBC( dir_list, nFunc );
  NBC_list[3] = new NodalBC( dir_list, nFunc );

  // Inflow BC info
  std::vector< Vector_3 > inlet_outvec( sur_file_in.size() );

  if(elemType == 501 || elemType == 502)
  {
    for(unsigned int ii=0; ii<sur_file_in.size(); ++ii)
      inlet_outvec[ii] = TET_T::get_out_normal( sur_file_in[ii], ctrlPts, IEN );    
  }
  else if(elemType == 601 || elemType == 602)
  {
    for(unsigned int ii=0; ii<sur_file_in.size(); ++ii)
      inlet_outvec[ii] = HEX_T::get_out_normal( sur_file_in[ii], ctrlPts, IEN );  
  }
  else
    SYS_T::print_fatal("Error: unknown element type occurs when obtaining the outward normal vector for the inflow boundary condition. \n");

  INodalBC * InFBC = new NodalBC_3D_inflow( sur_file_in, sur_file_wall,
      nFunc, inlet_outvec, elemType );

  InFBC -> resetSurIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations

  // Setup Elemental Boundary Conditions
  // Obtain the outward normal vector
  std::vector< Vector_3 > outlet_outvec( sur_file_out.size() );
  
  if(elemType == 501 || elemType == 502)
  {
    for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
      outlet_outvec[ii] = TET_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN );  
  }
  else if(elemType == 601 || elemType == 602)
  {
    for(unsigned int ii=0; ii<sur_file_in.size(); ++ii)
      outlet_outvec[ii] = HEX_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN );
  }
  else
    SYS_T::print_fatal("Error: unknown element type occurs when obtaining the outward normal vector for the elemental boundary conditions. \n");

  ElemBC * ebc = new ElemBC_3D_outflow( sur_file_out, outlet_outvec, elemType );

  ebc -> resetSurIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations

  // Setup weakly enforced Dirichlet BC on wall if wall_model_type > 0
  ElemBC * wbc = new ElemBC_3D_turbulence_wall_model( weak_list, wall_model_type, IEN, elemType );

  // Set up interface info
  std::vector<double> intervals_0 {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5};

  Interface_pair itf_0(fixed_interface_file[0], rotated_interface_file[0], "epart_000_itf.h5",
    fixed_nElem, fixed_nFunc, ctrlPts, IEN, elemType, intervals_0, 0);

  std::vector<double> intervals_12 {0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3};

  Interface_pair itf_1(fixed_interface_file[1], rotated_interface_file[1], "epart_001_itf.h5",
    fixed_nElem, fixed_nFunc, ctrlPts, IEN, elemType, intervals_12, Vector_3(0.5, 0.0, 0.0));

  Interface_pair itf_2(fixed_interface_file[2], rotated_interface_file[2], "epart_002_itf.h5",
    fixed_nElem, fixed_nFunc, ctrlPts, IEN, elemType, intervals_12, Vector_3(-0.5, 0.0, 0.0));

  std::vector<Interface_pair> interfaces {itf_0, itf_1, itf_2};
 
  // Start partition the mesh for each cpu_rank 

  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0; // total number of ghost nodes

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();

    IPart * part = new Part_FEM( mesh, global_part, mnindex, IEN,
        ctrlPts, rotated_tag, proc_rank, cpu_size, elemType, {0, dofNum, true, "NS"} );
    
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file );

    part -> print_part_loadbalance_edgecut();
    
    // Partition Nodal BC and write to h5 file
    NBC_Partition * nbcpart = new NBC_Partition(part, mnindex, NBC_list);
    
    nbcpart -> write_hdf5( part_file );

    // Partition Nodal Inflow BC and write to h5 file
    NBC_Partition_inflow * infpart = new NBC_Partition_inflow(part, mnindex, InFBC);
    
    infpart->write_hdf5( part_file );
    
    // Partition Elemental BC and write to h5 file
    EBC_Partition * ebcpart = new EBC_Partition_outflow(part, mnindex, ebc, NBC_list);

    ebcpart -> write_hdf5( part_file );

    // Partition Weak BC and write to h5 file
    EBC_Partition * wbcpart = new EBC_Partition_turbulence_wall_model(part, mnindex, wbc);

    wbcpart -> write_hdf5( part_file );

    // Partition sliding interface and write to h5 file
    Interface_Partition * itfpart = new Interface_Partition(part, mnindex, interfaces);

    itfpart -> write_hdf5( part_file );

    // Collect partition statistics
    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete infpart; delete ebcpart; delete wbcpart; delete itfpart;
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

  delete InFBC; delete ebc; delete wbc; delete mytimer;
  delete mnindex; delete global_part; delete mesh; delete IEN;

  return EXIT_SUCCESS;
}

// EOF