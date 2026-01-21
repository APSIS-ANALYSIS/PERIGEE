// ============================================================================
// preprocess.cpp
// 
// This is a preprocessor code for handling solid mechanics problems using the
// stabilized mixed formulation.
//
// Date Created: Jan 04 2026
// ============================================================================
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "NodalBC.hpp"
#include "ElemBC_3D.hpp"
#include "NBC_Partition.hpp"
#include "EBC_Partition.hpp"
#include "yaml-cpp/yaml.h"

int main( int argc, char * argv[] )
{
  // Set number of threads and  print info of OpenMP
  SYS_T::print_omp_info();
  SYS_T::set_omp_num_threads();

  // Clean the potentially pre-existing hdf5 files in the job folder
  SYS_T::execute("rm -rf preprocessor_cmd.h5");
  SYS_T::execute("rm -rf apart");
  SYS_T::execute("mkdir apart");

  // Define basic problem settins
  constexpr int dofNum = 7; // degree-of-freedom for the physical problem
  constexpr int dofMat = 4; // degree-of-freedom in the matrix problem

  // Yaml options
  const std::string yaml_file("preprocess.yml");

  // Check if the yaml file exist on disk
  SYS_T::file_check(yaml_file);

  YAML::Node paras = YAML::LoadFile( yaml_file );

  const std::string elemType_str = paras["elem_type"].as<std::string>();
  const std::string geo_file     = paras["geo_file"].as<std::string>();
  const std::string part_file    = paras["part_file"].as<std::string>();
  const int cpu_size             = paras["cpu_size"].as<int>();
  const int in_ncommon           = paras["in_ncommon"].as<int>();
  const bool isDualGraph         = paras["is_dualgraph"].as<bool>();
  const FEType elemType          = FE_T::to_FEType(elemType_str);
  const std::vector<std::string> sur_file_neu   = paras["sur_file_neu"].as<std::vector<std::string>>();
  const std::vector<std::string> sur_file_dir_x = paras["sur_file_dir_x"].as<std::vector<std::string>>();
  const std::vector<std::string> sur_file_dir_y = paras["sur_file_dir_y"].as<std::vector<std::string>>();
  const std::vector<std::string> sur_file_dir_z = paras["sur_file_dir_z"].as<std::vector<std::string>>();

  if(elemType!=FEType::Tet4 && elemType!=FEType::Tet10 && elemType!=FEType::Hex8 && elemType!=FEType::Hex27) SYS_T::print_fatal("ERROR: unknown element type %s.\n", elemType_str.c_str());

  // Print the command line arguments
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -elem_type: "<<elemType_str<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -sur_file_dir_x: ";
  for(const auto &fname : sur_file_dir_x) cout<<fname<<" ";
  cout<<endl;
  cout<<" -sur_file_dir_y: ";
  for(const auto &fname : sur_file_dir_y) cout<<fname<<" ";
  cout<<endl;
  cout<<" -sur_file_dir_z: ";
  for(const auto &fname : sur_file_dir_z) cout<<fname<<" ";
  cout<<endl;
  cout<<" -sur_file_neu: ";
  for(const auto &fname : sur_file_neu) cout<<fname<<" ";
  cout<<endl;
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
  SYS_T::file_check(geo_file);
  std::cout << geo_file << " found. \n";
  for( const auto &fname : sur_file_dir_x )
  { SYS_T::file_check( fname ); std::cout << fname << " found. \n"; }
  for( const auto &fname : sur_file_dir_y )
  { SYS_T::file_check( fname ); std::cout << fname << " found. \n"; }
  for( const auto &fname : sur_file_dir_z )
  { SYS_T::file_check( fname ); std::cout << fname << " found. \n"; }
  for( const auto &fname : sur_file_neu )
  { SYS_T::file_check( fname ); std::cout << fname << " found. \n"; }

  // Record the problem setting into a HDF5 file: preprocessor_cmd.h5
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_string("elemType", elemType_str);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("part_file", part_file);
  cmdh5w->write_intScalar("dof_num", dofNum);
  cmdh5w->write_intScalar("dof_mat", dofMat);

  delete cmdh5w; H5Fclose(cmd_file_id);

  // Read the volumetric mesh file from the vtu file: geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);

  IIEN * IEN = new IEN_FEM(nElem, vecIEN);
  VEC_T::clean( vecIEN ); // clean the vector

  const int nLocBas = FE_T::to_nLocBas(elemType);

  SYS_T::print_fatal_if( IEN->get_nLocBas() != nLocBas, "Error: the nLocBas from the given element type is %d and the mesh file is %d, which do not match. \n", nLocBas, IEN->get_nLocBas() );

  // Call METIS to partition the mesh
  IGlobal_Part * global_part = nullptr;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, nElem, nFunc, nLocBas, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( nElem, nFunc, "epart", "npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  // Generate the new nodal numbering
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, nFunc);
  mnindex->write_hdf5("node_mapping");

  // Setup Nodal i.e. Dirichlet type Boundary Conditions
  std::vector<INodalBC *> NBC_list( dofMat, nullptr );

  NBC_list[0] = new NodalBC( nFunc );
  NBC_list[1] = new NodalBC( sur_file_dir_x, nFunc );
  NBC_list[2] = new NodalBC( sur_file_dir_y, nFunc );
  NBC_list[3] = new NodalBC( sur_file_dir_z, nFunc );

  ElemBC * ebc = new ElemBC_3D( sur_file_neu, elemType );
  ebc -> resetSurIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations

  // Start partitioning the mesh for each cpu rank
  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0; // total number of ghost nodes

  auto mytimer = SYS_T::make_unique<SYS_T::Timer>();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();
    auto part = SYS_T::make_unique<Part_FEM>(
        nElem, nFunc, nLocBas, global_part, mnindex, IEN,
        ctrlPts, proc_rank, cpu_size, elemType,
        Field_Property(0, dofNum, true, "Solid") );

    part -> print_part_loadbalance_edgecut();
    
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file );

    // Partition Nodal BC and write to h5 file
    auto nbcpart = SYS_T::make_unique<NBC_Partition>(part.get(), mnindex, NBC_list);

    nbcpart -> write_hdf5( part_file );

    // Partition Elemental BC and write to h5 file
    auto ebcpart = SYS_T::make_unique<EBC_Partition>(part.get(), mnindex, ebc);

    ebcpart -> write_hdf5( part_file );

    // Collect partition statistics
    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
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
  for(auto &it_nbc : NBC_list ) delete it_nbc;

  delete ebc; delete global_part; delete mnindex; delete IEN;

  return EXIT_SUCCESS;
}

// EOF
