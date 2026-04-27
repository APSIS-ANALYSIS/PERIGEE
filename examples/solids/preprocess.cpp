// ============================================================================
// preprocess.cpp
// 
// This is a preprocessor code for handling solid mechanics problems using the
// stabilized mixed formulation.
//
// Date Created: Jan 04 2026
// ============================================================================
#include "IEN_FEM.hpp"
#include "VTK_Tools.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "NodalBC_Solid.hpp"
#include "Part_FEM.hpp"
#include "ElemBC_3D.hpp"
#include "NBC_Partition_Solid.hpp"
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
  const std::vector<bool> is_disp_driven_x =
    paras["is_disp_driven_x"] ? paras["is_disp_driven_x"].as<std::vector<bool>>()
                              : std::vector<bool>(sur_file_dir_x.size(), false);
  const std::vector<bool> is_disp_driven_y =
    paras["is_disp_driven_y"] ? paras["is_disp_driven_y"].as<std::vector<bool>>()
                              : std::vector<bool>(sur_file_dir_y.size(), false);
  const std::vector<bool> is_disp_driven_z =
    paras["is_disp_driven_z"] ? paras["is_disp_driven_z"].as<std::vector<bool>>()
                              : std::vector<bool>(sur_file_dir_z.size(), false);

  SYS_T::print_fatal_if( sur_file_dir_x.size() != is_disp_driven_x.size(),
      "Error: sur_file_dir_x and is_disp_driven_x size mismatch.\n" );
  SYS_T::print_fatal_if( sur_file_dir_y.size() != is_disp_driven_y.size(),
      "Error: sur_file_dir_y and is_disp_driven_y size mismatch.\n" );
  SYS_T::print_fatal_if( sur_file_dir_z.size() != is_disp_driven_z.size(),
      "Error: sur_file_dir_z and is_disp_driven_z size mismatch.\n" );

  if(elemType!=FEType::Tet4 && elemType!=FEType::Tet10 && elemType!=FEType::Hex8 && elemType!=FEType::Hex27) SYS_T::print_fatal("ERROR: unknown element type %s.\n", elemType_str.c_str());

  // Print the command line arguments
  std::cout<<"==== Command Line Arguments ===="<<std::endl;
  std::cout<<" -elem_type: "<<elemType_str<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -sur_file_dir_x: ";
  for(const auto &fname : sur_file_dir_x) std::cout<<fname<<" ";
  std::cout<<std::endl;
  std::cout<<" -is_disp_driven_x: ";
  for(const auto &flag : is_disp_driven_x) std::cout<<(flag ? "true" : "false")<<" ";
  std::cout<<std::endl;
  std::cout<<" -sur_file_dir_y: ";
  for(const auto &fname : sur_file_dir_y) std::cout<<fname<<" ";
  std::cout<<std::endl;
  std::cout<<" -is_disp_driven_y: ";
  for(const auto &flag : is_disp_driven_y) std::cout<<(flag ? "true" : "false")<<" ";
  std::cout<<std::endl;
  std::cout<<" -sur_file_dir_z: ";
  for(const auto &fname : sur_file_dir_z) std::cout<<fname<<" ";
  std::cout<<std::endl;
  std::cout<<" -is_disp_driven_z: ";
  for(const auto &flag : is_disp_driven_z) std::cout<<(flag ? "true" : "false")<<" ";
  std::cout<<std::endl;
  std::cout<<" -sur_file_neu: ";
  for(const auto &fname : sur_file_neu) std::cout<<fname<<" ";
  std::cout<<std::endl;
  std::cout<<" -part_file: "<<part_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  if(isDualGraph) std::cout<<" -isDualGraph: true \n";
  else std::cout<<" -isDualGraph: false \n";
  std::cout<<"---- Problem definition ----\n";
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;

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
  {
    auto cmdh5w = SYS_T::make_unique<HDF5_Writer>("preprocessor_cmd.h5");
    cmdh5w->write_intScalar("cpu_size", cpu_size);
    cmdh5w->write_intScalar("in_ncommon", in_ncommon);
    cmdh5w->write_string("elemType", elemType_str);
    cmdh5w->write_string("geo_file", geo_file);
    cmdh5w->write_string("part_file", part_file);
    cmdh5w->write_intScalar("dof_num", dofNum);
    cmdh5w->write_intScalar("dof_mat", dofMat);
  }
  
  // Read the volumetric mesh file from the vtu file: geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);

  auto IEN = SYS_T::make_unique<IEN_FEM>(nElem, vecIEN);
  VEC_T::clean( vecIEN ); // clean the vector

  const int nLocBas = FE_T::to_nLocBas(elemType);

  SYS_T::print_fatal_if( IEN->get_nLocBas() != nLocBas, "Error: the nLocBas from the given element type is %d and the mesh file is %d, which do not match. \n", nLocBas, IEN->get_nLocBas() );

  // Call METIS to partition the mesh
  std::unique_ptr<IGlobal_Part> global_part = nullptr;
  if(cpu_size > 1)
    global_part = SYS_T::make_unique<Global_Part_METIS>( cpu_size, in_ncommon,
        isDualGraph, nElem, nFunc, nLocBas, IEN.get(), "epart", "npart" );
  else if(cpu_size == 1)
    global_part = SYS_T::make_unique<Global_Part_Serial>( nElem, nFunc, "epart", "npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  // Generate the new nodal numbering
  auto mnindex = SYS_T::make_unique<Map_Node_Index>(global_part.get(), cpu_size, nFunc);
  mnindex->write_hdf5("node_mapping");

  // Setup Nodal Boundary Conditions, grouped by direction.
  std::vector<NodalBC_Solid *> solid_nbc_list_x {};
  std::vector<NodalBC_Solid *> solid_nbc_list_y {};
  std::vector<NodalBC_Solid *> solid_nbc_list_z {};

  solid_nbc_list_x.reserve( sur_file_dir_x.size() );
  solid_nbc_list_y.reserve( sur_file_dir_y.size() );
  solid_nbc_list_z.reserve( sur_file_dir_z.size() );

  for(std::size_t ii=0; ii<sur_file_dir_x.size(); ++ii)
    solid_nbc_list_x.push_back(
        new NodalBC_Solid( sur_file_dir_x[ii], nFunc, is_disp_driven_x[ii] ) );

  for(std::size_t ii=0; ii<sur_file_dir_y.size(); ++ii)
    solid_nbc_list_y.push_back(
        new NodalBC_Solid( sur_file_dir_y[ii], nFunc, is_disp_driven_y[ii] ) );

  for(std::size_t ii=0; ii<sur_file_dir_z.size(); ++ii)
    solid_nbc_list_z.push_back(
        new NodalBC_Solid( sur_file_dir_z[ii], nFunc, is_disp_driven_z[ii] ) );

  auto ebc = SYS_T::make_unique<ElemBC_3D>( sur_file_neu, elemType );
  ebc -> resetSurIEN_outwardnormal( IEN.get() ); // reset IEN for outward normal calculations

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
        nElem, nFunc, nLocBas, global_part.get(), mnindex.get(), IEN.get(),
        ctrlPts, proc_rank, cpu_size, elemType,
        Field_Property(0, dofNum, true, "Solid") );

    part -> print_part_loadbalance_edgecut();
    
    mytimer->Stop();
    std::cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file );

    // Partition Nodal BC and write to h5 file
    auto nbcpart = SYS_T::make_unique<NBC_Partition_Solid>(
        part.get(), mnindex.get(),
        solid_nbc_list_x, solid_nbc_list_y, solid_nbc_list_z,
        dofMat, nFunc );
    nbcpart->write_hdf5( part_file );

    // Partition Elemental BC and write to h5 file
    auto ebcpart = SYS_T::make_unique<EBC_Partition>(part.get(), mnindex.get(), ebc.get());

    ebcpart -> write_hdf5( part_file );

    // Collect partition statistics
    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
  }

  std::cout<<"\n===> Mesh Partition Quality: "<<std::endl;
  std::cout<<"The largest ghost / local node ratio is: "<<VEC_T::max(list_ratio_g2l)<<std::endl;
  std::cout<<"The smallest ghost / local node ratio is: "<<VEC_T::min(list_ratio_g2l)<<std::endl;
  std::cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<std::endl;
  std::cout<<"The maximum badnode number is: "<<VEC_T::max(list_nbadnode)<<std::endl;

  const int maxpart_nlocalnode = VEC_T::max(list_nlocalnode);
  const int minpart_nlocalnode = VEC_T::min(list_nlocalnode);

  std::cout<<"The maximum and minimum local node numbers are ";
  std::cout<<maxpart_nlocalnode<<"\t"<<minpart_nlocalnode<<std::endl;
  std::cout<<"The maximum / minimum of local node is: ";
  std::cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<std::endl;

  // Finalize the code and exit
  for(auto &it_nbc : solid_nbc_list_x ) delete it_nbc;
  for(auto &it_nbc : solid_nbc_list_y ) delete it_nbc;
  for(auto &it_nbc : solid_nbc_list_z ) delete it_nbc;

  return EXIT_SUCCESS;
}

// EOF
