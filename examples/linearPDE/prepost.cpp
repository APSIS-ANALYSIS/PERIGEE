// ============================================================================
// prepost.cpp
// 
// This is a prepost driver for handling linear differential equations.
//
// Author: Chi Ding
// Email:  am-dight@outlook.com 
// ============================================================================
#include "HDF5_Reader.hpp"
#include "VTK_Tools.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "yaml-cpp/yaml.h"

int main( int argc, char * argv[] )
{
  // Set number of threads and  print info of OpenMP
  SYS_T::print_omp_info();
  SYS_T::set_omp_num_threads();
  
  // Clean the existing part hdf5 files
  SYS_T::execute("rm -rf ppart");
  SYS_T::execute("mkdir ppart");

  // Define the partition file name
  const std::string part_file("./ppart/part");

  int cpu_size = 1;
  bool isDualGraph = true;
  
  // Read the problem setting recorded in the .h5 file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  std::string geo_file = cmd_h5r -> read_string("/", "geo_file");
  const std::string elemType_str  = cmd_h5r -> read_string("/","elemType");
  const int dofNum     = cmd_h5r -> read_intScalar("/","dof_num");
  const int dofMat     = cmd_h5r -> read_intScalar("/","dof_mat");
  int in_ncommon       = cmd_h5r -> read_intScalar("/","in_ncommon");
  const FEType elemType = FE_T::to_FEType(elemType_str);

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // The user can specify the new mesh partition options from the yaml file
  const std::string yaml_file("prepost.yml");

  SYS_T::file_check(yaml_file); std::cout << yaml_file << " found. \n";

  YAML::Node paras = YAML::LoadFile( yaml_file );

  cpu_size     = paras["cpu_size"].as<int>(); 
  in_ncommon   = paras["in_ncommon"].as<int>();
  isDualGraph  = paras["is_dualgraph"].as<bool>();

  cout << "==== Command Line Arguments ====" << endl;
  cout << " -cpu_size: "   << cpu_size   << endl;
  cout << " -in_ncommon: " << in_ncommon << endl;
  if(isDualGraph) cout << " -METIS_isDualGraph: true \n";
  else cout << " -METIS_isDualGraph: false \n";
  cout << "----------------------------------\n";
  cout << "-part_file: " << part_file << endl;
  cout << "-geo_file: "  << geo_file  << endl;
  cout << "-elemType: "  << elemType_str  << endl;
  cout << "-dof_num: "   << dofNum    << endl;
  cout << "-dof_mat: "   << dofMat    << endl; 

  // Read the volumetric mesh file from the vtu file: geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);

  IIEN * IEN = new IEN_FEM(nElem, vecIEN);
  VEC_T::clean( vecIEN ); // clean the vector

  const int nLocBas = FE_T::to_nLocBas(elemType);
  
  SYS_T::print_fatal_if( IEN->get_nLocBas() != nLocBas, "Error: the nLocBas from the Mesh %d and the IEN %d classes do not match. \n", nLocBas, IEN->get_nLocBas()); 

  // Call METIS to partition the mesh
  IGlobal_Part * global_part = nullptr;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, nElem, nFunc, nLocBas, IEN, "post_epart", "post_npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( nElem, nFunc, "post_epart", "post_npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, nFunc);
  mnindex->write_hdf5("post_node_mapping");

  cout<<"=== Start Partition ... \n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer -> Reset();
    mytimer -> Start();

    IPart * part = new Part_FEM( nElem, nFunc, nLocBas, global_part, mnindex, IEN,
        ctrlPts, proc_rank, cpu_size, elemType, {0, dofNum, true, "linearPDE"} );

    part -> write(part_file.c_str());
    mytimer -> Stop();
    cout << "-- proc " << proc_rank << " Time taken: " << mytimer -> get_sec() << " sec. \n";
    delete part;
  }

  delete mytimer; delete global_part; delete mnindex; delete IEN;
  return EXIT_SUCCESS;
}

// EOF
