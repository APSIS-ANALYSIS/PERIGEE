// ============================================================================
// prepost_ns_tets.cpp
//
// This is the partitioning routine for parallel postprocessors.
//
// Date: Jan. 24 2017
// ============================================================================
#include "HDF5_Reader.hpp"
#include "VTK_Tools.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "yaml-cpp/yaml.h"
#include <memory>

int main( int argc, char * argv[] )
{
  // Set number of threads and  print info of OpenMP
  SYS_T::print_omp_info();
  SYS_T::set_omp_num_threads();
  
  // Clean the existing part hdf5 files
  int sysret = system("rm -rf postpart_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "ERROR: system call failed. \n");

  // Read preprocessor command-line arguements recorded in the .h5 file

  auto cmd_h5r = SYS_T::make_unique<HDF5_Reader>( "preprocessor_cmd.h5" );

  std::string geo_file = cmd_h5r -> read_string("/", "geo_file");
  const std::string elemType_str = cmd_h5r -> read_string("/","elemType");
  const int dofNum = cmd_h5r -> read_intScalar("/","dofNum");
  int in_ncommon = cmd_h5r -> read_intScalar("/","in_ncommon");
  const FEType elemType = FE_T::to_FEType(elemType_str);

  // The user can specify the new mesh partition options from the yaml file
  const std::string yaml_file("ns_prepost.yml");

  SYS_T::file_check(yaml_file);

  YAML::Node paras = YAML::LoadFile( yaml_file );

  const int cpu_size = paras["cpu_size"].as<int>();
  in_ncommon = paras["in_ncommon"].as<int>();
  const bool isDualGraph = paras["is_dualgraph"].as<bool>();  
  const std::string part_file = paras["part_file"].as<std::string>();

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -METIS_isDualGraph: true \n";
  else cout<<" -METIS_isDualGraph: false \n";
  cout<<"----------------------------------\n";
  cout<<"part_file: "<<part_file<<endl;
  cout<<"geo_file: "<<geo_file<<endl;
  cout<<"elemType: "<<elemType_str<<endl;
  cout<<"dof_num: "<<dofNum<<endl;

  // Read the geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  // Check if the given geo file exist
  SYS_T::file_check( geo_file );

  VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);
  
  auto IEN = SYS_T::make_unique<IEN_FEM>(nElem, vecIEN);
  VEC_T::clean( vecIEN ); // clean the vector

  const int nLocBas = FE_T::to_nLocBas(elemType);
  
  SYS_T::print_fatal_if( IEN->get_nLocBas() != nLocBas, "Error: the nLocBas from the Mesh %d and the IEN %d classes do not match. \n", nLocBas, IEN->get_nLocBas());

  std::unique_ptr<IGlobal_Part> global_part = nullptr;
  if(cpu_size > 1)
    global_part = SYS_T::make_unique<Global_Part_METIS>( cpu_size, in_ncommon,
        isDualGraph, nElem, nFunc, nLocBas, IEN.get(), "post_epart", "post_npart" );
  else if(cpu_size == 1)
    global_part = SYS_T::make_unique<Global_Part_Serial>( nElem, nFunc, "post_epart", "post_npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  auto mnindex = SYS_T::make_unique<Map_Node_Index>(global_part.get(), cpu_size, nFunc);
  mnindex->write_hdf5("post_node_mapping");

  cout<<"=== Start Partition ... \n";
  int proc_size = cpu_size;
  auto mytimer = SYS_T::make_unique<SYS_T::Timer>();
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    mytimer->Reset(); mytimer->Start();
    auto part = SYS_T::make_unique<Part_FEM>( nElem, nFunc, nLocBas, global_part.get(), mnindex.get(), IEN.get(),
        ctrlPts, proc_rank, proc_size, elemType, {0, dofNum, true, "NS"} );
    part->write(part_file.c_str());
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";
  }

  // Clean memory
  cout<<"=== Clean memory. \n";
  return EXIT_SUCCESS;
}

// EOF
