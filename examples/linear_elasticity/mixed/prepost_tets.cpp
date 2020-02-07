// ==================================================================
// prepost_tets.cpp
//
// This is the partitioning routine for parallel postprocessors.
//
// Date: Jan. 24 2017
// ==================================================================
#include "HDF5_Reader.hpp"
#include "Tet_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet.hpp"

int main( int argc, char * argv[] )
{
  std::string geo_file, sur_file, pre_part_file;
  int dofNum, elemType, in_ncommon, probDim;
  
  std::string part_file("postpart");
  int cpu_size = 1;
  bool isDualGraph = true;
  bool isread_part = true;

  PetscMPIInt rank, size;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
   
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::print_fatal_if(size!=1,
      "ERROR: preprocessor is a serial program! \n");

  // Read preprocessor command-line arguements recorded in the .h5 file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", 
      H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  cmd_h5r -> read_string("/", "geo_file", geo_file);
  cmd_h5r -> read_string("/", "sur_file", sur_file);
  cmd_h5r -> read_string("/", "part_file", pre_part_file);
  elemType = cmd_h5r -> read_intScalar("/","elemType");
  dofNum = cmd_h5r -> read_intScalar("/","dofNum");
  probDim = cmd_h5r -> read_intScalar("/","probDim");
  in_ncommon = cmd_h5r -> read_intScalar("/","in_ncommon");

  delete cmd_h5r;
  
  // The user can specify the new mesh partition options
  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionBool("-isread_part", isread_part); 

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -METIS_isDualGraph: true \n";
  else cout<<" -METIS_isDualGraph: false \n";
  if(isread_part) cout<<" -isread_part: true \n";
  else cout<<" -isread_part: false \n";
  cout<<"----------------------------------\n";
  cout<<"geo_file: "<<geo_file<<endl;
  cout<<"sur_file: "<<sur_file<<endl;
  cout<<"pre_part_file: "<<pre_part_file<<endl;
  cout<<"probDim: "<<probDim<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  // Read the geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  // Check if the given geo file exist
  SYS_T::file_exist_check( geo_file.c_str() );

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);

  if(int(vecIEN.size()) != nElem * 4) SYS_T::print_fatal("Error: the IEN from geo_file does not match the given number of element. \n");

  if(int(ctrlPts.size()) != nFunc * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");

  std::cout<<"nElem: "<<nElem<<std::endl;
  std::cout<<"nFunc: "<<nFunc<<std::endl;

  IIEN * IEN = new IEN_Tetra_P1(nElem, vecIEN);

  VEC_T::clean(vecIEN);

  IMesh * mesh = new Mesh_Tet4(nFunc, nElem);
  mesh -> print_mesh_info();

  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "post_epart", "post_npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "post_epart", "post_npart" );
  else
  {
    cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<endl;
    exit(EXIT_FAILURE);
  }

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("post_node_mapping");

  cout<<"\n=== Start Partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_Tet( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );
    part->write(part_file.c_str());
    delete part;
  }

  // Clean memory
  cout<<endl<<"=== Clean memory. \n";
  delete mnindex; delete global_part; delete mesh; delete IEN;
  PetscFinalize();
  return 0;
}


// EOF
