// ==================================================================
// prepost_ns_tets.cpp
//
// This is the partitioning routine for parallel postprocessors.
//
// Date: Jan. 24 2017
// ==================================================================
#include "HDF5_Reader.hpp"
#include "Tet_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "Mesh_Tet10.hpp"
#include "IEN_Tetra_P1.hpp"
#include "IEN_Tetra_P2.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet.hpp"

int main( int argc, char * argv[] )
{
  // Clean the existing part hdf5 files
  int sysret = system("rm -rf postpart_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "ERROR: system call failed. \n");

  std::string geo_file;
  int dofNum, dofMat, elemType, in_ncommon;

  const std::string part_file("postpart");
  int cpu_size = 1;
  bool isDualGraph = true;

  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1, "ERROR: preprocessor is a serial program! \n");

  // Read preprocessor command-line arguements recorded in the .h5 file
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  cmd_h5r -> read_string("/", "geo_file", geo_file);
  elemType = cmd_h5r -> read_intScalar("/","elemType");
  dofNum = cmd_h5r -> read_intScalar("/","dofNum");
  dofMat   = cmd_h5r -> read_intScalar("/","dofMat");
  in_ncommon = cmd_h5r -> read_intScalar("/","in_ncommon");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // The user can specify the new mesh partition options
  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -METIS_isDualGraph: true \n";
  else cout<<" -METIS_isDualGraph: false \n";
  cout<<"----------------------------------\n";
  cout<<"part_file: "<<part_file<<endl;
  cout<<"geo_file: "<<geo_file<<endl;
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
  
  if(int(ctrlPts.size()) != nFunc * 3) SYS_T::print_fatal("ERROR: the ctrlPts from geo_file does not match the given number of nodes. \n");

  IIEN * IEN = nullptr;
  IMesh * mesh = nullptr;

  if(elemType == 501)
  {
    SYS_T::print_fatal_if(vecIEN.size() / nElem != 4, "ERROR: the mesh connectivity array size does not match with the element type 501. \n");

    IEN = new IEN_Tetra_P1(nElem, vecIEN);
    mesh = new Mesh_Tet4(nFunc, nElem);
  }
  else if(elemType == 502)
  {
    SYS_T::print_fatal_if(vecIEN.size() / nElem != 10, "ERROR: the mesh connectivity array size does not match with the element type 502. \n");

    IEN = new IEN_Tetra_P2(nElem, vecIEN);
    mesh = new Mesh_Tet10(nFunc, nElem);
  }
  else SYS_T::print_fatal("ERROR: element type %d not supported.\n", elemType);

  mesh -> print_info();

  IGlobal_Part * global_part = nullptr;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "post_epart", "post_npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "post_epart", "post_npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("post_node_mapping");

  cout<<"\n=== Start Partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  SYS_T::Timer * mytimer = new SYS_T::Timer();
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    mytimer->Reset(); mytimer->Start();
    IPart * part = new Part_Tet( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, dofMat, elemType,
        isPrintPartInfo );
    part->write(part_file.c_str());
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";
    delete part;
  }

  // Clean memory
  cout<<"\n=== Clean memory. \n";
  delete mnindex; delete global_part; delete mesh; delete IEN; delete mytimer;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
