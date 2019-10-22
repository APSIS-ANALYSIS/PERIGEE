// prepost_tet_h5.cpp
// Dec. 8 2017
#include "HDF5_Reader.hpp"
#include "Mesh_FEM.hpp"
#include "IEN_Gmsh.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"

using std::cout;
using std::endl;

int main( int argc, char * argv[] )
{
  // clean the previously generated h5 files
  int sysret = system("rm -rf postpart_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  std::string geo_file, pre_part_file;
  int dofNum, dofMat, elemType, in_ncommon, probDim;
  int degree;

  const std::string part_file("postpart");
  int cpu_size = 1;
  bool isDualGraph = true;
  bool isread_part = true;
 
  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::print_fatal_if(size!=1,"ERROR: prepost is a serial program!\n");

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5",
      H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  cmd_h5r -> read_string("/", "geo_file", geo_file);
  cmd_h5r -> read_string("/", "part_file", pre_part_file);
  elemType = cmd_h5r -> read_intScalar("/","elemType");
  degree = cmd_h5r -> read_intScalar("/", "degree");
  dofNum = cmd_h5r -> read_intScalar("/","dofNum");
  dofMat = cmd_h5r -> read_intScalar("/","dofMat");
  probDim = cmd_h5r -> read_intScalar("/", "probDim");
  in_ncommon = cmd_h5r -> read_intScalar("/","in_ncommon");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);
  SYS_T::GetOptionBool("-isread_part", isread_part);

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -METIS_isDualGraph: true \n";
  else cout<<" -METIS_isDualGraph: false \n";
  if(isread_part) cout<<" -isread_part: true \n";
  else cout<<" -isread_part: false \n";
  cout<<"----------------------------------\n";
  cout<<"geo_file: "<<geo_file<<endl;
  cout<<"pre_part_file: "<<pre_part_file<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;
  
  SYS_T::file_exist_check( geo_file.c_str() );
  int nFunc, nElem, nLocBas;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;
  
  hid_t gfile_t = H5Fopen( geo_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * g_h5r = new HDF5_Reader( gfile_t );
  
  nFunc = g_h5r -> read_intScalar("/", "num_node");
  nElem = g_h5r -> read_intScalar("/", "num_cell");
  nLocBas = g_h5r -> read_intScalar("/", "nLocBas");
  g_h5r -> read_intVector("/", "IEN", vecIEN);
  g_h5r -> read_doubleVector("/", "node", ctrlPts);
  
  delete g_h5r; H5Fclose( gfile_t );

  IIEN * IEN = new IEN_Gmsh(nElem, nLocBas, vecIEN);
  VEC_T::clean( vecIEN );

  IMesh * mesh = new Mesh_FEM(nFunc, nElem, nLocBas, degree);
  mesh -> print_mesh_info();

  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "post_epart", "post_npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "post_epart", "post_npart" );
  else
  {
    std::cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<std::endl;
    exit(EXIT_FAILURE);
  }

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, nFunc);
  mnindex->write_hdf5("post_node_mapping");

  cout<<"\n=== Start Partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_FEM( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, dofMat,
        probDim, elemType, isPrintPartInfo );
    part->write(part_file.c_str());
    delete part;
  }

  // clean memory
  cout<<"=== Clean memory. \n";
  delete mnindex; delete global_part; delete mesh; delete IEN;
  PetscFinalize();
  return EXIT_SUCCESS;
}


// EOF
