// ==================================================================
// prepost_tets.cpp
//
// This is the partitioning routine for parallel postprocessors.
//
// Date: July 29 2017
// ==================================================================
#include "HDF5_Reader.hpp"
#include "Tet_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet_FSI.hpp"

int main( int argc, char * argv[] )
{
  // clean the potentially pre-existing postpart h5 files
  SYS_T::execute("rm -rf ppart");
  SYS_T::execute("mkdir ppart");

  const std::string part_file("./ppart/postpart");
  
  int cpu_size = 1;
  bool isDualGraph = true;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor is a serial program! \n");

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string geo_file = cmd_h5r -> read_string("/", "geo_file");
  const int elemType = cmd_h5r -> read_intScalar("/","elemType");
  const int dofNum   = cmd_h5r -> read_intScalar("/","dofNum");
  const int dofMat   = cmd_h5r -> read_intScalar("/","dofMat");
  int in_ncommon = cmd_h5r -> read_intScalar("/","in_ncommon");

  delete cmd_h5r;
  H5Fclose(prepcmd_file);

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);
  
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -METIS_isDualGraph: true \n";
  else cout<<" -METIS_isDualGraph: false \n";
  cout<<"----------------------------------\n";
  cout<<"geo_file: "<<geo_file<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  int nFunc, nElem;
  std::vector<int> vecIEN, phy_tag;
  std::vector<double> ctrlPts;

  SYS_T::file_check( geo_file.c_str() );
  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN, phy_tag);

  if(int(vecIEN.size()) != nElem * 4) SYS_T::print_fatal("Error: the IEN from geo_file does not match the given number of element. \n");

  if(int(ctrlPts.size()) != nFunc * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");

  std::cout<<"nElem: "<<nElem<<std::endl;
  std::cout<<"nFunc: "<<nFunc<<std::endl;

  IIEN * IEN = new IEN_Tetra_P1(nElem, vecIEN);
  VEC_T::clean(vecIEN);

  std::vector<int> node_f, node_s; node_f.clear(); node_s.clear();

  for(int ee=0; ee<nElem; ++ee)
  {
    if( phy_tag[ee] == 0 )
    {
      for(int ii=0; ii<4; ++ii) node_f.push_back( IEN->get_IEN(ee, ii) );
    }
    else
    {
      for(int ii=0; ii<4; ++ii) node_s.push_back( IEN->get_IEN(ee, ii) );
    }
  }

  VEC_T::sort_unique_resize( node_f );
  VEC_T::sort_unique_resize( node_s );

  std::cout<<"Fluid domain number of nodes: "<<node_f.size()<<'\n';
  std::cout<<"Solid domain number of nodes: "<<node_s.size()<<'\n';

  IMesh * mesh = new Mesh_Tet4(nFunc, nElem);
  mesh -> print_info();

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
  SYS_T::Timer * mytimer = new SYS_T::Timer();
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    mytimer->Reset(); mytimer->Start();
    IPart * part = new Part_Tet_FSI( mesh, global_part, mnindex, IEN,
        ctrlPts, phy_tag, node_f, node_s, 
        proc_rank, proc_size, dofNum, dofMat, elemType, isPrintPartInfo );
    part->write(part_file.c_str());
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";
    delete part;
  }

  std::cout<<"=== Clean memory. \n";
  delete mytimer; delete mnindex; delete global_part; delete mesh; delete IEN;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
