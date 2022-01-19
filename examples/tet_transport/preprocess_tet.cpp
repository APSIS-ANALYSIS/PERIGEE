// ============================================================================
// preprocess_tet.cpp
//
// ============================================================================
#include "Mesh_Tet4.hpp"
#include "Mesh_Tet10.hpp"
#include "IEN_Tetra_P1.hpp"
#include "IEN_Tetra_P2.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Tet_Tools.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "ElemBC_3D_tet.hpp"
#include "Part_Tet.hpp"
#include "NBC_Partition.hpp"
#include "EBC_Partition.hpp"

int main( int argc, char * argv[] )
{
  SYS_T::execute("rm -rf preprocessor_cmd.h5");
  SYS_T::execute("rm -rf apart");
  SYS_T::execute("mkdir apart");
  
  // Define the partition file name 
  const std::string part_file("./apart/part");

  // Element option: 501 linear tet element
  int elemType = 501;

  // Default names for input geometry files
  std::string geo_file("./whole_vol.vtu");

  // Mesh partition setting
  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;
  
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor needs to be run in serial.\n");
  
  // Get the command line arguments
  SYS_T::GetOptionInt(   "-cpu_size",          cpu_size);
  SYS_T::GetOptionInt(   "-in_ncommon",        in_ncommon);
  SYS_T::GetOptionInt(   "-elem_type",         elemType);
  SYS_T::GetOptionString("-geo_file",          geo_file);
  
  if( elemType != 501 && elemType !=502 ) SYS_T::print_fatal("ERROR: unknown element type %d.\n", elemType);
  
  // Print the command line arguments
  std::cout << "==== Command Line Arguments ====" << std::endl;
  std::cout << " -elem_type: "   << elemType      << std::endl;
  std::cout << " -geo_file: "         << geo_file          << std::endl;
  std::cout << " -part_file: "        << part_file         << std::endl;
  std::cout << " -cpu_size: "         << cpu_size          << std::endl;
  std::cout << " -in_ncommon: "       << in_ncommon        << std::endl;
  std::cout << " -isDualGraph: true \n";
  std::cout << "====  Command Line Arguments/ ===="<<std::endl;
  
  // Check if the vtu geometry files exist on disk
  SYS_T::file_check(geo_file); std::cout<<geo_file<<" found. \n";

  // Record the problem setting into a HDF5 file: preprocessor_cmd.h5
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);
  
  // Read the volumetric mesh file from the vtu file: geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);
  
  IIEN * IEN = nullptr;
  IMesh * mesh = nullptr;

  if(elemType == 501)
  {
    SYS_T::print_fatal_if(vecIEN.size() / nElem != 4, "Error: the mesh connectivity array size does not match with the element type 501. \n");

    IEN = new IEN_Tetra_P1(nElem, vecIEN);
    mesh = new Mesh_Tet4(nFunc, nElem);
  }
  else
  {
    SYS_T::print_fatal_if(vecIEN.size() / nElem != 10, "Error: the mesh connectivity array size does not match with the element type 502. \n");

    IEN = new IEN_Tetra_P2(nElem, vecIEN);
    mesh = new Mesh_Tet10(nFunc, nElem);
  }

  VEC_T::clean( vecIEN ); // clean the vector

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

  // Setup Nodal (Dirichlet type) boundary condition(s)
  std::vector<std::string> dir_list;
  dir_list.push_back("bot_vol.vtp");

  INodalBC * nbc = new NodalBC_3D_vtp( dir_list, nFunc );
  
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.push_back( nbc );

  // Setup Elemental (Neumann type) boundary condition(s)
  std::vector<std::string> neu_list;
  neu_list.push_back("top_vol.vtp");
  ElemBC * ebc = new ElemBC_3D_tet( neu_list, elemType );
  
  ebc -> resetTriIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations
  
  // Start partition the mesh for each cpu_rank
  SYS_T::Timer * mytimer = new SYS_T::Timer();

  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0; // total number of ghost nodes

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();
    
    IPart * part = new Part_Tet( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, cpu_size, 1, 1, elemType,
        true );

    part -> print_part_loadbalance_edgecut();

    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file.c_str() );

    // Partition Nodal BC and write to h5 file
    NBC_Partition * nbcpart = new NBC_Partition(part, mnindex, NBC_list);

    nbcpart -> write_hdf5( part_file.c_str() );

    // Partition Elemental BC and write to h5 file
    EBC_Partition * ebcpart = new EBC_Partition(part, mnindex, ebc);

    ebcpart -> write_hdf5( part_file.c_str() );

    // Collect partition statistics
    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();

    delete part; delete ebcpart; delete nbcpart;
  }
  
  cout<<"\n===> Mesh Partition Quality: "<<endl;
  cout<<"The largest ghost / local node ratio is: ";
  cout<<*std::max_element(list_ratio_g2l.begin(), list_ratio_g2l.end())<<endl;

  cout<<"The smallest ghost / local node ratio is: ";
  cout<<*std::min_element(list_ratio_g2l.begin(), list_ratio_g2l.end())<<endl;

  cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<endl;

  cout<<"The maximum badnode number is: ";
  cout<<*std::max_element(list_nbadnode.begin(), list_nbadnode.end())<<endl;

  const int maxpart_nlocalnode = *std::max_element(list_nlocalnode.begin(),
      list_nlocalnode.end());
  const int minpart_nlocalnode = *std::min_element(list_nlocalnode.begin(),
      list_nlocalnode.end());

  cout<<"The maximum and minimum local node numbers are ";
  cout<<maxpart_nlocalnode<<"\t";
  cout<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of local node is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;

  // Finalize the code and exit
  for(auto it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;

  delete mytimer;
  delete ebc; delete global_part; delete mnindex; delete IEN; delete mesh;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
