// ==================================================================
// preprocess_tets.cpp
// ------------------------------------------------------------------
// This preprocess code is used for handling the 3D geometry described
// by tetrahedral elements. 
//
// Date: Dec. 18 2016
// Modified: May 22 2017
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_tet4.hpp"
#include "NBC_Partition_3D.hpp"
#include "NBC_Partition_3D_inflow.hpp"
#include "EBC_Partition_vtp.hpp"

int main( int argc, char * argv[] )
{
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  const int probDim = 3;
  const int dofNum = 4; // degree-of-freedom for the physical problem
  const int dofMat = 4; // degree-of-freedom in the matrix problem
  const int elemType = 501; // first order simplicial element

  std::string geo_file("./whole_vol.vtu");
  std::string part_file("part");

  int cpu_size = 1;
  int in_ncommon = 3;
  bool isDualGraph = true;

  PetscMPIInt rank, size;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionString("-geo_file", geo_file);

  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -part_file: "<<part_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<" -isDualGraph: true \n";
  std::cout<<"----------------------------------\n";
  std::cout<<" probDim: "<<probDim<<std::endl;
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;

  // ----- Write the input argument into a HDF5 file
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", 
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("cpu_size", cpu_size); 
  cmdh5w->write_intScalar("in_ncommon", in_ncommon); 
  cmdh5w->write_intScalar("dofNum", dofNum); 
  cmdh5w->write_intScalar("dofMat", dofMat); 
  cmdh5w->write_intScalar("probDim", probDim); 
  cmdh5w->write_intScalar("elemType", elemType); 
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; 
  H5Fclose(cmd_file_id);
  // ----- Finish writing

  // Read file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);

  std::cout<<"nElem: "<<nElem<<std::endl;
  std::cout<<"nFunc: "<<nFunc<<std::endl;

  IIEN * IEN = new IEN_Tetra_P1(nElem, vecIEN);

  VEC_T::clean( vecIEN );

  // check the tet mesh and print the quality of the mesh
  TET_T::tetmesh_check( ctrlPts, IEN, nElem, 3.5 ); 
  
  IMesh * mesh = new Mesh_Tet4(nFunc, nElem);
  mesh -> print_mesh_info();

  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  else
  {
    cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<endl;
    exit(EXIT_FAILURE);
  }

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");

  // ----------------------------------------------------------------
  // Setup boundary condition
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.resize( dofMat );

  std::vector<std::string> dir_list;
  dir_list.push_back("./inflow_vol.vtp");
  dir_list.push_back("./wall_vol.vtp");

  NBC_list[0] = new NodalBC_3D_vtp( nFunc );
  NBC_list[1] = new NodalBC_3D_vtp( dir_list, nFunc );
  NBC_list[2] = new NodalBC_3D_vtp( dir_list, nFunc );
  NBC_list[3] = new NodalBC_3D_vtp( dir_list, nFunc );

  std::vector<double> inflow_outward_normal;
  inflow_outward_normal.push_back(0.0);
  inflow_outward_normal.push_back(0.0);
  inflow_outward_normal.push_back(1.0);
  INodalBC * InFBC = new NodalBC_3D_inflow( "./inflow_vol.vtp", 
      "./wall_vol.vtp", nFunc, inflow_outward_normal );

  std::vector<std::string> ebclist;
  ebclist.clear();
  ElemBC * ebc = new ElemBC_3D_tet4( ebclist );

  // Correct the triangle's IEN so that we can get the outward normal 
  // direction easily. See the document for this reset function
  ebc -> resetTriIEN_outwardnormal( IEN );
  // ----------------------------------------------------------------

  const bool isPrintPartInfo = true;
  const int proc_size = cpu_size;

  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0;

  SYS_T::Timer * mytimer = new SYS_T::Timer();
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();
    IPart * part = new Part_Tet( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, dofMat, elemType,
        isPrintPartInfo );
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    part -> write( part_file.c_str() );

    part -> print_part_loadbalance_edgecut();

    // Partition Nodal BC
    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5(part_file.c_str());

    // Partition Inflow BC
    INBC_Partition * infpart = new NBC_Partition_3D_inflow(part, mnindex, InFBC);
    infpart->write_hdf5( part_file.c_str() );

    // Partition Elem BC
    IEBC_Partition * ebcpart = new EBC_Partition_vtp(part, mnindex, ebc);
    ebcpart -> write_hdf5(part_file.c_str());

    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete ebcpart; delete infpart;
  }

  VEC_T::write_int_h5("NumLocalNode","nln", list_nlocalnode);

  cout<<"\n===> Partition Quality: "<<endl;
  cout<<"The largest ghost / local node ratio is: ";
  cout<<*std::max_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;

  cout<<"The smallest ghost / local node ratio is: ";
  cout<<*std::min_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;

  cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<endl;

  cout<<"The maximum badnode number is: ";
  cout<<*std::max_element(&list_nbadnode[0], &list_nbadnode[cpu_size-1])<<endl;

  const int maxpart_nlocalnode = *std::max_element(&list_nlocalnode[0],
      &list_nlocalnode[cpu_size-1]);
  const int minpart_nlocalnode = *std::min_element(&list_nlocalnode[0],
      &list_nlocalnode[cpu_size-1]);

  cout<<"The maximum and minimum local node numbers are ";
  cout<<maxpart_nlocalnode<<"\t";
  cout<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of local node is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;

  // Free memory
  delete ebc; delete InFBC;
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;

  delete mnindex; delete global_part; delete mesh; delete IEN; 
  delete mytimer;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
