// ==================================================================
// preprocess_tet_h5.cpp
// ------------------------------------------------------------------
// This preprocess code is used to handle 3D geometry discretized
// by arbitrary order tetrahedral elements, which is stored in a h5
// file.
//
// Date created: Dec. 6 2017
// Author: Ju Liu
// ==================================================================
#include "IEN_Gmsh.hpp"
#include "Mesh_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Map_Node_Index.hpp"
#include "NodalBC_3D.hpp"
#include "ElemBC_3D.hpp"
#include "Part_FEM.hpp"
#include "NBC_Partition_3D.hpp"
#include "EBC_Partition_FEM.hpp"

using std::cout;
using std::endl;

int main( int argc, char * argv[] )
{
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  const int probDim = 3;
  const int dofNum = 7;
  const int dofMat = 4;
  const std::string part_file("part");

  std::string geo_file("./Gmsh_vol.h5");
  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  SYS_T::print_fatal_if(size!= 1,"ERROR: preprocessor is a serial program!\n");

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionString("-geo_file", geo_file);

  std::cout<<"==== Command Line Arguments ===="<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<" -isDualGraph: true \n";

  // ----------------------------------------------------------------
  // ------------- Read the info from the geo_file ------------------
  SYS_T::file_check(geo_file);
  hid_t gfile_t = H5Fopen( geo_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * g_h5r = new HDF5_Reader( gfile_t );

  int nFunc, nElem, nFace, elemType, nLocBas, start_index_3d;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  nFunc = g_h5r -> read_intScalar("/", "num_node");
  nElem = g_h5r -> read_intScalar("/", "num_cell");
  nFace = g_h5r -> read_intScalar("/", "num_face");
  elemType = g_h5r -> read_intScalar("/", "ele_type");
  nLocBas = g_h5r -> read_intScalar("/", "nLocBas");
  start_index_3d = g_h5r -> read_intScalar("/", "start_index_3d");
  g_h5r -> read_intVector("/", "IEN", vecIEN);
  g_h5r -> read_doubleVector("/", "node", ctrlPts);
  // ----------------------------------------------------------------
  // Boundary surface
  std::vector<std::string> face_name; face_name.resize(nFace);
  std::vector<int> nElem_2D; nElem_2D.resize(nFace);
  std::vector<int> nFunc_2D; nFunc_2D.resize(nFace);
  std::vector<int> eType_2D; eType_2D.resize(nFace);
  std::vector<int> nLocBas_2D; nLocBas_2D.resize(nFace);
  std::vector<std::vector<int> > IEN_glo_2D; IEN_glo_2D.resize(nFace);
  std::vector<std::vector<int> > IEN_loc_2D; IEN_loc_2D.resize(nFace);
  std::vector<std::vector<int> > bcpt_2D; bcpt_2D.resize(nFace);
  std::vector<std::vector<int> > face2elem; face2elem.resize(nFace);
  std::vector<std::vector<double> > pt_coor_2D; pt_coor_2D.resize(nFace);

  for(int ii=0; ii<nFace; ++ii)
  {
    const std::string g1d_name = std::to_string(ii);
    g_h5r->read_string( g1d_name.c_str(), "name", face_name[ii]);
    nElem_2D[ii] = g_h5r -> read_intScalar(g1d_name.c_str(), "num_cell");
    nFunc_2D[ii] = g_h5r -> read_intScalar(g1d_name.c_str(), "num_node");
    eType_2D[ii] = g_h5r -> read_intScalar(g1d_name.c_str(), "ele_type");
    nLocBas_2D[ii] = g_h5r -> read_intScalar(g1d_name.c_str(), "nLocBas");
    g_h5r -> read_intVector(g1d_name.c_str(), "IEN_glo", IEN_glo_2D[ii]);
    g_h5r -> read_intVector(g1d_name.c_str(), "IEN_loc", IEN_loc_2D[ii]);
    g_h5r -> read_intVector(g1d_name.c_str(), "pt_idx", bcpt_2D[ii]);
    g_h5r -> read_intVector(g1d_name.c_str(), "face2elem", face2elem[ii]);
    g_h5r -> read_doubleVector(g1d_name.c_str(), "pt_coor", pt_coor_2D[ii]);
  }

  delete g_h5r; H5Fclose( gfile_t );
  // ----------------------------------------------------------------

  std::cout<<"==== Data from "<<geo_file<<" ===="<<std::endl;
  std::cout<<" nFunc: "<<nFunc<<'\n';
  std::cout<<" nElem: "<<nElem<<'\n';
  std::cout<<" nFace: "<<nFace<<'\n';
  std::cout<<" nLocBas: "<<nLocBas<<'\n';
  std::cout<<" start_index_3d: "<<start_index_3d<<'\n';
  std::cout<<" nFunc_2D: "; VEC_T::print(nFunc_2D);
  std::cout<<" nFunc_2D: "; VEC_T::print(nFunc_2D);
  std::cout<<" eType_2D: "; VEC_T::print(eType_2D);
  std::cout<<" nLocBas_2D: "; VEC_T::print(nLocBas_2D);
  std::cout<<" 2D Face Name: "; VEC_T::print(face_name);

  // Check format
  if(int(vecIEN.size()) != nElem * nLocBas) SYS_T::print_fatal("Error: the IEN from geo_file does not match the given number of element. \n");

  if(int(ctrlPts.size()) != nFunc * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");

  if( start_index_3d != 0 ) SYS_T::print_fatal("Error: the start index is non-zero. Check the mesh. \n");

  for(int ii=0; ii<nFace; ++ii)
    if(eType_2D[ii] != eType_2D[0]) SYS_T::print_fatal("Error: the face element type is not uniform. Check the mesh. \n");

  // Define the interpolation degree
  int degree = -1, my_elem_type = -1;
  switch( elemType )
  {
    case 4:
      degree = 1;
      my_elem_type = 531;
      break;
    case 11:
      degree = 2;
      my_elem_type = 532;
      break;
    default:
      SYS_T::print_fatal("Error: the element type is not considered in the analysis code. \n");
  }

  std::cout<<" degree: "<<degree<<'\n';
  std::cout<<" elemType (Gmsh): "<<elemType<<'\n';
  std::cout<<" elemType (Mine): "<<my_elem_type<<'\n';

  // ----------------------------------------------------------------  
  // Write arguments into a hdf5 file
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5",
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("num_face", nFace);
  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("probDim", probDim);
  cmdh5w->write_intScalar("elemType", my_elem_type); 
  cmdh5w->write_intScalar("degree", degree);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);
  // ----------------------------------------------------------------  
  
  // Generate IEN
  IIEN * IEN = new IEN_Gmsh(nElem, nLocBas, vecIEN);
  VEC_T::clean( vecIEN );

  IMesh * mesh = new Mesh_FEM(nFunc, nElem, nLocBas, degree);

  mesh -> print_info();

  // Call METIS for partitioning
  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  else
  {
    std::cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<std::endl;
    exit(EXIT_FAILURE);
  }

  // Generate nodal index mapping
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");

  // Nodal BC 
  std::vector<INodalBC *> NBC_list; NBC_list.clear();
  NBC_list.resize( dofMat );

  NBC_list[0] = new NodalBC_3D( nFunc );
  
  std::vector<int> vnbclist; vnbclist.clear();
  vnbclist.push_back(0);
  vnbclist.push_back(5);
  NBC_list[1] = new NodalBC_3D( vnbclist, bcpt_2D, nFunc );

  vnbclist.clear();
  vnbclist.push_back(0);
  vnbclist.push_back(2);
  NBC_list[2] = new NodalBC_3D( vnbclist, bcpt_2D, nFunc);

  vnbclist.clear();
  vnbclist.push_back(1);
  NBC_list[3] = new NodalBC_3D( vnbclist, bcpt_2D, nFunc );
  
  //std::vector<int> vnbclist; vnbclist.clear();
  //vnbclist.push_back(0); vnbclist.push_back(1); vnbclist.push_back(2);
  //NBC_list[3] = new NodalBC_3D( vnbclist, bcpt_2D, nFunc );

  // Element BC
  std::vector<int> ebclist; ebclist.clear();
  ebclist.push_back( 0 );

  ElemBC * ebc = new ElemBC_3D( ebclist, nElem_2D, nFunc_2D,
      nLocBas_2D, pt_coor_2D, IEN_loc_2D, bcpt_2D, face2elem, IEN, ctrlPts );

  // ----------------------------------------------------------------
  // Partition the mesh & BC
  const bool isPrintPartInfo = true;
  const int proc_size = cpu_size;
  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;
  int sum_nghostnode = 0;
  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    mytimer->Reset(); mytimer->Start();
    IPart * part = new Part_FEM( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, dofMat,
        probDim, my_elem_type, isPrintPartInfo );
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    part -> write( part_file.c_str() );

    part -> print_part_loadbalance_edgecut();

    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());
    sum_nghostnode += part->get_nghostnode();

    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5(part_file.c_str());

    IEBC_Partition * ebcpart = new EBC_Partition_FEM(part, mnindex, ebc);
    ebcpart -> write_hdf5(part_file.c_str());

    delete part; delete nbcpart; delete ebcpart;
  }
  // ----------------------------------------------------------------

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

  // clean memory
  delete mytimer; delete ebc;
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;
  delete IEN; delete mesh; delete global_part; delete mnindex;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
