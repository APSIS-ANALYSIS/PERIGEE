// ==================================================================
// preprocess_fsi_beam.cpp
// ------------------------------------------------------------------
// This is a preprocessor code for the beam in flow FSI benchmark  
// in fully unstructured mesh.
//
// Date Created: Oct. 12 2017
// Author: Ju Liu
// ==================================================================
#include "Sys_Tools.hpp"
#include "Math_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet_FSI.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_vtu.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_tet4.hpp"
#include "NBC_Partition_3D.hpp"
#include "EBC_Partition_vtp.hpp"

int main( int argc, char * argv[] )
{
  // Clean up the h5 file from previous runs
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  const int probDim = 3;
  const int dofNum = 7;
  const int dofMat = 4;
  const int elemType = 501;

  std::string geo_file("./whole_vol.vtu");
  SYS_T::file_check(geo_file);

  std::string geo_fluid_file("./fluid.vtu");
  SYS_T::file_check(geo_fluid_file);

  std::string geo_solid_file("./solid.vtu");
  SYS_T::file_check(geo_solid_file);

  std::string part_file("part");

  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;

  PetscMPIInt rank, size;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  
  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<"----------------------------------\n";
  std::cout<<" geo_file: "<<geo_file<<std::endl;
  std::cout<<" geo_fluid_file: "<<geo_fluid_file<<std::endl;
  std::cout<<" geo_solid_file: "<<geo_solid_file<<std::endl;
  std::cout<<" probDim: "<<probDim<<std::endl;
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<" part_file: "<<part_file<<std::endl;
  std::cout<<" isDualGraph: true \n";
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;

  // Record info in h5 file
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
  cmdh5w->write_string("geo_fluid_file", geo_fluid_file);
  cmdh5w->write_string("geo_solid_file", geo_solid_file);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);

  // Read file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;
  std::vector<int> phy_tag;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN, phy_tag);

  for(unsigned int ii=0; ii<phy_tag.size(); ++ii)
  {
    if(phy_tag[ii] != 0 && phy_tag[ii] != 1)
      SYS_T::print_fatal("Error: FSI problem, the physical tag for element should be 0 for fluid domain or 1 for solid domain. \n");
  }

  std::cout<<"nElem: "<<nElem<<std::endl;
  std::cout<<"nFunc: "<<nFunc<<std::endl;

  // Generate IEN
  IIEN * IEN = new IEN_Tetra_P1(nElem, vecIEN);
  VEC_T::clean( vecIEN );

  // Generate the list of nodes for fluid and solid
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

  std::cout<<"\nFluid domain number of nodes: "<<node_f.size()<<'\n';
  std::cout<<"Solid domain number of nodes: "<<node_s.size()<<'\n';

  // Check the mesh
  TET_T::tetmesh_check( ctrlPts, IEN, nElem );

  // Generate the mesh
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

  // Setup the boundary conditions for the implicit soler
  cout<<"Boundary condition for the implicit solver: \n";
  std::vector<INodalBC *> NBC_list; NBC_list.clear();
  NBC_list.resize( dofMat );

  std::vector<std::string> dir_x_list; dir_x_list.clear();
  std::vector<std::string> dir_y_list; dir_y_list.clear();
  std::vector<std::string> dir_z_list; dir_z_list.clear();

  dir_x_list.push_back("flef_fluid.vtp");
  dir_x_list.push_back("cube_wall_fluid.vtp");
  dir_x_list.push_back("slef_solid.vtp");

  dir_y_list.push_back("flef_fluid.vtp");
  dir_y_list.push_back("ftop_fluid.vtp");
  dir_y_list.push_back("fbot_fluid.vtp");
  dir_y_list.push_back("cube_wall_fluid.vtp");
  dir_y_list.push_back("slef_solid.vtp");

  dir_z_list.push_back("flef_fluid.vtp");
  dir_z_list.push_back("ffro_fluid.vtp");
  dir_z_list.push_back("fbac_fluid.vtp");
  dir_z_list.push_back("cube_wall_fluid.vtp");
  dir_z_list.push_back("slef_solid.vtp");
  dir_z_list.push_back("sfro_solid.vtp");
  dir_z_list.push_back("sbac_solid.vtp");

  NBC_list[0] = new NodalBC_3D_vtp( nFunc );
  NBC_list[1] = new NodalBC_3D_vtp( dir_x_list, nFunc );
  NBC_list[2] = new NodalBC_3D_vtp( dir_y_list, nFunc );
  NBC_list[3] = new NodalBC_3D_vtp( dir_z_list, nFunc );

  //NBC_list[0] = new NodalBC_3D_vtu( geo_solid_file, nFunc );
  //NBC_list[1] = new NodalBC_3D_vtu( geo_solid_file, dir_x_list, nFunc );
  //NBC_list[2] = new NodalBC_3D_vtu( geo_solid_file, dir_y_list, nFunc );
  //NBC_list[3] = new NodalBC_3D_vtu( geo_solid_file, dir_z_list, nFunc );

  // Generate the mesh bc info
  cout<<"Boundary condition for the mesh motion: \n";
  std::vector<INodalBC *> meshBC_list;
  meshBC_list.clear(); meshBC_list.resize( 3 );

  std::vector<std::string> meshdir_vtp_list_x; meshdir_vtp_list_x.clear();
  meshdir_vtp_list_x.push_back("flef_fluid.vtp");
  meshdir_vtp_list_x.push_back("frig_fluid.vtp");
  meshdir_vtp_list_x.push_back("ftop_fluid.vtp");
  meshdir_vtp_list_x.push_back("fbot_fluid.vtp");
  meshdir_vtp_list_x.push_back("cube_wall_fluid.vtp");

  std::vector<std::string> meshdir_vtp_list_y; meshdir_vtp_list_y.clear();
  meshdir_vtp_list_y.push_back("ftop_fluid.vtp");
  meshdir_vtp_list_y.push_back("fbot_fluid.vtp");
  meshdir_vtp_list_y.push_back("flef_fluid.vtp");
  meshdir_vtp_list_y.push_back("frig_fluid.vtp");
  meshdir_vtp_list_y.push_back("cube_wall_fluid.vtp");

  std::vector<std::string> meshdir_vtp_list_z; meshdir_vtp_list_z.clear();
  meshdir_vtp_list_z.push_back("ftop_fluid.vtp");
  meshdir_vtp_list_z.push_back("fbot_fluid.vtp");
  meshdir_vtp_list_z.push_back("flef_fluid.vtp");
  meshdir_vtp_list_z.push_back("frig_fluid.vtp");
  meshdir_vtp_list_z.push_back("ffro_fluid.vtp");
  meshdir_vtp_list_z.push_back("fbac_fluid.vtp");
  meshdir_vtp_list_z.push_back("cube_wall_fluid.vtp");

  meshBC_list[0] = new NodalBC_3D_vtu( geo_solid_file, meshdir_vtp_list_x, nFunc );
  meshBC_list[1] = new NodalBC_3D_vtu( geo_solid_file, meshdir_vtp_list_y, nFunc );
  meshBC_list[2] = new NodalBC_3D_vtu( geo_solid_file, meshdir_vtp_list_z, nFunc );

  // Elemental BC
  cout<<"Elem boundary for the implicit solver: \n";
  std::vector<std::string> ebclist; ebclist.clear();
  ElemBC * ebc = new ElemBC_3D_tet4( ebclist );
  ebc -> resetTriIEN_outwardnormal( IEN );

  // Mesh EBC
  cout<<"Elem boundary for the mesh motion: \n";
  ebclist.clear();
  ElemBC * mesh_ebc = new ElemBC_3D_tet4( ebclist );

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

    IPart * part = new Part_Tet_FSI( mesh, global_part, mnindex, IEN,
        ctrlPts, phy_tag, node_f, node_s, proc_rank, proc_size, 
        dofNum, dofMat, elemType, isPrintPartInfo );

    mytimer -> Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    part -> write( part_file.c_str() );

    part -> print_part_loadbalance_edgecut();

    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5(part_file.c_str());

    INBC_Partition * mbcpart = new NBC_Partition_3D(part, mnindex, meshBC_list);
    mbcpart -> write_hdf5(part_file.c_str(), "/mesh_nbc");

    IEBC_Partition * ebcpart = new EBC_Partition_vtp(part, mnindex, ebc);
    ebcpart -> write_hdf5(part_file.c_str());

    IEBC_Partition * mebcpart = new EBC_Partition_vtp(part, mnindex, mesh_ebc);
    mebcpart-> write_hdf5(part_file.c_str(), "/mesh_ebc");

    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete mbcpart;
    delete ebcpart; delete mebcpart;
  }

  VEC_T::write_int_h5("NumLocalNode","nln", list_nlocalnode);

  // Print partition quality
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
  delete mesh_ebc; delete ebc;
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;
  for(it_nbc=meshBC_list.begin(); it_nbc != meshBC_list.end(); ++it_nbc) delete *it_nbc;

  delete mnindex; delete global_part; delete mesh; delete IEN; delete mytimer;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
