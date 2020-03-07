// ==================================================================
// preprocess_sv_fsi_tets.cpp
// ------------------------------------------------------------------
// This preprocess code is used for handling the 3D geometry described
// by tetrahdral elements and is tailored for the patient-specific
// vascular geometry generated from SimVascular.
//
// Users should call the sv_fsi_converter to convert the node and
// element indices before calling this routine.
//
// The users are also responsible for providing the proper file names
// for this routine to handle.
//
// Author: Ju Liu
// Date: April 3rd 2019
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet_FSI.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_vtu.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_tet4_outflow.hpp"
#include "NBC_Partition_3D.hpp"
#include "NBC_Partition_3D_inflow.hpp"
#include "EBC_Partition_vtp_outflow.hpp"

int main( int argc, char * argv[] )
{
  // Remove previously existing hdf5 files
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  // Define basic settings
  const int dofNum = 7; // degree-of-freedom for the physical problem
  const int dofMat = 4; // degree-of-freedom in the matrix problem
  const int elemType = 501; // first order simplicial element

  // Input files
  std::string geo_file("./whole_vol.vtu");

  std::string geo_f_file("./lumen_vol.vtu");
  std::string geo_s_file("./tissue_vol.vtu");

  std::string sur_f_file_in("./lumen_inlet_vol.vtp");
  std::string sur_f_file_wall("./lumen_wall_vol.vtp");
  std::string sur_f_file_out_base("./lumen_outlet_vol_");

  std::string sur_s_file_in("./tissue_inlet_vol.vtp");
  std::string sur_s_file_wall("./tissue_wall_vol.vtp");
  std::string sur_s_file_out_base("./tissue_outlet_vol_");

  int num_outlet = 1;

  const std::string part_file("part");

  int cpu_size = 1;
  int in_ncommon = 2;
  const bool isDualGraph = true;

  PetscMPIInt size;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);
  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-geo_f_file", geo_f_file);
  SYS_T::GetOptionString("-geo_s_file", geo_s_file);
  SYS_T::GetOptionString("-sur_f_file_in", sur_f_file_in);
  SYS_T::GetOptionString("-sur_f_file_wall", sur_f_file_wall);
  SYS_T::GetOptionString("-sur_f_file_out_base", sur_f_file_out_base);
  SYS_T::GetOptionString("-sur_s_file_in", sur_s_file_in);
  SYS_T::GetOptionString("-sur_s_file_wall", sur_s_file_wall);
  SYS_T::GetOptionString("-sur_s_file_out_base", sur_s_file_out_base);

  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -num_outlet: "<<num_outlet<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -geo_f_file: "<<geo_f_file<<std::endl;
  std::cout<<" -geo_s_file: "<<geo_s_file<<std::endl;
  std::cout<<" -sur_f_file_in: "<<sur_f_file_in<<std::endl;
  std::cout<<" -sur_f_file_wall: "<<sur_f_file_wall<<std::endl;
  std::cout<<" -sur_f_file_out_base: "<<sur_f_file_out_base<<std::endl;
  std::cout<<" -sur_s_file_in: "<<sur_s_file_in<<std::endl;
  std::cout<<" -sur_s_file_wall: "<<sur_s_file_wall<<std::endl;
  std::cout<<" -sur_s_file_out_base: "<<sur_s_file_out_base<<std::endl;
  std::cout<<" -part_file: "<<part_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<" -isDualGraph: true \n";
  std::cout<<"----------------------------------\n";
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;

  // Check if the geometrical file exist on disk
  SYS_T::file_check(geo_file); std::cout<<geo_file<<" found. \n";

  SYS_T::file_check(geo_f_file); std::cout<<geo_f_file<<" found. \n";

  SYS_T::file_check(geo_s_file); std::cout<<geo_s_file<<" found. \n";

  SYS_T::file_check(sur_f_file_in); std::cout<<sur_f_file_in<<" found. \n";

  SYS_T::file_check(sur_s_file_in); std::cout<<sur_s_file_in<<" found. \n";

  SYS_T::file_check(sur_f_file_wall); std::cout<<sur_f_file_wall<<" found. \n";

  SYS_T::file_check(sur_s_file_wall); std::cout<<sur_s_file_wall<<" found. \n";

  std::vector< std::string > sur_f_file_out;
  std::vector< std::string > sur_s_file_out;

  sur_f_file_out.resize( num_outlet );
  sur_s_file_out.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    std::ostringstream sf, ss;
    sf<<sur_f_file_out_base;
    ss<<sur_s_file_out_base;
    if( ii/10 == 0 )
    {
      sf << "00";
      ss << "00";
    }
    else if( ii/100 == 0 )
    {
      sf << "0";
      ss << "0";
    }

    sf<<ii<<".vtp";
    ss<<ii<<".vtp";

    sur_f_file_out[ii] = sf.str();
    sur_s_file_out[ii] = ss.str();

    SYS_T::file_check( sur_f_file_out[ii] );
    std::cout<<sur_f_file_out[ii]<<" found. \n";
    SYS_T::file_check( sur_s_file_out[ii] );
    std::cout<<sur_s_file_out[ii]<<" found. \n";
  } 

  // ----- Write the input argument into a HDF5 file
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5",
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("num_outlet", num_outlet);
  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("geo_f_file", geo_f_file);
  cmdh5w->write_string("geo_s_file", geo_s_file);
  cmdh5w->write_string("sur_f_file_in", sur_f_file_in);
  cmdh5w->write_string("sur_f_file_out_base", sur_f_file_out_base);
  cmdh5w->write_string("sur_f_file_wall", sur_f_file_wall);
  cmdh5w->write_string("sur_s_file_in", sur_s_file_in);
  cmdh5w->write_string("sur_s_file_out_base", sur_s_file_out_base);
  cmdh5w->write_string("sur_s_file_wall", sur_s_file_wall);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);
  // ----- Finish writing

  // Read the geometry file for the whole FSI domain
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<int> phy_tag;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN, phy_tag);

  for(unsigned int ii=0; ii<phy_tag.size(); ++ii)
  {
    if(phy_tag[ii] != 0 && phy_tag[ii] != 1) SYS_T::print_fatal("Error: FSI problem, the physical tag for element should be 0 (fluid domain) or 1 (solid domain).\n");
  }

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

  std::cout<<'\n'<<"Fluid domain number of nodes: "<<node_f.size()<<'\n';
  std::cout<<"Solid domain number of nodes: "<<node_s.size()<<'\n';
  
  // Check the mesh
  TET_T::tetmesh_check(ctrlPts, IEN, nElem, 3.5);

  // Generate the mesh
  IMesh * mesh = new Mesh_Tet4(nFunc, nElem);
  mesh -> print_mesh_info();

  // Partition
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

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");

  // ----------------------------------------------------------------
  // Setup boundary conditions
  // Physics NBC
  std::cout<<"Boundary condition for the implicit solver: \n";
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.resize( dofMat );

  std::vector<std::string> dir_list;
  dir_list.push_back( sur_f_file_in );
  dir_list.push_back( sur_s_file_in );
  for(int ii=0; ii<num_outlet; ++ii) dir_list.push_back( sur_s_file_out[ii] );

  NBC_list[0] = new NodalBC_3D_vtp( nFunc );
  NBC_list[1] = new NodalBC_3D_vtp( dir_list, nFunc );
  NBC_list[2] = new NodalBC_3D_vtp( dir_list, nFunc );
  NBC_list[3] = new NodalBC_3D_vtp( dir_list, nFunc );

  // Mesh NBC
  std::cout<<"Boundary condition for the mesh motion: \n";
  std::vector<INodalBC *> meshBC_list; meshBC_list.clear();
  meshBC_list.resize( 3 );
  std::vector<std::string> meshdir_vtp_list; meshdir_vtp_list.clear();
  meshdir_vtp_list.push_back( sur_f_file_in );
  for(int ii=0; ii<num_outlet; ++ii) 
    meshdir_vtp_list.push_back( sur_f_file_out[ii] );

  meshBC_list[0] = new NodalBC_3D_vtu( geo_s_file, meshdir_vtp_list, nFunc );
  meshBC_list[1] = new NodalBC_3D_vtu( geo_s_file, meshdir_vtp_list, nFunc );
  meshBC_list[2] = new NodalBC_3D_vtu( geo_s_file, meshdir_vtp_list, nFunc );

  // Generate inflow bc info
  std::vector<double> inflow_outward_vec; // inflow surface outward normal vec
  TET_T::get_out_normal( sur_f_file_in, ctrlPts, IEN, inflow_outward_vec );
  INodalBC * InFBC = new NodalBC_3D_inflow( sur_f_file_in, sur_f_file_wall,
      nFunc, inflow_outward_vec ); 

  // Elemental BC
  cout<<"Elem boundary for the implicit solver: \n";
  std::vector< std::vector<double> > outflow_outward_vec;
  outflow_outward_vec.resize( num_outlet );
  for(int ii=0; ii<num_outlet; ++ii)
    TET_T::get_out_normal( sur_f_file_out[ii], ctrlPts, IEN, outflow_outward_vec[ii] );

  ElemBC * ebc = new ElemBC_3D_tet4_outflow( sur_f_file_out, outflow_outward_vec );

  ebc -> resetTriIEN_outwardnormal( IEN );

  // Mesh EBC
  std::vector<std::string> mesh_ebclist;
  mesh_ebclist.clear();
  ElemBC * mesh_ebc = new ElemBC_3D_tet4( mesh_ebclist );
  // ----------------------------------------------------------------

  const bool isPrintPartInfo = true;

  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0;

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();

    IPart * part = new Part_Tet_FSI( mesh, global_part, mnindex, IEN,
        ctrlPts, phy_tag, node_f, node_s, 
        proc_rank, cpu_size, dofNum, dofMat, elemType, isPrintPartInfo );

    mytimer -> Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    part -> write( part_file.c_str() );
    part -> print_part_loadbalance_edgecut();

    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5(part_file.c_str()); 

    INBC_Partition * mbcpart = new NBC_Partition_3D(part, mnindex, meshBC_list);
    mbcpart -> write_hdf5(part_file.c_str(), "/mesh_nbc"); 

    INBC_Partition * infpart = new NBC_Partition_3D_inflow(part, mnindex, InFBC);
    infpart->write_hdf5( part_file.c_str() );

    IEBC_Partition * ebcpart = new EBC_Partition_vtp_outflow(part, mnindex, ebc, NBC_list);
    ebcpart -> write_hdf5(part_file.c_str());

    IEBC_Partition * mebcpart = new EBC_Partition_vtp(part, mnindex, mesh_ebc);
    mebcpart-> write_hdf5(part_file.c_str(), "/mesh_ebc");

    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete infpart; delete ebcpart;
    delete mbcpart; delete mebcpart;
  }

  // Print mesh partition statistics
  cout<<"\n===> Partition Quality: "<<endl;
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
  cout<<maxpart_nlocalnode<<'\t'<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of nlocalnode is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;

  // Clean memory
  for(auto it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc)
    delete *it_nbc;

  for(auto it_nbc=meshBC_list.begin(); it_nbc != meshBC_list.end(); ++it_nbc)
    delete *it_nbc;

  delete ebc; delete InFBC; delete mesh_ebc; delete mnindex; 
  delete global_part; delete mesh; delete IEN; delete mytimer;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
