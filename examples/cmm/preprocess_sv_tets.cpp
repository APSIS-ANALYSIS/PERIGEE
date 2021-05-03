// ==================================================================
// preprocess_sv_tets.cpp
//
// This is a preprocessor code for handling Navier-Stokes equations 
// discretized by tetradedral elements.
// 
// This preprocessor is similiar to that of the NS solver, please
// refer to the initialization of Elem_3D_tet_wall class for the
// specification of the wall properties, which will be utilized
// for FSI type simulations. Notice that the USERS may need to set
// proper file names for constructing Elem_3D_tet_wall and recompile
// the code.
//
// Date Created: Jan 01 2020
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "Mesh_Tet10.hpp"
#include "IEN_Tetra_P1.hpp"
#include "IEN_Tetra_P2.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_vtu.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "NodalBC_3D_ring.hpp"
#include "ElemBC_3D_tet_outflow.hpp"
#include "ElemBC_3D_tet_wall.hpp"
#include "NBC_Partition_3D_inflow.hpp"
#include "NBC_Partition_3D_ring.hpp"
#include "EBC_Partition_vtp_outflow.hpp"
#include "EBC_Partition_vtp_wall.hpp"

int main( int argc, char * argv[] )
{
  // Clean the potentially pre-existing hdf5 files in the job folder
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  // Define basic problem settins
  const int dofNum = 4; // degree-of-freedom for the physical problem
  const int dofMat = 4; // degree-of-freedom in the matrix problem
  const std::string part_file("part");
  
  // Element options: 501 linear tets, 502 quadratic tets
  int elemType = 501;
  int num_outlet = 1;
  
  // Default names for input geometry files
  std::string geo_file("./whole_vol.vtu");
  std::string sur_file_in("./inflow_vol.vtp");
  std::string sur_file_wall("./wall_vol.vtp");
  std::string sur_file_out_base("./outflow_vol_");

  // Mesh partition setting
  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor needs to be run in serial.\n");

  // Get the command line arguments
  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);
  SYS_T::GetOptionInt("-elem_type", elemType);
  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-sur_file_in", sur_file_in);
  SYS_T::GetOptionString("-sur_file_wall", sur_file_wall);
  SYS_T::GetOptionString("-sur_file_out_base", sur_file_out_base);

  if( elemType != 501 && elemType !=502 ) SYS_T::print_fatal("ERROR: unknown element type %d.\n", elemType);

  // Print the command line arguments
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -elem_type: "<<elemType<<endl;
  cout<<" -num_outlet: "<<num_outlet<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -sur_file_in: "<<sur_file_in<<endl;
  cout<<" -sur_file_wall: "<<sur_file_wall<<endl;
  cout<<" -sur_file_out_base: "<<sur_file_out_base<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  cout<<" -isDualGraph: true \n";
  cout<<"---- Problem definition ----\n";
  cout<<" dofNum: "<<dofNum<<endl;
  cout<<" dofMat: "<<dofMat<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // Check if the vtu geometry files exist on disk
  SYS_T::file_check(geo_file); cout<<geo_file<<" found. \n";

  // If quadratic, all mesh files will be in vtu format
  if(elemType == 502)
  {
    sur_file_in.erase( sur_file_in.end()-4, sur_file_in.end() );
    sur_file_in += ".vtu";
    sur_file_wall.erase( sur_file_wall.end()-4, sur_file_wall.end() );
    sur_file_wall += ".vtu";
  }

  SYS_T::file_check(sur_file_in); cout<<sur_file_in<<" found. \n";

  SYS_T::file_check(sur_file_wall); cout<<sur_file_wall<<" found. \n";

  // Generate the outlet file names and check existence
  std::vector< std::string > sur_file_out;
  sur_file_out.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    std::ostringstream ss;
    ss<<sur_file_out_base;
    if( ii/10 == 0 ) ss<<"00";
    else if( ii/100 == 0 ) ss<<"0";

    if(elemType == 501 ) ss<<ii<<".vtp";
    else ss<<ii<<".vtu";
 
    sur_file_out[ii] = ss.str(); // generate the outlet face file name
    
    SYS_T::file_check(sur_file_out[ii]);
    cout<<sur_file_out[ii]<<" found. \n";
  }

  // Record the problem setting into a HDF5 file: preprocessor_cmd.h5
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("num_outlet", num_outlet);
  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("sur_file_in", sur_file_in);
  cmdh5w->write_string("sur_file_out_base", sur_file_out_base);
  cmdh5w->write_string("sur_file_wall", sur_file_wall);
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

  // Inflow BC info
  Vector_3 inlet_outvec;
  TET_T::get_out_normal( sur_file_in, ctrlPts, IEN, inlet_outvec );
  
  INodalBC * InFBC = new NodalBC_3D_inflow( sur_file_in, sur_file_wall,
      nFunc, inlet_outvec, elemType );

  // Set up Elemental i.e. Neumann type Boundary Conditions
  // Obtain the outward normal vector
  std::vector< std::vector<double> > outflow_outward_vec;
  outflow_outward_vec.resize( sur_file_out.size() );
  for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
    TET_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN, outflow_outward_vec[ii] );

  std::vector< Vector_3 > outlet_outvec;
  outlet_outvec.resize( sur_file_out.size() );
  for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
    TET_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN, outlet_outvec[ii] );

  ElemBC * ebc = new ElemBC_3D_tet_outflow( sur_file_out, outlet_outvec, elemType );

  ebc -> resetTriIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations

  // Set up ring BC
  INodalBC * ring_bc = new NodalBC_3D_ring( sur_file_in, inlet_outvec,
       sur_file_wall, sur_file_out, outlet_outvec, nFunc, elemType );

  // Set up Nodal i.e. Dirichlet type Boundary Conditions. For CMM with prescribed inflow,
  // this includes all inlet interior nodes. To enable in-plane motion of the inlet & outlet
  // ring nodes, these ring nodes are included only for the dominant component of the
  // corresponding cap's unit normal.
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear(); NBC_list.resize( dofMat );

  if(elemType == 501)
  {
    NBC_list[0] = new NodalBC_3D_vtp( nFunc );
    NBC_list[1] = new NodalBC_3D_vtp( sur_file_in, inlet_outvec, sur_file_wall, sur_file_out, outflow_outward_vec, 1, 0, nFunc );
    NBC_list[2] = new NodalBC_3D_vtp( sur_file_in, inlet_outvec, sur_file_wall, sur_file_out, outflow_outward_vec, 1, 1, nFunc );
    NBC_list[3] = new NodalBC_3D_vtp( sur_file_in, inlet_outvec, sur_file_wall, sur_file_out, outflow_outward_vec, 1, 2, nFunc );
  }
  else
  {
    NBC_list[0] = new NodalBC_3D_vtu( nFunc );
    NBC_list[1] = new NodalBC_3D_vtu( sur_file_in, inlet_outvec, sur_file_wall, sur_file_out, outflow_outward_vec, 1, 0, nFunc );
    NBC_list[2] = new NodalBC_3D_vtu( sur_file_in, inlet_outvec, sur_file_wall, sur_file_out, outflow_outward_vec, 1, 1, nFunc );
    NBC_list[3] = new NodalBC_3D_vtu( sur_file_in, inlet_outvec, sur_file_wall, sur_file_out, outflow_outward_vec, 1, 2, nFunc );
  }

  // ----------------------------------------------------------------
  // Wall mesh for CMM-type model is set as an elemental bc.
  // Set the wall region, its corresponding centerline, and the thickness-to-radius ratio
  const std::string walls_combined = sur_file_wall;
  const std::string centerlines_combined = "centerlines.vtp";
  const double thickness2radius_combined = 0.2;

  // For variable wall properties:
  // If constructing wall properties with multiple spatial distributions,
  // provide three additional vectors of equal length: 
  //     1. wallsList:            surface vtp's, each a subset of the entire wall
  //     2. centerlinesList:      corresponding centerline vtp's
  //     3. thickness2radiusList: corresponding ratios
  // The background wall properties will first be prescribed to the entire wall
  // using centerlines_combined and thickness2radius_combined. Wall properties
  // in wallsList will then be overwritten using the corresponding lists.
  // ----------------------------------------------------------------
  // std::vector<std::string> wallsList; wallsList.clear();
  // std::vector<std::string> centerlinesList; centerlinesList.clear();
  // std::vector<double> thickness2radiusList; thickness2radiusList.clear();

  // if(elemType == 501) wallsList.push_back( "wall_aorta.vtp" );
  // else wallsList.push_back( "wall_aorta.vtu" );

  // centerlinesList.push_back( "centerlines_aorta.vtp" );
  // thickness2radiusList.push_back( 0.2 );

  // // Initialized with default fluid density 1.065
  // ElemBC * wall_bc = new ElemBC_3D_tet_wall( walls_combined, centerlines_combined,
  //                                            thickness2radius_combined, wallsList,
  //                                            centerlinesList, thickness2radiusList, elemType );

  ElemBC * wall_bc = nullptr;

  if( SYS_T::file_exist(centerlines_combined) )
    wall_bc = new ElemBC_3D_tet_wall( walls_combined, centerlines_combined,
        thickness2radius_combined, elemType );
  else
    wall_bc = new ElemBC_3D_tet_wall( walls_combined, 0.1, 1.0e6, elemType );

  wall_bc -> resetTriIEN_outwardnormal( IEN );
  // ----------------------------------------------------------------

  // Start partition the mesh for each cpu_rank 
  const bool isPrintPartInfo = true;

  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0; // total number of ghost nodes

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();

    IPart * part = new Part_Tet( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, cpu_size, dofNum, dofMat, elemType,
        isPrintPartInfo );
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file.c_str() );

    part -> print_part_loadbalance_edgecut();

    // Partition Nodal BC and write to h5 file
    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5( part_file.c_str() );

    // Partition Nodal Inflow BC and write to h5 file
    INBC_Partition * infpart = new NBC_Partition_3D_inflow(part, mnindex, InFBC);
    infpart->write_hdf5( part_file.c_str() );

    // Partition Nodal Ring BC and write to h5 file
    INBC_Partition * ringpart = new NBC_Partition_3D_ring(part, mnindex, ring_bc);
    ringpart->write_hdf5( part_file.c_str() );

    // Partition Elemental BC and write to h5 file
    IEBC_Partition * ebcpart = new EBC_Partition_vtp_outflow(part, mnindex, ebc, NBC_list);
    ebcpart -> write_hdf5( part_file.c_str() );

    // Partition Elemental Wall BC and write it to h5 file
    IEBC_Partition * wbcpart = new EBC_Partition_vtp_wall(part, mnindex, wall_bc );
    wbcpart -> write_hdf5( part_file.c_str() );

    // Collect partition statistics
    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete infpart; delete ringpart; delete ebcpart; delete wbcpart; 
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

  delete InFBC; delete ring_bc; delete ebc; delete mytimer; delete wall_bc;
  delete mnindex; delete global_part; delete mesh; delete IEN;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
