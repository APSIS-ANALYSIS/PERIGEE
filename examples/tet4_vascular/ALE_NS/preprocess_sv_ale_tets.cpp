// ==================================================================
// preprocess_sv_ale_tets.cpp
// ------------------------------------------------------------------
// This preprocess code is used for handling the 3D geometry described
// by tetrahedral elements.
//
// Users should call the sv converter tool to convert the node and
// element indices to make them start from zero rather than one.
//
// Then the user should give the proper input file name and number of 
// outlets from the command-line argument to specify the volumetric 
// mesh and boundary mesh. 
//
// Date: Dec. 18 2016
// Modified: May 22 2017
// Modified: Mar. 17 2019
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_tet4_outflow.hpp"
#include "NBC_Partition_3D.hpp"
#include "NBC_Partition_3D_inflow.hpp"
#include "EBC_Partition_vtp_outflow.hpp"

int main( int argc, char * argv[] )
{
  // Call system routine to clean the previously generated hdf5 files
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  // Declare basic parameters
  const int dofNum = 7; // degree-of-freedom for the physical problem
  const int dofMat = 4; // degree-of-freedom in the matrix problem
  const int elemType = 501; // first order simplicial element

  // volumetric, wall, and cap files
  std::string geo_file("./whole_vol.vtu");
  std::string sur_file_in("./inflow_vol.vtp");
  std::string sur_file_wall("./wall_vol.vtp");
  std::string sur_file_out_base("./outflow_vol_");
  
  int num_outlet = 1; // Number of outlet faces

  // List of outlet faces files. Length is num_outlet, and
  // the name will be sur_file_out_basexx.vtp. e.g. outlet_00.vtp.
  std::vector< std::string > sur_file_out;
  
  const std::string part_file("part"); // name of partition file
 
  // mesh partition settings 
  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;

  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  // Get command-line argument
  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);
  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-sur_file_in", sur_file_in);
  SYS_T::GetOptionString("-sur_file_wall", sur_file_wall);
  SYS_T::GetOptionString("-sur_file_out_base", sur_file_out_base);

  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -num_outlet: "<<num_outlet<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -sur_file_in: "<<sur_file_in<<std::endl;
  std::cout<<" -sur_file_wall: "<<sur_file_wall<<std::endl;
  std::cout<<" -sur_file_out_base: "<<sur_file_out_base<<std::endl;
  std::cout<<" -part_file: "<<part_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<" -isDualGraph: true \n";
  std::cout<<"----------------------------------\n";
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;
  
  // Check if the file exists
  SYS_T::file_check(geo_file);
  std::cout<<geo_file<<" found. \n";

  SYS_T::file_check(sur_file_in);
  std::cout<<sur_file_in<<" found. \n";
  
  SYS_T::file_check(sur_file_wall);
  std::cout<<sur_file_wall<<" found. \n";
 
  // Generate sur_file_out file name and check its file existance 
  sur_file_out.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    std::ostringstream ss;
    ss<<sur_file_out_base;
    if( ii/10 == 0 ) ss<<"00";
    else if( ii/100 == 0 ) ss<<"0";

    ss<<ii<<".vtp";
    sur_file_out[ii] = ss.str(); // generate the outlet face file name
    SYS_T::file_check(sur_file_out[ii]);
    std::cout<<sur_file_out[ii]<<" found. \n";
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
  cmdh5w->write_string("sur_file_in", sur_file_in);
  cmdh5w->write_string("sur_file_out_base", sur_file_out_base);
  cmdh5w->write_string("sur_file_wall", sur_file_wall);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);
  // ----- Finish writing

  // Read the volumetric mesh file from a vtu file
  int nFunc, nElem; // number of nodes and elements
  std::vector<int> vecIEN; // IEN stored in a std::vector
  std::vector<double> ctrlPts; // x-y-z coordinates of the nodes

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);

  IIEN * IEN = new IEN_Tetra_P1(nElem, vecIEN);

  VEC_T::clean( vecIEN );

  // check and display the tet mesh quality
  TET_T::tetmesh_check(ctrlPts, IEN, nElem, 3.5);

  // Create the global mesh
  IMesh * mesh = new Mesh_Tet4(nFunc, nElem);
  mesh -> print_mesh_info();

  // Call METIS to partition the mesh
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

  // Generate the new numbering of the nodes based on the partition
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");

  // ----------------------------------------------------------------
  // Setup boundary condition
  // Nodal boundary conditions
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.resize( dofMat );

  // obtain the vtp files for the surfaces
  std::vector<std::string> dir_list;
  dir_list.push_back( sur_file_in );
  dir_list.push_back( sur_file_wall );

  NBC_list[0] = new NodalBC_3D_vtp( nFunc );
  NBC_list[1] = new NodalBC_3D_vtp( dir_list, nFunc );
  NBC_list[2] = new NodalBC_3D_vtp( dir_list, nFunc );
  NBC_list[3] = new NodalBC_3D_vtp( dir_list, nFunc );

  // Generate the nodal BC for the ALE mesh motion
  std::vector<INodalBC *> meshBC_list;
  meshBC_list.clear();
  meshBC_list.resize( 3 );

  std::vector<std::string> meshdir_vtp_list;
  meshdir_vtp_list.clear();
  meshdir_vtp_list.push_back( sur_file_in );
  meshdir_vtp_list.push_back( sur_file_wall );
  
  for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
    meshdir_vtp_list.push_back( sur_file_out[ii] );

  meshBC_list[0] = new NodalBC_3D_vtp( meshdir_vtp_list, nFunc );
  meshBC_list[1] = new NodalBC_3D_vtp( meshdir_vtp_list, nFunc );
  meshBC_list[2] = new NodalBC_3D_vtp( meshdir_vtp_list, nFunc );

  // Genereate the inflow bc info
  std::vector<double> inflow_outward_vec; // inflow surface outward normal vec
  TET_T::get_out_normal( sur_file_in, ctrlPts, IEN, inflow_outward_vec );
  INodalBC * InFBC = new NodalBC_3D_inflow( sur_file_in, sur_file_wall, 
      nFunc, inflow_outward_vec );

  // Elemental boundary conditions
  std::vector< std::vector<double> > outflow_outward_vec;
  outflow_outward_vec.resize( sur_file_out.size() );
  for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
    TET_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN, outflow_outward_vec[ii] );

  ElemBC * ebc = new ElemBC_3D_tet4_outflow( sur_file_out, outflow_outward_vec );

  // Correct the triangle's IEN so that we can get the outward normal 
  // direction easily. See the document for this reset function
  ebc -> resetTriIEN_outwardnormal( IEN );
  
  // Elemental BC for the ALE mesh
  std::vector<std::string> mesh_ebclist;
  mesh_ebclist.clear();
  ElemBC * mesh_ebc = new ElemBC_3D_tet4( mesh_ebclist );
  // ----------------------------------------------------------------
  
  // flag for printing partition statistics
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

    // Partition Nodal BC
    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5(part_file.c_str());

    // Partition Nodal Inflow BC
    INBC_Partition * infpart = new NBC_Partition_3D_inflow(part, mnindex, InFBC);
    infpart->write_hdf5( part_file.c_str() );

    // Partition Elem BC
    IEBC_Partition * ebcpart = new EBC_Partition_vtp_outflow(part, mnindex, ebc, NBC_list);

    ebcpart -> write_hdf5(part_file.c_str());

    // Partition Mesh Nodal BC
    INBC_Partition * mbcpart = new NBC_Partition_3D(part, mnindex, meshBC_list);
    mbcpart -> write_hdf5(part_file.c_str(), "/mesh_nbc");

    // Partition Mesh Elem BC
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

  // Free memory
  delete ebc; delete InFBC; delete mesh_ebc;
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) 
    delete *it_nbc;

  for(it_nbc=meshBC_list.begin(); it_nbc != meshBC_list.end(); ++it_nbc) 
    delete *it_nbc;

  delete mnindex; delete global_part; delete mesh; delete IEN; delete mytimer;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
