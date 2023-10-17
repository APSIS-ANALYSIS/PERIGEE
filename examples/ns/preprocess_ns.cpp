// ==================================================================
// preprocess_ns.cpp
//
// This is a preprocessor code for handling Navier-Stokes equations 
// discretized by tetradedral elements.
//
// Date Created: Jan 01 2020
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "NBC_Partition.hpp"
#include "NBC_Partition_inflow.hpp"
#include "EBC_Partition_outflow.hpp"

int main( int argc, char * argv[] )
{
// Set number of threads and  print info of OpenMP
#ifdef _OPENMP
  omp_set_num_threads( omp_get_num_procs() );
  SYS_T::print_omp_info();
#endif

  // Clean the potentially pre-existing hdf5 files in the job folder
  SYS_T::execute("rm -rf part_p*.h5");
  SYS_T::execute("rm -rf preprocessor_cmd.h5");

  // Define basic problem settins
  constexpr int dofNum = 4; // degree-of-freedom for the physical problem
  constexpr int dofMat = 4; // degree-of-freedom in the matrix problem
  const std::string part_file("part");
  
  // Element options: 501 linear tets, 502 quadratic tets
  int elemType = 501;
  int num_outlet = 1;
  int num_inlet = 1;
  
  // Default names for input geometry files
  std::string geo_file("./whole_vol.vtu");
  std::string sur_file_in_base("./inflow_vol_");
  std::string sur_file_wall("./wall_vol.vtp");
  std::string sur_file_out_base("./outflow_vol_");

  // Mesh partition setting
  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./options_pre.yml");

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor needs to be run in serial.\n");

  if (is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // Get the command line arguments
  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionInt("-num_inlet", num_inlet);
  SYS_T::GetOptionInt("-num_outlet", num_outlet);
  SYS_T::GetOptionInt("-elem_type", elemType);
  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-sur_file_in_base", sur_file_in_base);
  SYS_T::GetOptionString("-sur_file_wall", sur_file_wall);
  SYS_T::GetOptionString("-sur_file_out_base", sur_file_out_base);

  if( elemType != 501 && elemType !=502 ) SYS_T::print_fatal("ERROR: unknown element type %d.\n", elemType);

  // Print the command line arguments
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -elem_type: "<<elemType<<endl;
  cout<<" -num_outlet: "<<num_outlet<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -sur_file_in_base: "<<sur_file_in_base<<endl;
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

  // If it is quadratic mesh, the mesh file will be all in vtu format
  if(elemType == 502)
  {
    sur_file_wall.erase( sur_file_wall.end()-4, sur_file_wall.end() );
    sur_file_wall += ".vtu";
  }

  SYS_T::file_check(sur_file_wall); cout<<sur_file_wall<<" found. \n";

  // Generate the inlet file names and check existance
  std::vector< std::string > sur_file_in;
  sur_file_in.resize( num_inlet );

  for(int ii=0; ii<num_inlet; ++ii)
  {  
    if(elemType == 501 )
      sur_file_in[ii] = SYS_T::gen_capfile_name( sur_file_in_base, ii, ".vtp" );   
    else
      sur_file_in[ii] = SYS_T::gen_capfile_name( sur_file_in_base, ii, ".vtu" ); 
  
    SYS_T::file_check(sur_file_in[ii]);
    cout<<sur_file_in[ii]<<" found. \n";
  }

  // Generate the outlet file names and check existance
  std::vector< std::string > sur_file_out;
  sur_file_out.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    if(elemType == 501 )
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtp" ); 
    else
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtu" ); 
 
    SYS_T::file_check(sur_file_out[ii]);
    cout<<sur_file_out[ii]<<" found. \n";
  }

  // Record the problem setting into a HDF5 file: preprocessor_cmd.h5
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("num_inlet", num_inlet);
  cmdh5w->write_intScalar("num_outlet", num_outlet);
  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("sur_file_in_base", sur_file_in_base);
  cmdh5w->write_string("sur_file_out_base", sur_file_out_base);
  cmdh5w->write_string("sur_file_wall", sur_file_wall);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);

  // Read the volumetric mesh file from the vtu file: geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;
  
  VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);
  
  IIEN * IEN = new IEN_FEM(nElem, vecIEN);
  VEC_T::clean( vecIEN ); // clean the vector
  
  IMesh * mesh = nullptr;

  switch( elemType )
  {
    case 501:
      mesh = new Mesh_Tet(nFunc, nElem, 1);
      break;
    case 502:
      mesh = new Mesh_Tet(nFunc, nElem, 2);
      break;
    default:
      SYS_T::print_fatal("Error: elemType %d is not supported.\n", elemType);
      break;
  }

  SYS_T::print_fatal_if( IEN->get_nLocBas() != mesh->get_nLocBas(), "Error: the nLocBas from the Mesh %d and the IEN %d classes do not match. \n", mesh->get_nLocBas(), IEN->get_nLocBas());

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

  // Setup Nodal i.e. Dirichlet type Boundary Conditions
  std::vector<INodalBC *> NBC_list( dofMat, nullptr );

  std::vector<std::string> dir_list {};
  for(int ii=0; ii<num_inlet; ++ii)
    dir_list.push_back( sur_file_in[ii] );
 
  dir_list.push_back( sur_file_wall );

  NBC_list[0] = new NodalBC( nFunc );
  NBC_list[1] = new NodalBC( dir_list, nFunc );
  NBC_list[2] = new NodalBC( dir_list, nFunc );
  NBC_list[3] = new NodalBC( dir_list, nFunc );

  // Inflow BC info
  std::vector< Vector_3 > inlet_outvec( sur_file_in.size() );
  for(unsigned int ii=0; ii<sur_file_in.size(); ++ii)
    inlet_outvec[ii] = TET_T::get_out_normal( sur_file_in[ii], ctrlPts, IEN );

  INodalBC * InFBC = new NodalBC_3D_inflow( sur_file_in, sur_file_wall,
      nFunc, inlet_outvec, elemType );

  InFBC -> resetSurIEN_outwardnormal( IEN ); // assign outward orientation for triangles


  // Setup Elemental Boundary Conditions
  // Obtain the outward normal vector
  std::vector< Vector_3 > outlet_outvec( sur_file_out.size() );
  for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
    outlet_outvec[ii] = TET_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN );

  ElemBC * ebc = new ElemBC_3D_outflow( sur_file_out, outlet_outvec, elemType );

  ebc -> resetSurIEN_outwardnormal( IEN ); // reset IEN for outward normal calculations
 
  // Start partition the mesh for each cpu_rank 

  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0; // total number of ghost nodes

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();
    IPart * part = new Part_FEM( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, cpu_size, dofNum, dofMat, elemType );
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file );

    part -> print_part_loadbalance_edgecut();
    
    // Partition Nodal BC and write to h5 file
    NBC_Partition * nbcpart = new NBC_Partition(part, mnindex, NBC_list);
    
    nbcpart -> write_hdf5( part_file );

    // Partition Nodal Inflow BC and write to h5 file
    NBC_Partition_inflow * infpart = new NBC_Partition_inflow(part, mnindex, InFBC);
    
    infpart->write_hdf5( part_file );
    
    // Partition Elemental BC and write to h5 file
    EBC_Partition * ebcpart = new EBC_Partition_outflow(part, mnindex, ebc, NBC_list);

    ebcpart -> write_hdf5( part_file );

    // Collect partition statistics
    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete infpart; delete ebcpart; 
  }

  cout<<"\n===> Mesh Partition Quality: "<<endl;
  cout<<"The largest ghost / local node ratio is: "<<VEC_T::max(list_ratio_g2l)<<endl;
  cout<<"The smallest ghost / local node ratio is: "<<VEC_T::min(list_ratio_g2l)<<endl;
  cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<endl;
  cout<<"The maximum badnode number is: "<<VEC_T::max(list_nbadnode)<<endl;

  const int maxpart_nlocalnode = VEC_T::max(list_nlocalnode); 
  const int minpart_nlocalnode = VEC_T::min(list_nlocalnode);

  cout<<"The maximum and minimum local node numbers are ";
  cout<<maxpart_nlocalnode<<"\t"<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of local node is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;

  // Finalize the code and exit
  for(auto it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;

  delete InFBC; delete ebc; delete mytimer;
  delete mnindex; delete global_part; delete mesh; delete IEN;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
