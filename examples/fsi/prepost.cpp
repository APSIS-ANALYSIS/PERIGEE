// ==================================================================
// prepost.cpp
//
// This is the partitioning routine for parallel postprocessors.
//
// Date: Jan 16 2022
// ==================================================================
#include "HDF5_Reader.hpp"
#include "Tet_Tools.hpp"
#include "Mesh_Tet.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM_FSI.hpp"

int main( int argc, char * argv[] )
{
  // clean the potentially pre-existing postpart h5 files
  if( SYS_T::directory_exist("apart") )
  {
    std::cout<<"Clean the folder ppart.\n";
    SYS_T::execute("rm -rf ppart");
  }

  SYS_T::execute("mkdir ppart");

  const std::string part_file_p("./ppart/postpart_p");
  const std::string part_file_v("./ppart/postpart_v");

  const int num_fields = 2; // Two fields : pressure + velocity/displacement
  const std::vector<int> dof_fields {1, 3}; // pressure 1 ; velocity/displacement 3

  int cpu_size = 1;
  bool isDualGraph = true;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  SYS_T::print_fatal_if(SYS_T::get_MPI_size() != 1, "ERROR: preprocessor is a serial program! \n");

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string geo_file = cmd_h5r -> read_string("/", "geo_file");
  const std::string sur_s_file_interior_wall = cmd_h5r -> read_string("/", "sur_s_file_interior_wall");
  const int elemType = cmd_h5r -> read_intScalar("/","elemType");
  int in_ncommon = cmd_h5r -> read_intScalar("/","in_ncommon");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);

  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -part_file_v: "<<part_file_v<<endl;
  cout<<" -part_file_p: "<<part_file_p<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -METIS_isDualGraph: true \n";
  else cout<<" -METIS_isDualGraph: false \n";
  cout<<"----------------------------------\n";
  cout<<"geo_file: "<<geo_file<<endl;
  cout<<"sur_s_file_interior_wall: "<<sur_s_file_interior_wall<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  // Read the geometry file for the whole FSI domain for the velocity /
  // displacement field
  int nFunc_v, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid( geo_file, nFunc_v, nElem, ctrlPts, vecIEN );
  
  const std::vector<int> phy_tag = VTK_T::read_int_CellData( geo_file, "Physics_tag" );

  // Generate IEN
  IIEN * IEN_v = new IEN_FEM( nElem, vecIEN );

  // --------------------------------------------------------------------------
  // The fluid-solid interface file will be read and the nodal index will be
  // mapped to a new value by the following rule. The ii-th node in the
  // interface wall node will be assgiend to nFunc_v + ii.
  // Read the F-S interface vtp file
  const std::vector<int> wall_node_id = VTK_T::read_int_PointData( sur_s_file_interior_wall, "GlobalNodeID" );

  const int nFunc_interface = static_cast<int>( wall_node_id.size() );
  const int nFunc_p = nFunc_v + nFunc_interface;

  // We will generate a new IEN array for the pressure variable by updating the
  // IEN for the solid element. If the solid element has node on the fluid-solid
  // interface, it will be mapped to the new index, that is nFunc + ii.
  std::vector<int> vecIEN_p ( vecIEN );

  for(int ee=0; ee<nElem; ++ee)
  {
    if( phy_tag[ee] == 1 )
    {
      // In solid element, loop over its IEN and correct if the node is on the
      // interface
      for(int ii=0; ii<4; ++ii)
      {
        const int pos = VEC_T::get_pos( wall_node_id, vecIEN_p[ee*4 +ii] );
        if( pos >=0 ) vecIEN_p[ee*4+ii] = nFunc_v + pos;
      }
    }
  }

  IIEN * IEN_p = new IEN_FEM( nElem, vecIEN_p );

  VEC_T::clean( vecIEN ); VEC_T::clean( vecIEN_p );
  // --------------------------------------------------------------------------

  // Generate the list of nodes for fluid and solid
  std::vector<int> v_node_f, v_node_s; v_node_f.clear(); v_node_s.clear();

  for(int ee=0; ee<nElem; ++ee)
  {
    if( phy_tag[ee] == 0 )
    {
      for(int ii=0; ii<4; ++ii) v_node_f.push_back( IEN_v->get_IEN(ee, ii) );
    }
    else
    {
      for(int ii=0; ii<4; ++ii) v_node_s.push_back( IEN_v->get_IEN(ee, ii) );
    }
  }

  VEC_T::sort_unique_resize( v_node_f ); VEC_T::sort_unique_resize( v_node_s );

  std::vector<int> p_node_f, p_node_s; p_node_f.clear(); p_node_s.clear();

  for(int ee=0; ee<nElem; ++ee)
  {
    if( phy_tag[ee] == 0 )
    {
      for(int ii=0; ii<4; ++ii) p_node_f.push_back( IEN_p->get_IEN(ee, ii) );
    }
    else
    {
      for(int ii=0; ii<4; ++ii) p_node_s.push_back( IEN_p->get_IEN(ee, ii) );
    }
  }

  VEC_T::sort_unique_resize( p_node_f ); VEC_T::sort_unique_resize( p_node_s );

  // Generate the mesh for kinematics
  IMesh * mesh_v = new Mesh_Tet(nFunc_v, nElem, 1);

  // Generate the mesh for pressure (discontinuous over interface)
  IMesh * mesh_p = new Mesh_Tet(nFunc_p, nElem, 1);

  std::vector<IMesh const *> mlist;
  mlist.push_back(mesh_p); mlist.push_back(mesh_v);

  mlist[0]->print_info();
  mlist[1]->print_info();

  std::cout<<"Fluid domain: "<<v_node_f.size()<<" nodes.\n";
  std::cout<<"Solid domain: "<<v_node_s.size()<<" nodes.\n";
  std::cout<<"Fluid-Solid interface: "<<nFunc_interface<<" nodes.\n";

  std::vector<IIEN const *> ienlist;
  ienlist.push_back(IEN_p); ienlist.push_back(IEN_v);

  // Partition the mesh
  IGlobal_Part * global_part = nullptr;
  if(cpu_size > 1)
  {
    global_part = new Global_Part_METIS( num_fields, cpu_size, in_ncommon, isDualGraph,
        mlist, ienlist, "post_epart", "post_npart" );
  }
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( num_fields, mlist, "post_epart", "post_npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  // Re-ordering nodal indices
  Map_Node_Index * mnindex_p = new Map_Node_Index(global_part, cpu_size, nFunc_p, 0);
  Map_Node_Index * mnindex_v = new Map_Node_Index(global_part, cpu_size, nFunc_v, 1);

  mnindex_p -> write_hdf5("post_node_mapping_p");
  mnindex_v -> write_hdf5("post_node_mapping_v");

  // Generate a list of local node number
  std::vector<int> list_nn_v(cpu_size), list_nn_p(cpu_size);
  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    // list stores the number of velo/pres nodes in each cpu
    list_nn_p[proc_rank] = 0; list_nn_v[proc_rank] = 0;
    for(int nn=0; nn<mesh_p -> get_nFunc(); ++nn)
    {
      if(global_part->get_npart(nn,0) == proc_rank) list_nn_p[proc_rank] += 1;
    }

    for(int nn=0; nn<mesh_v -> get_nFunc(); ++nn)
    {
      if(global_part->get_npart(nn,1) == proc_rank) list_nn_v[proc_rank] += 1;
    }
  }

  // Now generate the mappings from the gird pt idx to the matrix row idx
  // This is needed because we will have a matrix that has a special structure
  // due to the use of mix fem.
  std::vector<int> start_idx_v(cpu_size), start_idx_p(cpu_size);
  start_idx_v[0] = 0;
  start_idx_p[0] = 3 * list_nn_v[0];
  for(int ii = 1; ii < cpu_size; ++ii )
  {
    start_idx_v[ii] = start_idx_v[ii-1] + list_nn_v[ii-1]*3 + list_nn_p[ii-1];
    start_idx_p[ii] = start_idx_v[ii  ] + list_nn_v[ii  ]*3;
  }

  // mapper maps from the new grid point index to the matrix row index
  std::vector< std::vector<int> > mapper_p, mapper_v;
  mapper_p.resize(1); mapper_v.resize(3);

  for(int ii=0; ii<cpu_size; ++ii)
  {
    for(int jj=0; jj<list_nn_v[ii]; ++jj)
    {
      mapper_v[0].push_back( start_idx_v[ii] + jj * 3     );
      mapper_v[1].push_back( start_idx_v[ii] + jj * 3 + 1 );
      mapper_v[2].push_back( start_idx_v[ii] + jj * 3 + 2 );
    }

    for(int jj=0; jj<list_nn_p[ii]; ++jj)
      mapper_p[0].push_back( start_idx_p[ii] + jj );
  }

  // Start partition
  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();
    
    IPart * part_p = new Part_FEM_FSI( mesh_p, global_part, mnindex_p, IEN_p,
        ctrlPts, phy_tag, p_node_f, p_node_s,
        proc_rank, cpu_size, elemType, 0, dof_fields[0], start_idx_p[proc_rank], false );

    part_p -> write( part_file_p );
    delete part_p;

    IPart * part_v = new Part_FEM_FSI( mesh_v, global_part, mnindex_v, IEN_v,
        ctrlPts, phy_tag, v_node_f, v_node_s,
        proc_rank, cpu_size, elemType, 1, dof_fields[1], start_idx_v[proc_rank], true );

    part_v -> write( part_file_v );
    delete part_v;

    mytimer -> Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";
  }

  // Clean up Memory
  delete mnindex_p; delete mnindex_v; delete mesh_p; delete mesh_v;
  delete IEN_p; delete IEN_v; delete mytimer; delete global_part;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
