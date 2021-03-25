#include "PDNSolution_Wall_Disp.hpp"

PDNSolution_Wall_Disp::PDNSolution_Wall_Disp( 
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const int &type, const bool &isprint )
: PDNSolution( pNode, 3 ), is_print( isprint )
{
  switch( type )
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: in PDNSolution_Wall_Disp, No such type of initial condition. \n");
      break;
  }
}

// ==== WOMERSLEY CHANGES BEGIN ====
PDNSolution_Wall_Disp::PDNSolution_Wall_Disp(
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const ALocal_EBC * const &ebc_wall_part,
    const double &rho,
    const int &type, const bool &isprint )
: PDNSolution( pNode, 3 ), is_print( isprint )
{

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 1:
      Init_womersley( agmi_ptr, pNode, ebc_wall_part, rho );
      break;
    case 2:
      Init_womersley_dot( agmi_ptr, pNode, ebc_wall_part, rho );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Wall_Disp : No such type of initial condition.\n");
      break;
  }
}
// ==== WOMERSLEY CHANGES END ====


PDNSolution_Wall_Disp::~PDNSolution_Wall_Disp()
{}


void PDNSolution_Wall_Disp::Init_zero( const APart_Node * const &pNode_ptr )
{
  int location[3];
  const double value[3] = {0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);

  GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
    SYS_T::commPrint("                       disp_y = 0.0 \n");
    SYS_T::commPrint("                       disp_z = 0.0 \n");
  }
}


// ==== WOMERSLEY CHANGES BEGIN ====
void PDNSolution_Wall_Disp::Init_womersley(
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const APart_Node * const &pNode_ptr,
    const ALocal_EBC * const &ebc_wall_part,
    const double &rho )
{
  int location[3];
  double value[3] = {0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_snode = ebc_wall_part->get_num_local_node(ebc_id);

  const double R     = 0.3;                                                  // pipe radius
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);

  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
  const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

  // g1 = 2 * besselJ(1, Lambda) / ( Lambda * besselJ(0, Lambda) )
  const std::complex<double> g1(0.403118398110473,      -0.319596413162681); 

  std::vector<int> analysis_old2new;
  analysis_old2new.resize( agmi_ptr->get_nFunc() );
  ReadNodeMapping("node_mapping.h5", "old_2_new", agmi_ptr->get_nFunc(), &analysis_old2new[0] );

  for(int ii=0; ii<num_snode; ++ii)
  {
    const int global_node = ebc_wall_part->get_local_global_node(ebc_id, ii);
    const int pos = analysis_old2new[global_node];

    location[0] = pos * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    const double x = ebc_wall_part->get_local_pt_xyz(ebc_id, 3*ii);
    const double y = ebc_wall_part->get_local_pt_xyz(ebc_id, 3*ii + 1);
    const double z = ebc_wall_part->get_local_pt_xyz(ebc_id, 3*ii + 2);

    // axial wall disp
    const double xi = std::real( i1 * B1 / (rho * c1 * omega) * (G1 - 1.0) * exp(-i1*omega*z/c1) );

    // radial wall disp
    const double eta = std::real( B1 * R / (2.0 * rho * c1 * c1) * (1.0 - G1 * g1) * exp(-i1*omega*z/c1) );

    // polar to cartesian transformation
    const double theta = std::atan2(y, x);
    
    value[0] = eta * std::cos(theta);
    value[1] = eta * std::sin(theta);
    value[2] = xi;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }
  
  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();
}


void PDNSolution_Wall_Disp::Init_womersley_dot(
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const APart_Node * const &pNode_ptr,
    const ALocal_EBC * const &ebc_wall_part,
    const double &rho )
{
  int location[3];
  double value[3] = {0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_snode = ebc_wall_part->get_num_local_node(ebc_id);

  const double R     = 0.3;                                                  // pipe radius
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);

  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
  const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

  // g1 = 2 * besselJ(1, Lambda) / ( Lambda * besselJ(0, Lambda) )
  const std::complex<double> g1(0.403118398110473,      -0.319596413162681); 

  std::vector<int> analysis_old2new;
  analysis_old2new.resize( agmi_ptr->get_nFunc() );
  ReadNodeMapping("node_mapping.h5", "old_2_new", agmi_ptr->get_nFunc(), &analysis_old2new[0] );

  for(int ii=0; ii<num_snode; ++ii)
  {
    const int global_node = ebc_wall_part->get_local_global_node(ebc_id, ii);
    const int pos = analysis_old2new[global_node];

    location[0] = pos * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    const double x = ebc_wall_part->get_local_pt_xyz(ebc_id, 3*ii);
    const double y = ebc_wall_part->get_local_pt_xyz(ebc_id, 3*ii + 1);
    const double z = ebc_wall_part->get_local_pt_xyz(ebc_id, 3*ii + 2);

    // axial dot wall disp
    const double dot_xi = std::real( -B1 / (rho * c1) * (G1 - 1.0) * exp(-i1*omega*z/c1) );

    // radial dot wall disp
    const double dot_eta = std::real( i1 * omega * B1 * R / (2.0 * rho * c1 * c1) * (1.0 - G1 * g1) * exp(-i1*omega*z/c1) );

    // polar to cartesian transformation
    const double theta = std::atan2(y, x);
    
    value[0] = dot_eta * std::cos(theta);
    value[1] = dot_eta * std::sin(theta);
    value[2] = dot_xi;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

}


// Read the node mappings
void PDNSolution_Wall_Disp::ReadNodeMapping(
    const char * const &node_mapping_file,
    const char * const &mapping_type, const int &node_size,
    int * const &nodemap )
{
  hid_t file_id = H5Fopen(node_mapping_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t data_id = H5Dopen(file_id, mapping_type, H5P_DEFAULT);

  hid_t data_space = H5Dget_space( data_id );
  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  if( data_rank != 1)
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: the node mapping file has wrong format. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  hsize_t * data_dims = new hsize_t [1];

  H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  hsize_t dSize = data_dims[0];

  if( int(dSize) != node_size )
  {
    PetscPrintf(PETSC_COMM_SELF, "Error: the allocated array has wrong size! \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space,
      H5P_DEFAULT, nodemap );

  delete [] data_dims;
  H5Sclose( mem_space );
  H5Sclose(data_space);
  H5Dclose(data_id);
  H5Fclose(file_id);
}

// EOF
