// ==================================================================
// ns_tet_driver.cpp
//
// Tetrahedral element based finite element code for 3D Navier-Stokes
// equations using Variational Multiscale Formulation and Generalized
// alpha time stepping.
//
// Author: Ju Liu, liujuy@gmail.com
// Date: Feb. 6 2020
// ==================================================================
#include "HDF5_Writer.hpp"
#include "AGlobal_Mesh_Info.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node_Rotated.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "ALocal_InflowBC.hpp"
#include "ALocal_RotatedBC.hpp"
#include "ALocal_Interface.hpp"
#include "Sliding_Interface_Tools.hpp"
#include "Matrix_Free_Tools.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "QuadPts_UserDefined_Triangle.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "CVFlowRate_Cosine2Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"
#include "PLocAssem_VMS_NS_GenAlpha.hpp"
#include "PLocAssem_VMS_NS_GenAlpha_WeakBC.hpp"
#include "PLocAssem_VMS_NS_GenAlpha_Interface.hpp"
#include "PGAssem_NS_FEM.hpp"
#include "PTime_NS_Solver.hpp"
#include "SI_rotation_info.hpp"

int main(int argc, char *argv[])
{
  // Coefficient for weak bc
  double C_bI = 4.0;

  // Number of quadrature points for tets and triangles
  // Suggested values: 5 / 4 for linear, 17 / 13 for quadratic
  int nqp_tet = 5, nqp_tri = 4;

  // Number of quadrature points for hexs and quadrangles
  // Suggested values: 2 / 2 for linear, 4 / 4 for quadratic
  int nqp_vol_1D = 2, nqp_sur_1D = 2;

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // fluid properties
  double fluid_density = 1.065;
  double fluid_mu = 3.5e-2;
  double c_tauc = 1.0; // scaling factor for tau_c, take 0.0, 0.125, or 1.0
  double c_ct = 4.0; // C_T parameter for defining tau_M

  // inflow file
  std::string inflow_file("inflow_fourier_series.txt");

  double inflow_thd_time = 0.01; // prescribed time for inflow to reach steadness
  double inflow_tgt_rate = 1.0; // prescribed flow rate at steady state

  // Turbulence intensity for the purtabation at inlets, 3% ==> 0.03
  double inflow_TI_perturbation = 0.0;

  // LPN file
  std::string lpn_file("lpn_rcr_input.txt");

  // back flow stabilization
  double bs_beta = 0.2;

  // generalized-alpha rho_inf
  double genA_rho_inf = 0.5;

  // part file location
  std::string part_file("part");

  // nonlinear solver parameters
  double nl_rtol = 1.0e-3; // convergence criterion relative tolerance
  double nl_atol = 1.0e-6; // convergence criterion absolute tolerance
  double nl_dtol = 10.0;   // divergence criterion
  int nl_maxits = 20;      // maximum number if nonlinear iterations
  int nl_refreq = 4;       // frequency of tangent matrix renewal
  int nl_threshold = 4;    // threshold of tangent matrix renewal

  // time stepping parameters
  double initial_time = 0.0; // time of the initial condition
  double initial_step = 0.1; // time step size
  int initial_index = 0;     // indiex of the initial condition
  double final_time = 1.0;   // final time
  std::string sol_bName("SOL_"); // base name of the solution file
  int ttan_renew_freq = 1;   // frequency of tangent matrix renewal
  int sol_record_freq = 1;   // frequency of recording the solution

  // Restart options
  bool is_restart = false;
  int restart_index = 0;     // restart solution time index
  double restart_time = 0.0; // restart time
  double restart_step = 1.0e-3; // restart simulation time step size
  std::string restart_name = "SOL_"; // restart solution base name

  // Angular velocity
  double angular_velo = 250 * MATH_T::PI; //(rad/s)
  double angular_thd_time = 1.5; // prescribed time for rotating part to reach angular velocity // 1.5 1.0

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./runscript.yml");

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_perigee_art();

  // ===== Yaml Arguments =====
  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  if (is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionReal("-angular_velo", angular_velo);
  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
  SYS_T::GetOptionInt("-nqp_vol_1d", nqp_vol_1D);
  SYS_T::GetOptionInt("-nqp_sur_1d", nqp_sur_1D);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionReal("-bs_beta", bs_beta);
  SYS_T::GetOptionReal("-rho_inf", genA_rho_inf);
  SYS_T::GetOptionReal("-fl_density", fluid_density);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionReal("-c_tauc", c_tauc);
  SYS_T::GetOptionReal("-c_ct", c_ct);
  SYS_T::GetOptionString("-inflow_file", inflow_file);
  SYS_T::GetOptionReal("-inflow_thd_time", inflow_thd_time);
  SYS_T::GetOptionReal("-inflow_tgt_rate", inflow_tgt_rate);
  SYS_T::GetOptionReal("-inflow_TI_perturbation", inflow_TI_perturbation);
  SYS_T::GetOptionReal("-angular_velo", angular_velo);
  SYS_T::GetOptionReal("-angular_thd_time", angular_thd_time);
  SYS_T::GetOptionString("-lpn_file", lpn_file);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionReal("-nl_rtol", nl_rtol);
  SYS_T::GetOptionReal("-nl_atol", nl_atol);
  SYS_T::GetOptionReal("-nl_dtol", nl_dtol);
  SYS_T::GetOptionInt("-nl_maxits", nl_maxits);
  SYS_T::GetOptionInt("-nl_refreq", nl_refreq);
  SYS_T::GetOptionInt("-nl_threshold", nl_threshold);
  SYS_T::GetOptionReal("-init_time", initial_time);
  SYS_T::GetOptionReal("-fina_time", final_time);
  SYS_T::GetOptionReal("-init_step", initial_step);
  SYS_T::GetOptionInt("-init_index", initial_index);
  SYS_T::GetOptionInt("-ttan_freq", ttan_renew_freq);
  SYS_T::GetOptionInt("-sol_rec_freq", sol_record_freq);
  SYS_T::GetOptionString("-sol_name", sol_bName);
  SYS_T::GetOptionBool("-is_restart", is_restart);
  SYS_T::GetOptionInt("-restart_index", restart_index);
  SYS_T::GetOptionReal("-restart_time", restart_time);
  SYS_T::GetOptionReal("-restart_step", restart_step);
  SYS_T::GetOptionString("-restart_name", restart_name);
  SYS_T::GetOptionReal("-C_bI", C_bI);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint("-angular_velo", angular_velo);
  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-nqp_vol_1d", nqp_vol_1D);
  SYS_T::cmdPrint("-nqp_sur_1d", nqp_sur_1D);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-bs_beta:", bs_beta);
  SYS_T::cmdPrint("-rho_inf:", genA_rho_inf);
  SYS_T::cmdPrint("-fl_density:", fluid_density);
  SYS_T::cmdPrint("-fl_mu:", fluid_mu);
  SYS_T::cmdPrint("-c_tauc:", c_tauc);
  SYS_T::cmdPrint("-c_ct:", c_ct);

  // if inflow file exists, print the file name
  // otherwise, print the parameter for linear2steady inflow setting
  if( SYS_T::file_exist( inflow_file ) )
    SYS_T::cmdPrint("-inflow_file:", inflow_file);
  else
  {
    SYS_T::cmdPrint("-inflow_thd_time:", inflow_thd_time);
    SYS_T::cmdPrint("-inflow_tgt_rate:", inflow_tgt_rate);
  }
  SYS_T::cmdPrint("-inflow_TI_perturbation:", inflow_TI_perturbation);

  SYS_T::cmdPrint("-lpn_file:", lpn_file);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-nl_rtol:", nl_rtol);
  SYS_T::cmdPrint("-nl_atol:", nl_atol);
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol);
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq);
  SYS_T::cmdPrint("-nl_threshold:", nl_threshold);
  SYS_T::cmdPrint("-init_time:", initial_time);
  SYS_T::cmdPrint("-init_step:", initial_step);
  SYS_T::cmdPrint("-init_index:", initial_index);
  SYS_T::cmdPrint("-fina_time:", final_time);
  SYS_T::cmdPrint("-ttan_freq:", ttan_renew_freq);
  SYS_T::cmdPrint("-sol_rec_freq:", sol_record_freq);
  SYS_T::cmdPrint("-sol_name:", sol_bName);
  if(is_restart)
  {
    SYS_T::commPrint("-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else SYS_T::commPrint("-is_restart: false \n");

  // ===== Record important solver options =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5",
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar("fl_density", fluid_density);
    cmdh5w->write_doubleScalar("fl_mu", fluid_mu);
    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_intScalar("sol_record_freq", sol_record_freq);
    cmdh5w->write_string("lpn_file", lpn_file);
    cmdh5w->write_doubleScalar("angular_velo", angular_velo);

    if( SYS_T::file_exist( inflow_file ) )
      cmdh5w->write_string("inflow_file", inflow_file);
    else
    {
      cmdh5w->write_doubleScalar("inflow_thd_time", inflow_thd_time );
      cmdh5w->write_doubleScalar("inflow_tgt_rate", inflow_tgt_rate );
    }
    cmdh5w->write_doubleScalar("inflow_TI_perturbation", inflow_TI_perturbation);
    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Data from Files =====
  // Read the info of rotation axis from h5file
  const std::string fName = SYS_T::gen_partfile_name( part_file, rank );
  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  const std::string gname("/rotation");
  const Vector_3 point_rotated( h5r->read_Vector_3( gname.c_str(), "point_rotated" ) );
  const Vector_3 angular_direction( h5r->read_Vector_3( gname.c_str(), "angular_direction" ) );
  delete h5r; H5Fclose( file_id );

  // Control points' xyz coordinates
  FEANode * fNode = new FEANode(part_file, rank);

  // Local sub-domain's IEN array
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // Global mesh info
  AGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info(part_file,rank);

  // Mesh partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  // Local sub-domain's element indices
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // Local sub-domain's nodal bc
  ALocal_NBC * locnbc = new ALocal_NBC(part_file, rank);

  // Local sub-domain's inflow bc
  ALocal_InflowBC * locinfnbc = new ALocal_InflowBC(part_file, rank);

  // Local sub-domain's rotated bc
  ALocal_RotatedBC * locrotnbc = new ALocal_RotatedBC(part_file, rank);

  // Local sub-domain's elemental bc
  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);

  // Local sub_domain's weak bc
  ALocal_WeakBC * locwbc = new ALocal_WeakBC(part_file, rank);
  locwbc -> print_info();

  // Interfaces info
  ALocal_Interface * locitf = new ALocal_Interface(part_file, rank);
  SYS_T::commPrint("Interfaces: %d\n", locitf->get_num_itf());
  locitf -> print_info();

  SI_T::SI_solution * SI_sol = new SI_T::SI_solution(part_file, rank);

  // Local sub-domain's nodal indices
  APart_Node * pNode = new APart_Node_Rotated(part_file, rank);

  SI_rotation_info * sir_info = new SI_rotation_info(angular_velo, angular_thd_time, point_rotated, angular_direction);

  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");

  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Inflow flow rate =====
  SYS_T::commPrint("===> Setup inflow flow rate. \n");

  ICVFlowRate * inflow_rate_ptr = nullptr;

  // If inflow file exist, load it
  // otherwise, call the linear incremental flow rate to reach a steady flow
  // if( SYS_T::file_exist( inflow_file ) )
  //   inflow_rate_ptr = new CVFlowRate_Unsteady( inflow_file.c_str() );
  // else
  inflow_rate_ptr = new CVFlowRate_Cosine2Steady( inflow_thd_time, inflow_TI_perturbation, inflow_file );

  inflow_rate_ptr->print_info();

  // ===== Finite Element Container & Quadrature rules =====
  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  FEAElement * anchor_elementv = nullptr;
  FEAElement * opposite_elementv = nullptr;

  SYS_T::commPrint("===> Build quadrature rules. \n");
  const int nqp_vol { (GMIptr->get_elemType() == FEType::Tet4 || GMIptr->get_elemType() == FEType::Tet10) ? nqp_tet : (nqp_vol_1D * nqp_vol_1D * nqp_vol_1D) };
  const int nqp_sur { (GMIptr->get_elemType() == FEType::Tet4 || GMIptr->get_elemType() == FEType::Tet10) ? nqp_tri : (nqp_sur_1D * nqp_sur_1D) };

  IQuadPts * quadv = nullptr;
  IQuadPts * quads = nullptr;
  IQuadPts * free_quad = nullptr;

  if( GMIptr->get_elemType() == FEType::Tet4 )
  {
    if( nqp_tet > 5 ) SYS_T::commPrint("Warning: the tet element is linear and you are using more than 5 quadrature points.\n");
    if( nqp_tri > 4 ) SYS_T::commPrint("Warning: the tri element is linear and you are using more than 4 quadrature points.\n");

    elementv = new FEAElement_Tet4( nqp_vol ); // elem type Tet4
    elements = new FEAElement_Triangle3_3D_der0( nqp_sur );
    anchor_elementv = new FEAElement_Tet4( nqp_sur );
    opposite_elementv = new FEAElement_Tet4( 1 );
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
    quads = new QuadPts_Gauss_Triangle( nqp_sur );
    free_quad = new QuadPts_UserDefined_Triangle();
  }
  else if( GMIptr->get_elemType() == FEType::Tet10 )
  {
    SYS_T::print_fatal_if( nqp_tet < 29, "Error: not enough quadrature points for tets.\n" );
    SYS_T::print_fatal_if( nqp_tri < 13, "Error: not enough quadrature points for triangles.\n" );

    elementv = new FEAElement_Tet10( nqp_vol ); // elem type Tet10
    elements = new FEAElement_Triangle6_3D_der0( nqp_sur );
    anchor_elementv = new FEAElement_Tet10( nqp_sur );
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
    quads = new QuadPts_Gauss_Triangle( nqp_sur );
  }
  else if( GMIptr->get_elemType() == FEType::Hex8 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 2, "Error: not enough quadrature points for hex.\n" );
    SYS_T::print_fatal_if( nqp_sur_1D < 1, "Error: not enough quadrature points for quad.\n" );

    elementv = new FEAElement_Hex8( nqp_vol ); // elem type Hex8
    elements = new FEAElement_Quad4_3D_der0( nqp_sur );
    anchor_elementv = new FEAElement_Hex8( nqp_sur );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else if( GMIptr->get_elemType() == FEType::Hex27 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 4, "Error: not enough quadrature points for hex.\n" );
    SYS_T::print_fatal_if( nqp_sur_1D < 3, "Error: not enough quadrature points for quad.\n" );

    elementv = new FEAElement_Hex27( nqp_vol ); // elem type Hex27
    elements = new FEAElement_Quad9_3D_der0( nqp_sur );
    anchor_elementv = new FEAElement_Hex27( nqp_sur );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  SI_T::SI_quad_point * SI_qp = new SI_T::SI_quad_point(locitf, nqp_sur);

  // ===== Generate a sparse matrix for the enforcement of essential BCs
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  pmat->gen_perm_bc(pNode, locnbc);

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local Assembly routine =====
  IPLocAssem * locAssem_ptr = nullptr;
  // if( locwbc->get_wall_model_type() == 0 )
  // {
  //   locAssem_ptr = new PLocAssem_VMS_NS_GenAlpha(
  //     tm_galpha_ptr, elementv->get_nLocBas(),
  //     quadv->get_num_quadPts(), elements->get_nLocBas(),
  //     fluid_density, fluid_mu, bs_beta, GMIptr->get_elemType(), point_rotated, angular_velo, c_ct, c_tauc );
  // }
  // else if( locwbc->get_wall_model_type() == 1 )
  // {
  //   locAssem_ptr = new PLocAssem_VMS_NS_GenAlpha_WeakBC(
  //     tm_galpha_ptr, elementv->get_nLocBas(),
  //     quadv->get_num_quadPts(), elements->get_nLocBas(),
  //     fluid_density, fluid_mu, bs_beta, GMIptr->get_elemType(), point_rotated, angular_velo, c_ct, c_tauc, C_bI );
  // }
  // else SYS_T::print_fatal("Error: Unknown wall model type.\n");

  locAssem_ptr = new PLocAssem_VMS_NS_GenAlpha_Interface(
    tm_galpha_ptr, elementv->get_nLocBas(),
    quadv->get_num_quadPts(), elements->get_nLocBas(),
    fluid_density, fluid_mu, bs_beta, GMIptr->get_elemType(), angular_velo, point_rotated, angular_direction, c_ct, c_tauc, C_bI );

  // ===== Initial condition =====
  PDNSolution * base = new PDNSolution_NS( pNode, fNode, locinfnbc, 1 );

  PDNSolution * sol = new PDNSolution_NS( pNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_NS( pNode, 0 );

  PDNSolution * disp_mesh = new PDNSolution_V( pNode );

  PDNSolution * velo_mesh = new PDNSolution_V( pNode );

  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    // Read sol file
    SYS_T::file_check(restart_name);
    sol->ReadBinary(restart_name);

    // generate the corresponding dot_sol file name
    std::string restart_dot_name = "dot_";
    restart_dot_name.append(restart_name);

    // generate the corresponding mdisp file name
    std::string restart_mdisp_name = "DISP_";
    restart_mdisp_name.append(std::to_string(900000000 + initial_index));

    // generate the corresponding mvelo file name
    std::string restart_mvelo_name = "MVELO_";
    restart_mvelo_name.append(std::to_string(900000000 + initial_index));

    // Read dot_sol file
    SYS_T::file_check(restart_dot_name);
    dot_sol->ReadBinary(restart_dot_name);

    // Read mdisp file
    SYS_T::file_check(restart_mdisp_name);
    disp_mesh->ReadBinary(restart_mdisp_name);

    // Read mvelo file
    SYS_T::file_check(restart_mvelo_name);
    velo_mesh->ReadBinary(restart_mvelo_name);

    SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
    SYS_T::commPrint("     restart_name: %s \n", restart_name.c_str());
    SYS_T::commPrint("     restart_dot_name: %s \n", restart_dot_name.c_str());
    SYS_T::commPrint("     restart_mdisp_name: %s \n", restart_mdisp_name.c_str());
    SYS_T::commPrint("     restart_mvelo_name: %s \n", restart_mvelo_name.c_str());
    SYS_T::commPrint("     restart_time: %e \n", restart_time);
    SYS_T::commPrint("     restart_index: %d \n", restart_index);
    SYS_T::commPrint("     restart_step: %e \n", restart_step);
  }

  SI_sol->update_node_sol(sol);
  SI_sol->update_node_mdisp(disp_mesh);
  SI_sol->update_node_mvelo(velo_mesh);
  
  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ===== LPN models =====
  IGenBC * gbc = nullptr;

  if( GENBC_T::get_genbc_file_type( lpn_file ) == 1  )
    gbc = new GenBC_Resistance( lpn_file );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 2  )
    gbc = new GenBC_RCR( lpn_file, 1000, initial_step );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 3  )
    gbc = new GenBC_Inductance( lpn_file );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 4  )
    gbc = new GenBC_Coronary( lpn_file, 1000, initial_step, initial_index );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 5  )
    gbc = new GenBC_Pressure( lpn_file, initial_time );
  else
    SYS_T::print_fatal( "Error: GenBC input file %s format cannot be recongnized.\n", lpn_file.c_str() );

  gbc -> print_info();

  // Make sure the gbc number of faces matches that of ALocal_EBC
  SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc->get_num_ebc(),
      "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Global assembly =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  SI_qp->search_all_opposite_point(anchor_elementv, opposite_elementv, elements, quads, free_quad, locitf, SI_sol);

  IPGAssem * gloAssem_ptr = new PGAssem_NS_FEM( locAssem_ptr, elements, anchor_elementv, opposite_elementv, quads, free_quad,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc, locitf, SI_sol, SI_qp, gbc, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr,
      elements, quads, locIEN, pNode, locnbc, locebc, gbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Initialize the shell matrix =====
  Mat shell_mat;
  PetscInt local_row_size, local_col_size;
  MatGetLocalSize( gloAssem_ptr->K, &local_row_size, &local_col_size );
  
  MatCreateShell( PETSC_COMM_WORLD, local_row_size, local_col_size,
    PETSC_DETERMINE, PETSC_DETERMINE, (void *)gloAssem_ptr, &shell_mat);

  MatShellSetOperation(shell_mat, MATOP_MULT, (void(*)(void))MF_T::MF_MatMult);

  // ===== Initialize the dot_sol vector by solving mass matrix =====
  if( is_restart == false )
  {
    SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
    PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "mass_", "mass_" );

    KSPSetType(lsolver_acce->ksp, KSPGMRES);
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
        KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );

    SI_qp->search_all_opposite_point(anchor_elementv, opposite_elementv, elements, quads, free_quad, locitf, SI_sol);

    gloAssem_ptr->Assem_mass_residual( sol, disp_mesh, locElem, locAssem_ptr, elementv,
        elements, anchor_elementv, opposite_elementv, quadv, quads, free_quad, locIEN, fNode,
        locnbc, locebc, locwbc, locitf, SI_sol, SI_qp );

    lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_sol );

    dot_sol -> ScaleValue(-1.0);

    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained. \n");
    lsolver_acce -> print_info();
    delete lsolver_acce;
    SYS_T::commPrint(" The mass matrix lsolver is destroyed.\n");
  }

  // ===== Linear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

  // ===== Nonlinear solver context =====
  PNonlinear_NS_Solver * nsolver = new PNonlinear_NS_Solver( pNode, fNode,
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq, nl_threshold );

  nsolver->print_info();

  // ===== Temporal solver context =====
  PTime_NS_Solver * tsolver = new PTime_NS_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );

  tsolver->print_info();

  // ===== Outlet data recording files =====
  for(int ff=0; ff<locebc->get_num_ebc(); ++ff)
  {
    const double dot_face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        dot_sol, locAssem_ptr, elements, quads, locebc, ff );

    const double face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        sol, locAssem_ptr, elements, quads, locebc, ff );

    const double face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure(
        sol, locAssem_ptr, elements, quads, locebc, ff );

    // set the gbc initial conditions using the 3D data
    gbc -> reset_initial_sol( ff, face_flrate, face_avepre, timeinfo->get_time(), is_restart );

    const double dot_lpn_flowrate = dot_face_flrate;
    const double lpn_flowrate = face_flrate;
    const double lpn_pressure = gbc -> get_P( ff, dot_lpn_flowrate, lpn_flowrate, timeinfo->get_time() );

    // Create the txt files and write the initial flow rates
    if(rank == 0)
    {
      std::ofstream ofile;

      // If this is NOT a restart run, generate a new file, otherwise append to
      // existing file
      if( !is_restart )
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::app );

      // If this is NOT a restart, then record the initial values
      if( !is_restart )
      {
        ofile<<"Time index"<<'\t'<<"Time"<<'\t'<<"dot Flow rate"<<'\t'<<"Flow rate"<<'\t'<<"Face averaged pressure"<<'\t'<<"Reduced model pressure"<<'\n';
        ofile<<timeinfo->get_index()<<'\t'<<timeinfo->get_time()<<'\t'<<dot_face_flrate<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
      }

      ofile.close();
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Inlet data recording files =====
  for(int ff=0; ff<locinfnbc->get_num_nbc(); ++ff)
  {
    const double inlet_face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        sol, locAssem_ptr, elements, quads, locinfnbc, ff );

    const double inlet_face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure(
        sol, locAssem_ptr, elements, quads, locinfnbc, ff );

    if( rank == 0 )
    {
      std::ofstream ofile;
      if( !is_restart )
        ofile.open( locinfnbc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( locinfnbc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::app );

      if( !is_restart )
      {
        ofile<<"Time index"<<'\t'<<"Time"<<'\t'<<"Flow rate"<<'\t'<<"Face averaged pressure"<<'\n';
        ofile<<timeinfo->get_index()<<'\t'<<timeinfo->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
      }

      ofile.close();
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

  tsolver->TM_NS_GenAlpha(is_restart, base, dot_sol, sol, disp_mesh, velo_mesh,
      tm_galpha_ptr, timeinfo, inflow_rate_ptr, pNode, locElem, locIEN, fNode,
      locnbc, locinfnbc, locrotnbc, locebc, gbc, locwbc, locitf, sir_info, SI_sol, SI_qp,
      pmat, elementv, elements, anchor_elementv, opposite_elementv,
      quadv, quads, free_quad, locAssem_ptr, gloAssem_ptr, lsolver, nsolver, shell_mat);

  // ===== Print complete solver info =====
  lsolver -> print_info();

  MatDestroy(&shell_mat);

  // ===== Clean Memory =====
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic; delete sir_info; delete locrotnbc;
  delete locElem; delete locnbc; delete locebc; delete locwbc; delete pNode; delete locinfnbc; delete locitf; delete SI_sol; delete SI_qp;
  delete tm_galpha_ptr; delete pmat; delete elementv; delete elements; delete anchor_elementv; delete opposite_elementv;
  delete quads; delete quadv; delete free_quad; delete inflow_rate_ptr; delete gbc; delete timeinfo;
  delete locAssem_ptr; delete base; delete sol; delete dot_sol; delete disp_mesh; delete velo_mesh;
  delete gloAssem_ptr; delete lsolver; delete nsolver; delete tsolver;

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
