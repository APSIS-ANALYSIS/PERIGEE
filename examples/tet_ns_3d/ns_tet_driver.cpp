// ==================================================================
// ns_tet_driver.cpp
//
// Tetrahedral element based finite element code for 3D Navier-Stokes
// equations using Variational Multiscale Formulation and Generalized
// alpha time stepping.
//
// Date: Feb. 6 2020
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"


#include "ALocal_Elem.hpp"
#include "ALocal_NodalBC.hpp"
#include "APart_Node.hpp"
#include "Matrix_PETSc.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "PDNTimeStep.hpp"
#include "PLocAssem_Tet_VMS_NS_GenAlpha.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points for tets and triangles
  int nqp_tet = 5, nqp_tri = 4;
  
  // Estimate of the nonzero per row for the sparse matrix 
  int nz_estimate = 60;

  // fluid properties
  double fluid_density = 1.065;
  double fluid_mu = 3.5e-2;
  
  // inflow file
  std::string inflow_file("inflow_fourier_series.txt");

  // LPN file
  std::string lpn_file("lpn_rcr_input.txt");

  // back flow stabilization
  double bs_beta = 0.2;

  // part file location
  std::string part_file("part");

  // determine if we want to print solver info
  bool is_LS_info = false;

  // nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits = 20;
  int nl_refreq = 4;

  // time stepping parameters
  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 1.0;
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

  // Restart options
  bool is_restart = false;
  int restart_index = 0;
  double restart_time = 0.0;
  double restart_step = 1.0e-3;
  std::string restart_name = "SOL_";

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionReal("-bs_beta", bs_beta);
  SYS_T::GetOptionReal("-fl_density", fluid_density);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionString("-inflow_file", inflow_file);
  SYS_T::GetOptionString("-lpn_file", lpn_file);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionBool("-is_ls_info", is_LS_info);
  SYS_T::GetOptionReal("-nl_rtol", nl_rtol);
  SYS_T::GetOptionReal("-nl_atol", nl_atol);
  SYS_T::GetOptionReal("-nl_dtol", nl_dtol);
  SYS_T::GetOptionInt("-nl_maxits", nl_maxits);
  SYS_T::GetOptionInt("-nl_refreq", nl_refreq);
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

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-bs_beta:", bs_beta);
  SYS_T::cmdPrint("-fl_density:", fluid_density);
  SYS_T::cmdPrint("-fl_mu:", fluid_mu);
  SYS_T::cmdPrint("-inflow_file:", inflow_file);
  SYS_T::cmdPrint("-lpn_file:", lpn_file);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-nl_rtol:", nl_rtol);
  SYS_T::cmdPrint("-nl_atol:", nl_atol);
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol);
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq);
  SYS_T::cmdPrint("-init_time:", initial_time);
  SYS_T::cmdPrint("-init_step:", initial_step);
  SYS_T::cmdPrint("-init_index:", initial_index);
  SYS_T::cmdPrint("-fina_time:", final_time);
  SYS_T::cmdPrint("-ttan_freq:", ttan_renew_freq);
  SYS_T::cmdPrint("-sol_rec_freq:", sol_record_freq);
  SYS_T::cmdPrint("-sol_name:", sol_bName);
  if(is_restart)
  {
    PetscPrintf(PETSC_COMM_WORLD, "-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else PetscPrintf(PETSC_COMM_WORLD, "-is_restart: false \n");

  if(is_LS_info) PetscPrintf(PETSC_COMM_WORLD, "-is_ls_info: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-is_ls_info: false \n");

  // ===== Data from Files =====
  // Control points' xyz coordinates
  FEANode * fNode = new FEANode(part_file, rank);

  // Local IEN array
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  // Mesh partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  // Local element indices
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // Local nodal bc
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);

  // Local inflow bc
  ALocal_Inflow_NodalBC * locinfnbc = new ALocal_Inflow_NodalBC(part_file, rank);

  // Local elemental bc
  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);

  // Nodal indices in the subdomain
  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");

  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  PetscPrintf(PETSC_COMM_WORLD,
      "===> %d processor(s) are assigned for FEM analysis. \n", size);
  
  // ===== Inflow flow rate =====
  SYS_T::commPrint("===> Setup inflow flow rate. \n");

  ICVFlowRate * inflow_rate_ptr = new CVFlowRate_Unsteady( inflow_file.c_str() );

  inflow_rate_ptr->print_info();

  // ===== GenBC =====
  IGenBC * gbc = nullptr;
  
  if( SYS_T::get_genbc_file_type( lpn_file.c_str() ) == 1  )
    gbc = new GenBC_Resistance( lpn_file.c_str() );
  else if( SYS_T::get_genbc_file_type( lpn_file.c_str() ) == 2  )
    gbc = new GenBC_RCR( lpn_file.c_str(), 1000, initial_step );
  else
    SYS_T::print_fatal( "Error: GenBC input file %s format cannot be recongnized.\n", lpn_file.c_str() );

  gbc -> print_info();
  
  // Make sure the gbc number of faces mathes that of ALocal_EBC
  SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc->get_num_ebc(),
      "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Quadrature rules =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  // ===== Finite Element Container =====
  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( quadv-> get_num_quadPts() );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0(
      quads-> get_num_quadPts() );

  // ===== Generate a sparse matrix for strong enforcement of essential BCs
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);
  pmat->gen_perm_bc(pNode, locnbc);

  // ===== Generalized-alpha ====
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  const double genA_spectrium = 0.5;
  const bool genA_is2ndSystem = false;
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_spectrium, genA_is2ndSystem);
  tm_galpha_ptr->print_info();

  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ===== Local Assembly routine =====
  IPLocAssem * locAssem_ptr = new PLocAssem_Tet_VMS_NS_GenAlpha(
      tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas(),
      fluid_density, fluid_mu, bs_beta );

  // ===== Clean Memory =====
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete locnbc; delete locebc; delete pNode; delete locinfnbc;
  delete tm_galpha_ptr; delete pmat; delete elementv; delete elements;
  delete quads; delete quadv; delete inflow_rate_ptr; delete gbc; delete timeinfo;
  delete locAssem_ptr;

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
