// ==================================================================
// heat_2d_genAlpha.cpp
// ------------------------------------------------------------------
// Objects:
// Perform parallel 2d nonlinear heat equation FEM/IGA dynamic 
// analysis.
//
// Input:
// Partition files in hdf5 format.
//
// Output:
// finite element solution vector in PETSc binary format.
//
// Date:
// April 12 2014
//
// Author:
// Ju Liu, Ph.D. candidate, the University of Texas at Austin
// ==================================================================
#include <cmath>
#include "HDF5_PartReader.hpp"
#include "AExtractor_2D_NURBS_xy.hpp"
#include "ALocal_meshSize_2D_NURBS.hpp"
#include "AGlobal_Mesh_Info_1Patch_NURBS_2D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_BC_2D.hpp"
#include "QuadPts_Gauss.hpp"
#include "FEAElement_NURBS_2D_der2.hpp"
#include "PDNSolution_heatEqn.hpp"
#include "PLocAssem_NLHeat_2D_GenAlpha.hpp"
#include "PTime_Solver.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature pts
  int nqpx = 3; int nqpy = 3;

  // Partition file base name, use default "part"
  std::string part_file("part");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-5;
  double nl_dtol = 0.9;
  int nl_maxits = 20;
  int nl_refreq = 4;

  // Time stepping parameters
  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 1.0;

  // Time solver parameters
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

  // ======= PETSc Initialization =======
  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ======= Command-Line arguments =======
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-nqpx", nqpx);
  SYS_T::GetOptionInt("-nqpy", nqpy);
  SYS_T::GetOptionString("-part_file", part_file);

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

  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-nqpx:", nqpx); SYS_T::cmdPrint("-nqpy:",nqpy);
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

  // ======= Generate Main Data Structure =======
  SYS_T::commPrint("===> Reading mesh files ... \n");
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);

  // GMDS.1 Control Points
  FEANode * fNode = new FEANode(part_file, rank);

  // GMDS.2 Local mesh size
  IALocal_meshSize * locmSize = new ALocal_meshSize_2D_NURBS(h5reader);    

  // GMDS.3 Extractiion operator for local elements
  IAExtractor * fExt = new AExtractor_2D_NURBS_xy(h5reader);

  // GMSD.4 LIEN
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // GMSD.5 Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_1Patch_NURBS_2D(h5reader);

  // GMSD.6 Partition parameters
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  // GMSD.7 Local element info
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // GMSD.8 BC
  IALocal_BC * locBC = new ALocal_BC_2D(h5reader);

  // GMSD.9 Local node info
  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::synPrint(" done.\t", rank);
  delete h5reader;

  if( size != PartBasic->get_cpu_size() )
  {
    PetscPrintf(PETSC_COMM_WORLD,
        "Error: Assigned CPU number does not match the partition number. \n");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ======= Generate quadrature data =======
  SYS_T::commPrint("===> Build quadrature rules ... \n");
  IQuadPts * quad_x = new QuadPts_Gauss(nqpx);
  IQuadPts * quad_y = new QuadPts_Gauss(nqpy);

  SYS_T::commPrint("===> Build quadrature weights ... \n");
  AInt_Weight * Int_w = new AInt_Weight(quad_x, quad_y);

  SYS_T::commPrint("===> Build univariate Bezier elements ... \n");
  BernsteinBasis_Array Bena_x(GMIptr->get_xdegree(), quad_x);
  BernsteinBasis_Array Bena_y(GMIptr->get_ydegree(), quad_y);

  SYS_T::commPrint("===> Build shape functions ... \n");
  std::vector<FEAElement *> elemArray;
  elemArray.resize(locElem->get_nlocalele());
  double feaelement_memsize = 0.0;
  clock_t elem_timer = clock();

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    FEAElement * elem_ptr = new FEAElement_NURBS_2D_der2( ee, locmSize,
        &Bena_x, &Bena_y, fNode, fExt, locIEN );
    elemArray[ee] = elem_ptr;
    feaelement_memsize += elem_ptr->get_memory_usage();
  }
  elem_timer = clock() - elem_timer;
  MPI_Barrier(PETSC_COMM_WORLD);
  SYS_T::synPrintElementInfo(locElem->get_nlocalele(), feaelement_memsize,
      (double)elem_timer/(double)CLOCKS_PER_SEC, rank);

  // ======= Finite Element Analysis =======
  // FEA.1 Initial solution
  PDNSolution * disp = new PDNSolution_heatEqn(pNode, locBC, 0);
  PDNSolution * velo = new PDNSolution_heatEqn(pNode, locBC, 0);

  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // FEA.2 Local assembly setup
  SYS_T::commPrint("===> Genereate the Generalized-alpha time scheme ... \n");
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(0.5);

  SYS_T::commPrint("===> Initialize local assembly routine ... \n");
  IPLocAssem * locAssem_ptr = new PLocAssem_NLHeat_2D_GenAlpha(
      tm_galpha_ptr, GMIptr->get_nLocBas(), Int_w->get_num() );

  // FEA.3 Globaly assembly setup
  int vpetsc_type = 0;
  PGAssem * gloAssem_ptr = new PGAssem(locAssem_ptr, GMIptr, pNode, vpetsc_type);

  // FEA.4 Estimate nonzero structure
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN, pNode, locBC );
  gloAssem_ptr->Fix_nonzero_err_str();
  SYS_T::commPrint("===> Matrix nonzero structure fixed ... \n");

  // FEA.5 Linear solver setup
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();
  SYS_T::commPrint("===> PETSc linear solver setted up:\n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
  lsolver->Info();
  SYS_T::commPrint("----------------------------------------------------------- \n");

  // FEA.6 Solve for consistent initial condition
  gloAssem_ptr->Clear_KG();
  gloAssem_ptr->Assem_mass_residual( disp, locElem, locAssem_ptr, locIEN, pNode,
      fNode, Int_w, elemArray, locBC );

  lsolver->Solve( gloAssem_ptr->K, gloAssem_ptr->G, velo);
  SYS_T::commPrint("\n===> Initial solution's time derivative obtained. \n");

  // FEA.7 Nonlinear solver steup
  PNonlinear_Solver * nsolver = new PNonlinear_Solver(nl_rtol, nl_atol,
      nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->Info();

  // FEA.8 Time solver steup
  PTime_Solver * tsolver = new PTime_Solver( sol_bName, sol_record_freq,
      ttan_renew_freq, final_time );

  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->Info();

  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_generalized_alpha(
      velo, disp, timeinfo, tm_galpha_ptr, locElem, locIEN, pNode,
      fNode, locBC, Int_w, elemArray, locAssem_ptr, gloAssem_ptr,
      lsolver, nsolver );

  // ======= PETSc Finalize =======
  SYS_T::commPrint("===> Clean memory ...\n");
  delete fNode;   delete locmSize;  delete fExt;
  delete locIEN;  delete GMIptr;    delete PartBasic;
  delete locElem; delete locBC;     delete pNode;

  delete quad_x; delete quad_y; delete Int_w;

  delete disp; delete velo; delete timeinfo; delete tm_galpha_ptr;
  delete locAssem_ptr; delete gloAssem_ptr;

  delete lsolver; delete nsolver; delete tsolver;

  std::vector<FEAElement *>::iterator it_elema;
  for(it_elema = elemArray.begin(); it_elema != elemArray.end(); ++it_elema)
    delete *it_elema;

  PetscFinalize();
  return 0;
}

// EOF
