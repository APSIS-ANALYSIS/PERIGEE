// ==================================================================
// FEM_analysis_main.cpp
// ------------------------------------------------------------------
// Objects:
// Perform parallel 3d dynamic nonlinear finite element (isogeometric)
// analysis.
//
// Input:
// partition file in hdf5 format.
//
// Output:
// finite element solution vector.
//
// Note:
// This code relies on the PETSc library.
//
// Date:
// Oct. 15 2013
//
// Author:
// Ju Liu, Ph.D. candidate, the University of Texas at Austin
// ==================================================================
#include <cmath>
#include <iomanip>
#include "Vec_Tools.hpp"
#include "QuadPts_Gauss.hpp"
#include "HDF5_PartReader.hpp"
#include "HDF5_Writer.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEANode.hpp"
#include "AExtractor_3D_NURBS_xyz.hpp"
#include "FEAElement_NURBS_3D_der1_v3.hpp"
#include "AGlobal_Mesh_Info_1Patch_NURBS_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_BC_3D.hpp"
#include "IALocal_meshSize.hpp"
#include "ALocal_meshSize_3D_NURBS.hpp"
#include "AInt_Weight.hpp"
#include "APart_Node.hpp"
#include "PDNSolution.hpp"
#include "PDNTimeStep.hpp"
#include "PDNSolution_heatEqn.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "PLocAssem_NLHeat_3D_GenAlpha.hpp"
#include "PGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "PNonlinear_Solver.hpp"
#include "PTime_Solver.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points
  int nqpx = 3; int nqpy = 3; int nqpz = 3;

  // partition file base name
  std::string part_file("part");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 0.9;
  int nl_maxits = 20;
  int nl_refreq = 4;
  
  // Time step initailization
  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 1.0;

  // Time solver parameters
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

  PetscMPIInt rank, size;
  // ======= PETSc Initialize =======
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ======= Read Command Line Arguments =======
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  
  SYS_T::GetOptionInt("-nqpx", nqpx);
  SYS_T::GetOptionInt("-nqpy", nqpy);
  SYS_T::GetOptionInt("-nqpz", nqpz);
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
  SYS_T::cmdPrint("-nqpz:", nqpz); SYS_T::cmdPrint("-nl_rtol:", nl_rtol); 
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
  // 1.1 Get control points
  FEANode * fNode = new FEANode(part_file, rank);

  // 1.2 Get mesh size for each local element
  IALocal_meshSize * locmSize = new ALocal_meshSize_3D_NURBS(h5reader);

  // 1.3 Get extraction operator for each local elements
  IAExtractor * fExt = new AExtractor_3D_NURBS_xyz(h5reader);

  // 1.4 Get LIEN for each local elements
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // 1.5 Get Global Mesh Info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_1Patch_NURBS_3D(h5reader);

  // 1.6 Get partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  // 1.7 Get local element info
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // 1.8 Get local BC info
  IALocal_BC * locBC = new ALocal_BC_3D(h5reader);

  // 1.9 Get node partition info
  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::synPrint(" done\t", rank);
  delete h5reader;

  if(size != PartBasic->get_cpu_size())
    MPI_Abort(PETSC_COMM_WORLD, 1);

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for FEM analysis. ", size);

  // ======= Generate Finite Element =======
  SYS_T::commPrint("\n===> Build quadrature rules ... \n");
  IQuadPts * quad_z = new QuadPts_Gauss(nqpz); 
  IQuadPts * quad_y = new QuadPts_Gauss(nqpy);  
  IQuadPts * quad_x = new QuadPts_Gauss(nqpx);

  SYS_T::commPrint("===> Build quadrature weight ... \n");
  AInt_Weight * Int_w = new AInt_Weight(quad_x, quad_y, quad_z);

  SYS_T::commPrint("===> Build univariate Bezier elements ... \n");
  BernsteinBasis_Array Bena_x(GMIptr->get_xdegree(), quad_x);
  BernsteinBasis_Array Bena_y(GMIptr->get_ydegree(), quad_y);
  BernsteinBasis_Array Bena_z(GMIptr->get_zdegree(), quad_z);

  SYS_T::commPrint("===> Build shape functions ... \n");
  std::vector<FEAElement *> elemArray; 
  elemArray.resize(locElem->get_nlocalele());
  double feaelement_memsize = 0.0; clock_t elem_timer = clock();

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    FEAElement * elem_ptr = new FEAElement_NURBS_3D_der1_v3(ee, locmSize,
        &Bena_x, &Bena_y, &Bena_z, fNode, fExt, locIEN);
    elemArray[ee] = elem_ptr;
    feaelement_memsize += elem_ptr->get_memory_usage(); 
  }
  elem_timer = clock() - elem_timer;
  MPI_Barrier(PETSC_COMM_WORLD);
  SYS_T::synPrintElementInfo(locElem->get_nlocalele(), feaelement_memsize,
      (double)elem_timer/(double)CLOCKS_PER_SEC, rank);
  // ---------------------------------------

  // ======= Finite Element Analysis =======
  // 2.1 Solution Initialization
  PDNSolution * disp = new PDNSolution_heatEqn(pNode, locBC, 0);
  PDNSolution * velo = new PDNSolution_heatEqn(pNode, locBC, 0);

  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // 2.2 Local Assembly pointer
  SYS_T::commPrint("===> Genereate the Generalized-alpha time scheme ... \n");
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(0.5);

  SYS_T::commPrint("===> Initialize local assembly routine ... \n");
  IPLocAssem * locAssem_ptr = new PLocAssem_NLHeat_3D_GenAlpha(
      tm_galpha_ptr, GMIptr->get_nLocBas(), Int_w->get_num() );

  // 2.3 Global Assembly pointer
  int vpetsc_type = 0; // petsc version controller
  PGAssem * gloAssem_ptr = new PGAssem(locAssem_ptr, GMIptr, pNode, vpetsc_type);

  // 2.4 Estimate the matrix structure
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN, pNode, locBC );
  gloAssem_ptr->Fix_nonzero_err_str();
  SYS_T::commPrint("===> Matrix nonzero structure fixed ... \n");

  // 2.5 Setup linear solver context
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();
  SYS_T::commPrint("===> PETSc linear solver setted up:\n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
  lsolver->Info();
  SYS_T::commPrint("----------------------------------------------------------- \n");

  // 2.6 Assembly mass matrix and solve for consistent initial solution
  gloAssem_ptr->Clear_KG();
  gloAssem_ptr->Assem_mass_residual( disp, locElem, locAssem_ptr, locIEN, pNode,
      fNode, Int_w, elemArray, locBC );

  lsolver->Solve( gloAssem_ptr->K, gloAssem_ptr->G, velo); 
  SYS_T::commPrint("initial solution's time derivative obtained. \n"); 

  // 2.7 Setup nonlinear solver context
  PNonlinear_Solver * nsolver = new PNonlinear_Solver(nl_rtol, nl_atol,
      nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->Info();

  // 2.8 Setup time marching context
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
  SYS_T::commPrint("\n===> Clean memory ... \n");
  delete fNode;
  delete fExt;
  delete GMIptr;
  delete PartBasic;
  delete locIEN;
  delete locElem;
  delete locBC;
  delete locmSize;
  delete quad_z; delete quad_y; delete quad_x;
  delete pNode;
  delete Int_w;
  delete disp;
  delete velo;
  delete timeinfo;
  delete tm_galpha_ptr;
  delete locAssem_ptr;
  delete gloAssem_ptr;
  delete lsolver;
  delete nsolver;
  delete tsolver;
  std::vector<FEAElement *>::iterator it_elema;
  for(it_elema = elemArray.begin(); it_elema != elemArray.end(); ++it_elema)
    delete *it_elema;
  // ---------------------------------------

  SYS_T::commPrint("===> Exit program. \n");
  PetscFinalize();
  return 0;
}
//EOF
