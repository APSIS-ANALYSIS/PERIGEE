#ifndef SOLIDS_INITHELPERS_HPP
#define SOLIDS_INITHELPERS_HPP
// ============================================================================
// InitHelpers.hpp
//
// Initial-condition and restart helpers for mixed solid dynamics.
// ============================================================================
#include "PDNSolution.hpp"
#include "PGAssem_Solid_FEM.hpp"
#include "PLinear_Solver_PETSc.hpp"

namespace SOLID_INIT
{
  inline void zero_solution( PDNSolution * const &sol )
  {
    VecSet(sol->solution, 0.0);
    sol->Assembly_GhostUpdate();
  }

  inline void read_binary_checked( PDNSolution * const &sol,
      const std::string &name )
  {
    SYS_T::file_check(name);
    sol->ReadBinary(name);
  }

  inline void initialize_solution_state( const APart_Node * const pNode,
      const bool is_restart, const int restart_index,
      const double restart_time, const double restart_step,
      const std::string &restart_disp_name,
      const std::string &restart_velo_name,
      const std::string &restart_pres_name,
      std::unique_ptr<PDNSolution> &disp,
      std::unique_ptr<PDNSolution> &velo,
      std::unique_ptr<PDNSolution> &pres,
      std::unique_ptr<PDNSolution> &dot_disp,
      std::unique_ptr<PDNSolution> &dot_velo,
      std::unique_ptr<PDNSolution> &dot_pres,
      int &initial_index, double &initial_time, double &initial_step )
  {
    disp = SYS_T::make_unique<PDNSolution>(
        PDNSolution::Gen_zero( pNode, 3 ) );
    velo = SYS_T::make_unique<PDNSolution>(
        PDNSolution::Gen_zero( pNode, 3 ) );
    pres = SYS_T::make_unique<PDNSolution>(
        PDNSolution::Gen_zero( pNode, 1 ) );

    dot_disp = SYS_T::make_unique<PDNSolution>(
        PDNSolution::Gen_zero( pNode, 3 ) );
    dot_velo = SYS_T::make_unique<PDNSolution>(
        PDNSolution::Gen_zero( pNode, 3 ) );
    dot_pres = SYS_T::make_unique<PDNSolution>(
        PDNSolution::Gen_zero( pNode, 1 ) );

    if(is_restart)
    {
      initial_index = restart_index;
      initial_time  = restart_time;
      initial_step  = restart_step;

      read_binary_checked( disp.get(), restart_disp_name );
      read_binary_checked( velo.get(), restart_velo_name );
      read_binary_checked( pres.get(), restart_pres_name );

      const std::string restart_dot_disp_name = "dot_" + restart_disp_name;
      const std::string restart_dot_velo_name = "dot_" + restart_velo_name;
      const std::string restart_dot_pres_name = "dot_" + restart_pres_name;

      read_binary_checked( dot_disp.get(), restart_dot_disp_name );
      read_binary_checked( dot_velo.get(), restart_dot_velo_name );
      read_binary_checked( dot_pres.get(), restart_dot_pres_name );

      SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
      SYS_T::commPrint("     restart_disp_name: %s \n",
          restart_disp_name.c_str());
      SYS_T::commPrint("     restart_velo_name: %s \n",
          restart_velo_name.c_str());
      SYS_T::commPrint("     restart_pres_name: %s \n",
          restart_pres_name.c_str());
      SYS_T::commPrint("     restart_dot_disp_name: %s \n",
          restart_dot_disp_name.c_str());
      SYS_T::commPrint("     restart_dot_velo_name: %s \n",
          restart_dot_velo_name.c_str());
      SYS_T::commPrint("     restart_dot_pres_name: %s \n",
          restart_dot_pres_name.c_str());
      SYS_T::commPrint("     restart_time: %e \n", restart_time);
      SYS_T::commPrint("     restart_index: %d \n", restart_index);
      SYS_T::commPrint("     restart_step: %e \n", restart_step);
    }
  }

  inline void initialize_dot_solution( PGAssem_Solid_FEM * const gloAssem,
      const IS &is_velo, const IS &is_pres,
      PDNSolution * const &dot_disp,
      PDNSolution * const &dot_velo,
      PDNSolution * const &dot_pres,
      PDNSolution * const &disp,
      PDNSolution * const &velo,
      PDNSolution * const &pres,
      const bool is_restart )
  {
    if(is_restart) return;

    SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
    auto lsolver_acce = SYS_T::make_unique<PLinear_Solver_PETSc>(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "mass_", "mass_" );

    KSPSetType(lsolver_acce->ksp, KSPGMRES);
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
        KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );

    gloAssem->Assem_mass_residual( disp, velo, pres );

    Vec dot_vp;
    VecDuplicate(gloAssem->G, &dot_vp);

    lsolver_acce->Solve( gloAssem->K, gloAssem->G, dot_vp );
    VecScale(dot_vp, -1.0);

    Vec sol_v, sol_p;
    VecGetSubVector(dot_vp, is_velo, &sol_v);
    VecGetSubVector(dot_vp, is_pres, &sol_p);

    zero_solution( dot_velo );
    zero_solution( dot_pres );
    dot_velo->PlusAX( sol_v, 1.0 );
    dot_pres->PlusAX( sol_p, 1.0 );

    VecRestoreSubVector(dot_vp, is_velo, &sol_v);
    VecRestoreSubVector(dot_vp, is_pres, &sol_p);
    VecDestroy(&dot_vp);

    dot_disp->Copy( velo );

    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained. \n");
    lsolver_acce->print_info();
    SYS_T::commPrint(" The mass matrix lsolver is destroyed.\n");
  }
}

#endif
