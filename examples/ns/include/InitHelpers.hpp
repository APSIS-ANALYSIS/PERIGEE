#ifndef NS_INITHELPERS_HPP
#define NS_INITHELPERS_HPP

#include "IPGAssem.hpp"
#include "PDNSolution_NS.hpp"
#include "PLinear_Solver_PETSc.hpp"

namespace NS_INIT
{
  inline void initialize_solution_state(const APart_Node * const pNode,
      const bool is_restart, const int restart_index, const double restart_time,
      const double restart_step, const std::string &restart_name,
      std::unique_ptr<PDNSolution> &sol, std::unique_ptr<PDNSolution> &dot_sol,
      int &initial_index, double &initial_time, double &initial_step)
  {
    sol = SYS_T::make_unique<PDNSolution_NS>(pNode, 0);
    dot_sol = SYS_T::make_unique<PDNSolution_NS>(pNode, 0);

    if(is_restart)
    {
      initial_index = restart_index;
      initial_time  = restart_time;
      initial_step  = restart_step;

      SYS_T::file_check(restart_name);
      sol->ReadBinary(restart_name);

      const std::string restart_dot_name = "dot_" + restart_name;
      SYS_T::file_check(restart_dot_name);
      dot_sol->ReadBinary(restart_dot_name);

      SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
      SYS_T::commPrint("     restart_name: %s \n", restart_name.c_str());
      SYS_T::commPrint("     restart_dot_name: %s \n", restart_dot_name.c_str());
      SYS_T::commPrint("     restart_time: %e \n", restart_time);
      SYS_T::commPrint("     restart_index: %d \n", restart_index);
      SYS_T::commPrint("     restart_step: %e \n", restart_step);
    }
  }

  inline void initialize_dot_solution(IPGAssem * const gloAssem,
      PDNSolution * const sol, PDNSolution * const dot_sol, const bool is_restart)
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

    gloAssem->Assem_mass_residual(sol);
    lsolver_acce->Solve(gloAssem->K, gloAssem->G, dot_sol);
    dot_sol->ScaleValue(-1.0);

    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained. \n");
    lsolver_acce->print_info();
    SYS_T::commPrint(" The mass matrix lsolver is destroyed.\n");
  }
}

#endif
