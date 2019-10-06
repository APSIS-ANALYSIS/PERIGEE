#ifndef PNONLINEAR_SEG_SOLVER_HPP
#define PNONLINEAR_SEG_SOLVER_HPP
// ==================================================================
// PNonlinear_Seg_Solver.hpp
// 
// This class implements the segregated algorithm for elastodynamics.
// This file resides in the tet4_seg_elastodynamic project folder.
//
// Author: Ju Liu
// Date: May 23 2017
// ==================================================================
#include "ICVFlowRate.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_P_V_Mixed_3D.hpp"
#include "PDNSolution_U_Mixed_3D.hpp"
#include "Seg_Sol_Tools.hpp"

class PNonlinear_Seg_Solver
{
  public:
    PNonlinear_Seg_Solver( const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const double &input_nrtol, 
        const double &input_natol, const double &input_ndtol,
        const int &input_max_iteration, const int &input_renew_freq );

    ~PNonlinear_Seg_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    // --------------------------------------------------------------
    // GenAlpha_Seg_solve_ALE_NS:
    // This is a ALE solver for fluid dynamics.
    //
    // This solver solves the Navier-Stokes using full Generalized
    // alpha method (by "full", I mean that the pressure is evaluated
    // at time step n+alpha_f as well.)
    //
    // The mesh update is solved in a segregated algorithm, meaning
    // that the mesh is updated by solving a harmonic extension problem
    // only. Then we continue to solve NS to achieve nonlinear
    // convergence.
    // --------------------------------------------------------------
    void GenAlpha_Seg_solve_ALE_NS(
        const bool &new_tangent_flag,
        const bool &is_ale_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &sol_base,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ICVFlowRate * const flr_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_Inflow_NodalBC * const &infnbc_part,
        const ALocal_NodalBC * const &nbc_mesh_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_mesh_part,
        const IGenBC * const &gbc,
        const Matrix_PETSc * const &bc_mat,
        const Matrix_PETSc * const &bc_mesh_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_solid_ptr,
        IPLocAssem * const &lassem_mesh_ptr,
        IPGAssem * const &gassem_ptr,
        IPGAssem * const &gassem_mesh_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        bool &conv_flag, int &nl_counter ) const;


    // --------------------------------------------------------------
    // GenAlpha_Seg_solve_FSI:
    // This is the coupled FSI solver.
    // 
    // The solid and fluid are discretized by VMS and generalized-alpha.
    //
    // The mesh update is solved in a segregated algorithm, meaning
    // that the mesh is updated by solving a harmonic extension problem
    // only. Then we continue to solve NS to achieve nonlinear
    // convergence.
    // --------------------------------------------------------------
    void GenAlpha_Seg_solve_FSI(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &sol_base,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ICVFlowRate * const flr_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_Inflow_NodalBC * const &infnbc_part,
        const ALocal_NodalBC * const &nbc_mesh_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_mesh_part,
        const Matrix_PETSc * const &bc_mat,
        const Matrix_PETSc * const &bc_mesh_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_fluid_ptr,
        IPLocAssem * const &lassem_solid_ptr,
        IPLocAssem * const &lassem_mesh_ptr,
        IPGAssem * const &gassem_ptr,
        IPGAssem * const &gassem_mesh_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol;
    const double na_tol;
    const double nd_tol;
    const int nmaxits;
    const int nrenew_freq;

    // vector container for the step update in the smaller matrix problem
    PDNSolution * dot_P_V_step;

    // vector container for the mesh motion solution
    PDNSolution * mesh_disp;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {PetscPrintf(PETSC_COMM_WORLD,
        "  === NR ite: %d, r_error: %e, a_error: %e \n",
        count, rel_err, abs_err);}

    void rescale_inflow_value( const double &stime,
        const ALocal_Inflow_NodalBC * const &infbc,
        const ICVFlowRate * const &flrate,
        const PDNSolution * const &sol_base,
        PDNSolution * const &sol ) const;
};

#endif
