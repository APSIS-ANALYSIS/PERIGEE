#ifndef PNONLINEAR_TET4_NS_3D_SOLVER_HPP
#define PNONLINEAR_TET4_NS_3D_SOLVER_HPP
// ==================================================================
// PNonlinear_Tet4_NS_3D_Solver.hpp
//
// This class implements the Tet4 element-based NS 3D Solver, nonlinear
// solver part.
//
// Author: Ju Liu
// Date Created: June 19 2017
// ==================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "PGAssem_v360_FEM.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_Tet4_NS_3D.hpp"
#include "ALocal_Inflow_NodalBC.hpp"

class PNonlinear_Tet4_NS_3D_Solver
{
  public:
    PNonlinear_Tet4_NS_3D_Solver( 
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const double &input_nrtol,
        const double &input_natol, const double &input_ndtol,
        const int &input_max_iteration, const int &input_renew_freq );

    ~PNonlinear_Tet4_NS_3D_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    // --------------------------------------------------------------
    // Generalized-alpha for velocity and Newton for pressure
    // CMAME 197 2007 173-201
    // --------------------------------------------------------------
    void GenAlpha_VMS_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        bool &conv_flag, int &nl_counter ) const;


    // --------------------------------------------------------------
    // This function passes a base sol vector that defines the value
    // on the inflow surface. This is used to adjust the inflow bc
    // values w.r.t. time as a function of time.
    // --------------------------------------------------------------
    void GenAlpha_VMS_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &sol_base,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_Inflow_NodalBC * const &infnbc_part,
        const ALocal_EBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol;
    const double na_tol;
    const double nd_tol;
    const int nmaxits;
    const int nrenew_freq;

    // update solution
    PDNSolution * velo_step;

    // dot_sol at alpha_m
    PDNSolution * dot_sol_alpha;

    // sol at alpha_f
    PDNSolution * sol_alpha;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {PetscPrintf(PETSC_COMM_WORLD, 
        "  === NR ite: %d, r_error: %e, a_error: %e \n",
        count, rel_err, abs_err);}

    void rescale_inflow_value( const double &stime,
        const ALocal_Inflow_NodalBC * const &infbc,
        const PDNSolution * const &sol_base,
         PDNSolution * const &sol ) const;
};

#endif
