#ifndef PNONLINEAR_CMM_SOLVER_HPP
#define PNONLINEAR_CMM_SOLVER_HPP
// ==================================================================
// PNonlinear_CMM_Solver.hpp
// 
// Parallel nonlinear solver for Navier-Stokes equations. 
//
// Author: Ju Liu
// Date: Feb 11 2020
// ==================================================================
#include "ICVFlowRate.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc_CMM.hpp"
#include "PDNSolution_NS.hpp"
#include "PDNSolution_Wall_Disp.hpp"

class PNonlinear_CMM_Solver
{
  public:
    PNonlinear_CMM_Solver( const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const double &input_nrtol, const double &input_natol, 
        const double &input_ndtol, const int &input_max_iteration, 
        const int &input_renew_freq, 
        const int &input_renew_threshold = 4,
        const bool &prestress_flag = false,
        const double &ps_disp_atol = 1.0e-6 );

    ~PNonlinear_CMM_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    // --------------------------------------------------------------
    // GenAlpha_Solve_CMM:
    // This is a solver for CMM-FSI using the 2nd-order Generalized.
    // alpha method.
    // --------------------------------------------------------------
    void GenAlpha_Solve_CMM(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &sol_base,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        const PDNSolution * const &pre_dot_sol_wall_disp,
        const PDNSolution * const &pre_sol_wall_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ICVFlowRate * const flr_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_Inflow_NodalBC * const &infnbc_part,
        const ALocal_Ring_NodalBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        ALocal_EBC * const &ebc_wall_part,
        const IGenBC * const &gbc,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementw,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        PDNSolution * const &dot_sol_wall_disp,
        PDNSolution * const &sol_wall_disp,
        bool &prestress_conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    // flag for whether the wall prestress is being solved for and updated 
    const bool solve_prestress;

    // tolerance for displacement L2 norm when solving for wall prestress 
    const double prestress_tol;

    // vector container for the step update in the smaller matrix problem
    PDNSolution * dot_step;

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

    // CMM segregated algorithm: nodal update of sol_wall_disp and 
    // dot_sol_wall_disp based on the dot_velo increment 
    // The dot_step vector is assumed to take the following form:
    //     [dot_pres_step, dot_velo_x_step, dot_velo_y_step, dot_velo_z_step] 
    // The wall_data is a vector of three degrees of freedom:
    //     [wall_data_x, wall_data_y, wall_data_z]
    // and it can be either sol_wall_disp or dot_sol_wall_disp
    // The update will occur in the following manner:
    //     wall_data_x += val * dot_velo_x_step
    //     wall_data_y += val * dot_velo_y_step
    //     wall_data_z += val * dot_velo_z_step
    void update_wall( const double &val,
        const PDNSolution * const &dot_step,
        const ALocal_EBC * const &ebc_wall_part, 
        PDNSolution * const &wall_data ) const;


    // Check whether the ring BC constraints are properly satisfied
    // by printing their evaluations 
    void compute_ringbc_constraints(
        const PDNSolution * const &sol,
        const PDNSolution * const &sol_wall_disp,
        const ALocal_Ring_NodalBC * const &ringnbc_part ) const;
};

#endif
