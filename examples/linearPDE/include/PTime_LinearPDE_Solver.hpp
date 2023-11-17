#ifndef PTIME_LINEARPDE_SOLVER_HPP
#define PTIME_LINEARPDE_SOLVER_HPP
// ==================================================================
// PTime_LinearPDE_Solver.hpp
//
// Parallel time solver for linear PDE.
//
// Date: Oct. 26 2023
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_LinearPDE_Solver.hpp"

class PTime_LinearPDE_Solver
{
  public:
    PTime_LinearPDE_Solver( const std::string &input_name,
        const int &input_record_freq, const int &input_renew_tang_freq,
        const double &input_final_time );

    ~PTime_LinearPDE_Solver();

    void print_info() const;

    void TM_GenAlpha_Transport(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_dot_sol,
        const PDNSolution * const &init_sol,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PNonlinear_LinearPDE_Solver * const &nsolver_ptr ) const;
    
    void TM_GenAlpha_Elastodynamics(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_dot_disp,
        const PDNSolution * const &init_dot_velo,
        const PDNSolution * const &init_disp,
        const PDNSolution * const &init_velo,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PNonlinear_LinearPDE_Solver * const &nsolver_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const std::string &middle_name,
        const int &counter ) const;

    std::string Name_dot_Generator( const std::string &middle_name,
        const int &counter ) const;
};

#endif
