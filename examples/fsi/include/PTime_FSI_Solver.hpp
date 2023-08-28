#ifndef PTIME_FSI_SOLVER_HPP
#define PTIME_FSI_SOLVER_HPP
// ============================================================================
// PTime_FSI_Solver.hpp
//
// Parallel time solver for time integration.
//
// Date: Jan. 13 2022
// Author: Ju Liu
// ============================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_FSI_Solver.hpp"

class PTime_FSI_Solver
{
  public:
    PTime_FSI_Solver( const std::string &input_name, const int &input_record_freq, 
        const int &input_renew_tang_freq, const double &input_final_time );

    ~PTime_FSI_Solver();

    void print_info() const;

    void TM_FSI_GenAlpha(
        const bool &restart_init_assembly_flag,
        const IS &is_v,
        const IS &is_p,
        const PDNSolution * const &sol_base,
        const PDNSolution * const &init_dot_disp,
        const PDNSolution * const &init_dot_velo,
        const PDNSolution * const &init_dot_pres,
        const PDNSolution * const &init_disp,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_pres,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ICVFlowRate * const flr_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_InflowBC * const &infnbc,
        const ALocal_NBC * const &nbc_mesh,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const ALocal_EBC * const &ebc_mesh,
        IGenBC * const &gbc,
        const Matrix_PETSc * const &bc_mat,
        const Matrix_PETSc * const &bc_mesh_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const Tissue_prestress * const &ps_ptr,
        IPLocAssem_2x2Block * const &lassem_fluid_ptr,
        IPLocAssem_2x2Block * const &lassem_solid_ptr,
        IPLocAssem * const &lassem_mesh_ptr,
        IPGAssem * const &gassem_ptr,
        IPGAssem * const &gassem_mesh_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
        const PNonlinear_FSI_Solver * const &nsolver_ptr ) const;

    void TM_FSI_Prestress(
        const bool &is_record_sol_flag,
        const double &prestress_tol,
        const IS &is_v,
        const IS &is_p,
        const PDNSolution * const &init_dot_disp,
        const PDNSolution * const &init_dot_velo,
        const PDNSolution * const &init_dot_pres,
        const PDNSolution * const &init_disp,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_pres,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        Tissue_prestress * const &ps_ptr,
        IPLocAssem_2x2Block * const &lassem_solid_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        const PNonlinear_FSI_Solver * const &nsolver_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const std::string &middle_name, 
        const int &counter ) const;

    std::string Name_dot_Generator( const std::string &middle_name, 
        const int &counter ) const;

    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;

    void Nullify_solid_dof( const APart_Node * const &pnode,
        const int &in_dof, PDNSolution * const &sol ) const;
};

#endif
