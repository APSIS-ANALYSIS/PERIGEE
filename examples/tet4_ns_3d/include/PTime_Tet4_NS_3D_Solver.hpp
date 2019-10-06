#ifndef PTIME_TET4_NS_3D_SOLVER_HPP
#define PTIME_TET4_NS_3D_SOLVER_HPP
// ==================================================================
// PTime_Tet4_NS_3D_Solver.hpp
//
// Parallel time solver for Tet4-based 3D NS solver.
//
// Author: Ju Liu
// Date Created: June 19 2017
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_Tet4_NS_3D_Solver.hpp"

class PTime_Tet4_NS_3D_Solver
{
  public:
    PTime_Tet4_NS_3D_Solver( const std::string &input_name,
        const int &input_record_freq, const int &input_renew_tang_freq,
        const double &input_final_time );

    ~PTime_Tet4_NS_3D_Solver();

    void print_info() const;

    void TM_VMS_GenAlpha(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
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
        PNonlinear_Tet4_NS_3D_Solver * const &nsolver_ptr ) const;

    // --------------------------------------------------------------
    // Solver with inflow nodal bc information passed in, so I can
    // vary the inflow condition.
    // --------------------------------------------------------------
    void TM_VMS_GenAlpha(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &inflow_base,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
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
        PNonlinear_Tet4_NS_3D_Solver * const &nsolver_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const int &counter ) const;
};

#endif
