#ifndef PNONLINEAR_FSI_SOLVER_HPP
#define PNONLINEAR_FSI_SOLVER_HPP
// ============================================================================
// PNonlinear_FSI_Solver.hpp
//
// This is the nonlienar solver for FSI problems.
//
// Author: Ju Liu
// Date: Jan 4 2022
// ============================================================================
#include "ICVFlowRate.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_V.hpp"
#include "PDNSolution_P.hpp"

class PNonlinear_FSI_Solver
{
  public:
    PNonlinear_FSI_Solver( const double &input_nrtol, const double &input_natol,
        const double &input_ndtol, const int &input_max_iteration,
        const int &input_renew_freq, const int &input_renew_thred );

    ~PNonlinear_FSI_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    void GenAlpha_Seg_solve_FSI(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const IS &is_v,
        const IS &is_p,
        const PDNSolution * const &pre_dot_disp,
        const PDNSolution * const &pre_dot_velo,
        const PDNSolution * const &pre_dot_pres,
        const PDNSolution * const &pre_disp,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_pres,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ICVFlowRate * const flr_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &feanode_ptr,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        ALocal_InflowBC * const &infnbc_part,
        const ALocal_NBC * const &nbc_mesh,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_mesh_part,
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
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &dot_pres,
        PDNSolution * const &disp,
        PDNSolution * const &velo,
        PDNSolution * const &pres,
        bool &conv_flag, int &nl_counter,
        const std::vector<double> &inflow_area,
        const std::vector<Vector_3> &inlet_centroid,
        const std::vector<std::vector<Vector_3>> &inlet_ring ) const;

    void GenAlpha_Seg_solve_Prestress(
        const bool &new_tangent_flag,
        const double &prestress_tol,
        const double &curr_time,
        const double &dt,
        const IS &is_v,
        const IS &is_p,
        const PDNSolution * const &pre_dot_disp,
        const PDNSolution * const &pre_dot_velo,
        const PDNSolution * const &pre_dot_pres,
        const PDNSolution * const &pre_disp,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_pres,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &feanode_ptr,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
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
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &dot_pres,
        PDNSolution * const &disp,
        PDNSolution * const &velo,
        PDNSolution * const &pres,
        bool &prestress_conv_flag, int &nl_counter ) const;
  
  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_thred;

    void Print_convergence_info( const int &count, const double &rel_err,
        const double &abs_err ) const
    {
      SYS_T::commPrint( "  === NR ite: %d, r_error: %e, a_error: %e \n", count, rel_err, abs_err);
    }

    // This function will add input * val to the solution vector in output, and
    // the ghost values in output will be updated.
    // User is responsible for making sure that the layouts of input and
    // output->solution are identical.
    void update_solid_kinematics( const double &val,
        const APart_Node * const &pnode,
        const Vec &input, 
        PDNSolution * const &output ) const;

    void rescale_inflow_value( const double &stime,
        ALocal_InflowBC * const &infbc,
        const ICVFlowRate * const &flrate,
        const PDNSolution * const &sol_base,
        PDNSolution * const &sol ) const;
};

#endif
