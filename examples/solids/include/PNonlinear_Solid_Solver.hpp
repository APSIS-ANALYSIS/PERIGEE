#ifndef PNONLINEAR_SOLID_SOLVER_HPP
#define PNONLINEAR_SOLID_SOLVER_HPP
// ============================================================================
// PNonlinear_Solid_Solver.hpp
//
// Nonlinear solver for hyperelastic solid with mixed u-p formulation.
//
// Date: Jan 29 2026
// ==========================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution.hpp"
#include "APart_Node.hpp"

class PNonlinear_Solid_Solver
{
  public:
    PNonlinear_Solid_Solver(
        std::unique_ptr<IPGAssem> in_gassem,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
        std::unique_ptr<Matrix_PETSc> in_bc_mat,
        std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
        std::unique_ptr<APart_Node> in_pnode,
        const double &input_nrtol, const double &input_natol,
        const double &input_ndtol, const int &input_max_iteration,
        const int &input_renew_freq, const int &input_renew_threshold );

    ~PNonlinear_Solid_Solver() = default;

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    void print_lsolver_info() const {lsolver->print_info();}

    void GenAlpha_Seg_solve_Solid(
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
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &dot_pres,
        PDNSolution * const &disp,
        PDNSolution * const &velo,
        PDNSolution * const &pres,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    const std::unique_ptr<IPGAssem> gassem;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver;
    const std::unique_ptr<Matrix_PETSc> bc_mat;
    const std::unique_ptr<TimeMethod_GenAlpha> tmga;
    const std::unique_ptr<const APart_Node> pnode;

    void Print_convergence_info( const int &count,
        const double &rel_err, const double &abs_err ) const
    {
      SYS_T::commPrint("  === NR ite: %d, r_error: %e, a_error: %e \n", count, rel_err, abs_err);
    }

    void update_solid_kinematics( const double &val,
        const Vec &input,
        PDNSolution * const &output ) const;

    PNonlinear_Solid_Solver() = delete;
};

#endif
