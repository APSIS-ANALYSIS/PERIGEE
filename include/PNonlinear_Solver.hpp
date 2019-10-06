#ifndef PNONLINEAR_SOLVER_HPP
#define PNONLINEAR_SOLVER_HPP
// ==================================================================
// PNonlinear_Solver.hpp
//
// This class implements the nonlinear solver procedure.
//
// Date: Dec 9 2013
// ==================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "IALocal_BC.hpp"
#include "IPLocAssem.hpp"
#include "PGAssem.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"

class PNonlinear_Solver
{
  public:
    PNonlinear_Solver(
        const double &input_nrtol, const double &input_natol,
        const double &input_ndtol,
        const int &input_max_iteration, const int &input_renew_freq );

    ~PNonlinear_Solver();

    int get_non_max_its() const {return nmaxits;}

    void Info() const;

    //! Generalized-alpha nonlinear solver
    //  which is used in each time step of generalized-alpha time 
    //  scheme
    //  -------------------------------------------------------------
    //  \para bool new_tangent_flag: true->assembly tangent matrix
    //  \para double curr_time: t
    //  \para double dt: dt
    //  \para PDNSolution pre_disp: d_n
    //  \para PDNSolution pre_velo: v_n
    //  \para TimeMethod_GenAlpha tmga_ptr: pointer to time method
    //  \para ALocal_Elem alelem_ptr: local element info
    //  \para ALocal_IEN lien_ptr: IEN arrays
    //  \para APart_Node anode_ptr: local_to_global info
    //  \para FEANode feanode_ptr: control points
    //  \para IALocal_BC bc_part: boundary conditions
    //  \para AInt_Weight wei_ptr: quadrature weights
    //  \para vector<FEAElement*>: basis function at quadrature pts
    //  \para IPLocAssem lassem_ptr: interface for local assembly 
    //  \para PGAssem gassem_ptr: global assembly method
    //  \para PLinear_Solver_PETSc: linear solver method
    //  \para output PDNSolution disp: d_n+1
    //  \para output PDNSolution velo: v_n+1
    //  \para output bool conv_flag: true if nl iteration converged
    //  \para output int non_ite_counter: number of nl iterations 
    //  -------------------------------------------------------------
    void Gen_alpha_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const IALocal_BC * const &bc_part,
        const AInt_Weight * const &wei_ptr,
        const std::vector<FEAElement *> &ele_ptr,
        IPLocAssem * const &lassem_ptr,
        PGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &velo,
        PDNSolution * const &disp,
        bool &conv_flag,
        int &nl_counter ) const;


    // ! NewtonRaphson_solve
    //   The Newton-Raphson nonlinear solver, given the PDNSolution
    //   curr, which is the solution at time step n, the nonlinear
    //   solver gives the PDNSolution next, which is the solution at
    //   time step n+1. conv_flag is the bool which tells the
    //   convergence status of this solver, nl_counter gives the number
    //   of nonlinear iteration.
    //   ------------------------------------------------------------
    //  \para step: the PDNSolution class that is used as the nonlinear
    //              solver's incremental step. This should be allocated
    //              before this function call, and delete after this 
    //              function call.
    //   ------------------------------------------------------------
    void NewtonRaphson_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &curr,
        PDNSolution * const &step,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const IALocal_BC * const &bc_part,
        const AInt_Weight * const &wei_ptr,
        const std::vector<FEAElement *> &ele_ptr,
        IPLocAssem * const &lassem_ptr,
        PGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &next,
        bool &conv_flag,
        int &nl_counter ) const;


    // ! Gen_alpha_VMS_solve
    //   This is a time integration scheme like the generalized-alpha
    //   method, which is applied for the VMS-Navier-Stokes solver.
    //   The method can be found on pp 182, CMAME 197 2007.
    //   The key difference is that the pressure variable is not updated
    //   in the generalized alpha fashion.
    //   It is assumed that the last dof at each node is the pressure node.
    void Gen_alpha_VMS_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const IALocal_BC * const &bc_part,
        const AInt_Weight * const &wei_ptr,
        const std::vector<FEAElement *> &ele_ptr,
        IPLocAssem * const &lassem_ptr,
        PGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &velo,
        PDNSolution * const &disp,
        bool &conv_flag,
        int &nl_counter ) const;
  

    // ------------------------------------------------------------------------
    // ! Gen_alpha_VMS_noCache_3D_solve
    //   Gen_alpha_VMS_noCache_solve function is a variant of the original
    //   Gen_alpha_VMS_solve function. 
    //   This function does NOT require cached quadrature info, i.e. 
    //   vector<FEAElement *>. Instead, the following input needs to be
    //   passed into the function to calculate the quadrature info in each
    //   assembly call.
    //   1. IALocal_meshSize : object saving the mesh size
    //   2. BernsteinBasis_Array bs, bt, bu : pre-evaluated Berntesin
    //                         polynomial at Gauss quadrature points for
    //                         volumetric integration
    //   3. IAExtractor : Bezier extraction operator
    //   This function is written for 3D problems (input bs bt bu).
    //   See the comment of Gen_alpha_VMS_solve for more info.
    // ------------------------------------------------------------------------
    void Gen_alpha_VMS_noCache_3D_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const IALocal_BC * const &bc_part,
        const AInt_Weight * const &wei_ptr,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &velo,
        PDNSolution * const &disp,
        bool &conv_flag,
        int &nl_counter ) const;



    // ------------------------------------------------------------------------
    // ! NewtonRaphson_solve
    //   This is a variant of the NewtonRaphson_solve function for the no-cache
    //   style.
    //   1. We pass-in the element basis memory layout, extraction operator,
    //   and the Bezier element quadrature info for volumetric integration.
    //   2. The boundary condition is split into nodal_BC and elementa_BC.   
    // ------------------------------------------------------------------------
    void NewtonRaphson_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &curr,
        PDNSolution * const &step,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part,
        const AInt_Weight * const &wei_ptr,
        FEAElement * const &element,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const IALocal_meshSize * const &mSize,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &next,
        bool &conv_flag, int &nl_counter ) const;


    // ------------------------------------------------------------------------
    // ! NewtonRaphson_solve
    //   This is a variant of the NewtonRaphson_solve function for the no-cache
    //   style. In addition to the previous 2 points of modifications, this
    //   Newton-Raphson solver enables a thrid point:
    //   3. A Matrix_PETSc object is added as a parameter passed in the function
    //   call that enforces the strong satisfaction of essential boundary
    //   conditions.
    // ------------------------------------------------------------------------
    void NewtonRaphson_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &curr,
        PDNSolution * const &step,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        const AInt_Weight * const &wei_ptr,
        FEAElement * const &element,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const IALocal_meshSize * const &mSize,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &next,
        bool &conv_flag, int &nl_counter ) const;

    
    
    
    // --------------------------------------------------------------
    // This is a new Generalized-alpha method's nonlinear solver.
    // This solver is a no-cache style, which requires the input of
    // an element container (FEAElement), the volumetric Bernstein
    // basis function bs/bt/bu, the extraction operator suite--
    // extractor, and the mSize, containing the parametric element 
    // size.  A Matrix_PETSc::bc_mat is passed to enforce the strong
    // essential conditions.
    //
    // The design of this solver is for the second-order, elastodynamic
    // problems.
    //
    // Date: July 7 2016 
    // --------------------------------------------------------------
    void GenAlpha_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_acce,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        const AInt_Weight * const &wei_ptr,
        FEAElement * const &element,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const IALocal_meshSize * const &mSize,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &acce,
        PDNSolution * const &velo,
        PDNSolution * const &disp,
        bool &conv_flag, int &nl_counter ) const;


    // --------------------------------------------------------------
    // This is a Generalized-alpha solver.
    // This solver is no-cache type, which requires the input of an 
    // element container, the volumetric Bernstein basis functions, 
    // the extraction operator, and mSize. bc_mat is a sparse matrix 
    // used to enforce the essential boundary conditions.
    //
    // The design of this solver is for the first-order time derivative
    // type equations, like the Navier-Stokes.
    //
    // Date: Dec. 11 2016
    // --------------------------------------------------------------
    void GenAlpha_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        const AInt_Weight * const &wei_ptr,
        FEAElement * const &element,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const IALocal_meshSize * const &mSize,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &velo,
        PDNSolution * const &disp,
        bool &conv_flag, int &nl_counter ) const;


    // --------------------------------------------------------------
    // This is a Generalized-alpha solver.
    // This solve is no-cache type for classical finite element
    // methods (instead of IGA methods).
    //
    // The design of this solver is for first-order time derivative
    // type equations.
    //
    // Date: Jan. 24 2017
    // --------------------------------------------------------------
    void GenAlpha_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_disp,
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
        PDNSolution * const &velo,
        PDNSolution * const &disp,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol;
    const double na_tol;
    const double nd_tol;
    const int nmaxits;
    const int nrenew_freq;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {PetscPrintf(PETSC_COMM_WORLD,
        "  === NR ite: %d, r_error: %e, a_error: %e \n", 
        count, rel_err, abs_err);}

};

#endif
