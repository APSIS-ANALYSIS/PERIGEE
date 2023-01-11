#ifndef PTIME_SOLVER_HPP
#define PTIME_SOLVER_HPP
// ==================================================================
// PTime_Solver.hpp
// 
// Parallel Time marching solver. We implement various different time
// integration schemes in this class.
//
// Date: Dec 9th 2013
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_Solver.hpp"

class PTime_Solver
{
  public:
    PTime_Solver(const std::string &input_name, const int &input_record_freq,
        const int &input_renew_tang_freq, const double &input_final_time );

    ~PTime_Solver();

    // ! print the key parameters of the time solver on screen
    void Info() const;

    // ! Perform generalized-alpha time marching
    void TM_generalized_alpha(
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
        PDNTimeStep * const &time_info,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;
    
    
    // ! Perform time marching using Newton-Raphson nonlinear solver
    void TM_NewtonRaphson(
        const PDNSolution * const &init_disp,
        PDNTimeStep * const &time_info,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;

    // ! Perform the VMS generalized-alpha time marching scheme
    //   See CMAME 197 (2007), pp 182
    void TM_VMS_generalized_alpha(
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
        PDNTimeStep * const &time_info,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;
  

    // ------------------------------------------------------------------------
    // ! Perform the VMS generalized-alpha time marching scheme
    //   with NO cached quadrature info for 3D problems.
    //   
    //   TM_VMS_generalized_alpha_noCache_3D is a variant of the original
    //   TM_VMS_generalized_alpha function.
    //
    //   This function does NOT require vector<FEAElement *> as input. Instead,
    //   the following data structures are needed for evaluating basis functions
    //   in each local assembly call:
    //
    //   1. IALocal_meshSize : object saving the mesh size
    //   2. BernsteinBasis_Array bs, bt, bu : pre-evaluated Berntesin
    //      polynomial at Gauss quadrature points for volumetric integration
    //   3. IAExtractor : Bezier extraction operator
    //   This function is written for 3D problems (input bs bt bu). 
    //   
    //   See CMAME 197 (2007), pp 182 for algorithm details.
    // ------------------------------------------------------------------------
    void TM_VMS_generalized_alpha_noCache_3D(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
        PDNTimeStep * const &time_info,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;

    
    // ------------------------------------------------------------------------
    // ! Time marching using Newton-Raphson nonlinear solver
    //   quadrature info is calculated inside the local assembly routine
    // ------------------------------------------------------------------------
    void TM_NewtonRaphson(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_disp,
        PDNTimeStep * const &time_info,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;

    
    // ------------------------------------------------------------------------
    // ! Time marching using Newton-Raphson nonlinear solver
    //   quadrature info is calculated inside the local assembly routine
    //
    //   A Matrix_PETSc object is passed into the solver with the objective of
    //   strongly enforcing the essential boundary conditions.
    // ------------------------------------------------------------------------
    void TM_NewtonRaphson(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_disp,
        PDNTimeStep * const &time_info,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;


    void TM_GenAlpha(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_acce,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;


    // --------------------------------------------------------------
    // This is the generalized-alpha method for first-order system.
    // --------------------------------------------------------------
    void TM_GenAlpha(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        PDNTimeStep * const &time_info,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;


    // --------------------------------------------------------------
    // This is the generalized-alpha method for first-order system,
    // with classical FEM method.
    //
    // Date: Jan 24 2017
    // --------------------------------------------------------------
    void TM_GenAlpha(
        const bool &restart_init_assembly_flag,
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
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
        PNonlinear_Solver * const &nsolver_ptr
        ) const;


  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const int &counter ) const;
};

#endif
