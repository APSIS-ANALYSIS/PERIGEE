#ifndef PGASSEM_BLOCK_NS_FEM_HERK_HPP
#define PGASSEM_BLOCK_NS_FEM_HERK_HPP
// ==================================================================
// PGAssem_Block_NS_FEM_HERK.hpp
//
// Parallel golbal assembly based on PETSc, using AIJ matrix format.
// The assembly routine is designed for classical C0 FEM method, which
// means we do not need extraction operators and local mesh sizes.
//
// The assembly is for the NS equations written in VMS and HERK formualtion
// of NS equations. The input solution vectors contains
//  [ velocity/dot velocity; pressure ].
//
// Author: Yujie Sun 
// Date Created: Mar. 4 2025
// ==================================================================
#include "PETSc_Tools.hpp"
#include "PDNSolution_NS.hpp"
#include "PDNSolution_V.hpp"
#include "PDNSolution_P.hpp"
#include "PLocAssem_Block_VMS_NS_HERK.hpp"
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"
#include "AGlobal_Mesh_Info.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_IEN.hpp"

class PGAssem_Block_NS_FEM_HERK
{
  public:
    Mat subK[6];

    Vec G;
    Vec subG[2];

    // Constructor
    PGAssem_Block_NS_FEM_HERK(
        std::unique_ptr<ALocal_IEN> in_locien,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode,
        std::unique_ptr<ALocal_NBC> in_nbc,
        std::unique_ptr<ALocal_EBC> in_ebc,
        std::unique_ptr<PLocAssem_Block_VMS_NS_HERK> in_locassem,
        const int &in_nz_estimate=60 );

    // Destructor
    ~PGAssem_Block_NS_FEM_HERK();

    // ------------------------------------------------------------------------
    // ! Flag : Fix nonzero structure
    //          Add or insert in new location is ignored. Set after
    //         the first MatAssemblyEnd().
    // ------------------------------------------------------------------------
    void Fix_nonzero_str()
    {
      MatSetOption(subK[0], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(subK[1], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(subK[2], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(subK[3], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(subK[4], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    }

    // ------------------------------------------------------------------------
    // ! Flag : New allocation error
    //          Add or insert in new locations will generate an error message.
    //          Supports AIJ and BAIJ formats.
    // ------------------------------------------------------------------------
    void Fix_nonzero_err_str()
    {
      MatSetOption(subK[0], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(subK[1], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(subK[2], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(subK[3], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(subK[4], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    }
    
    // ------------------------------------------------------------------------
    // ! Flag : Ignore new allocation
    //          Add or insert in a new allocation will NOT generate error.
    // ------------------------------------------------------------------------
    void Release_nonzero_err_str()
    {
      MatSetOption(subK[0], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(subK[1], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(subK[2], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(subK[3], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(subK[4], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    // ------------------------------------------------------------------------
    // ! Flag : Keep nonzero pattern of the matrix K
    // ------------------------------------------------------------------------
    void Keep_nonzero_pattern()
    {
      MatSetOption(subK[0], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(subK[1], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(subK[2], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(subK[3], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(subK[4], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }

    void Clear_G() {VecSet(G, 0.0);}

    void Clear_subKG()
    {
      MatZeroEntries(subK[0]);
      MatZeroEntries(subK[1]);
      MatZeroEntries(subK[2]);
      MatZeroEntries(subK[3]);
      MatZeroEntries(subK[4]);
      VecSet(G, 0.0);
    }

    void Clear_subK5()
    {
      MatZeroEntries(subK[5]);
    }
    
    // Nonzero pattern estimate for the NS equations
    void Assem_nonzero_estimate();
    
    // Assembly the residual vector and tangent matrix for the sub-step of HERK
    void Assem_tangent_matrix(
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double &curr_time,
        const double &dt );

    void Assem_residual_substep(
        const int &substep_index,
        PDNSolution ** const &cur_velo_sols,
        PDNSolution ** const &cur_pres_sols,
        PDNSolution ** const &pre_velo_sols,
        PDNSolution * const &pre_velo,
        PDNSolution ** const &pre_pres_sols,
        PDNSolution * const &pre_velo_before,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double &curr_time,
        const double &dt );

    void Assem_residual_finalstep(
        PDNSolution ** const &cur_velo_sols,
        PDNSolution * const &cur_velo,
        PDNSolution ** const &cur_pres_sols,
        PDNSolution ** const &pre_velo_sols,
        PDNSolution * const &pre_velo,
        PDNSolution ** const &pre_pres_sols,
        PDNSolution * const &pre_velo_before,    
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double &curr_time,
        const double &dt );

    void Assem_residual_presstage(
        PDNSolution * const &cur_dot_velo,
        PDNSolution ** const &cur_velo_sols,
        PDNSolution * const &cur_velo,
        PDNSolution ** const &cur_pres_sols,
        PDNSolution * const &pre_velo,
        PDNSolution * const &cur_pres,    
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double &curr_time,
        const double &dt );

    // Assembly the residual vector and tangent matrix for the sub-step of HERK
    void Assem_tangent_residual_substep(
        const int &substep_index,
        PDNSolution ** const &cur_velo_sols,
        PDNSolution ** const &cur_pres_sols,
        PDNSolution ** const &pre_velo_sols,
        PDNSolution * const &pre_velo,
        PDNSolution ** const &pre_pres_sols,
        PDNSolution * const &pre_velo_before,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double &curr_time,
        const double &dt );

    // Assembly the residual vector and tangent matrix for the final step of HERK
    void Assem_tangent_residual_finalstep(
        PDNSolution ** const &cur_velo_sols,
        PDNSolution * const &cur_velo,
        PDNSolution ** const &cur_pres_sols,
        PDNSolution ** const &pre_velo_sols,
        PDNSolution * const &pre_velo,
        PDNSolution ** const &pre_pres_sols,
        PDNSolution * const &pre_velo_before,    
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double &curr_time,
        const double &dt );

    // Assembly the residual vector and tangent matrix for the pres stage of HERK 
    void Assem_tangent_residual_presstage(
        PDNSolution * const &cur_dot_velo,
        PDNSolution ** const &cur_velo_sols,
        PDNSolution * const &cur_velo,
        PDNSolution ** const &cur_pres_sols,
        PDNSolution * const &pre_velo,
        PDNSolution * const &cur_pres,    
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double &curr_time,
        const double &dt );

    void Update_tangent_alpha_RK( const double &aa )
    {tangent_alpha_RK = aa;}

    double Get_tangent_alpha_RK()
    {return tangent_alpha_RK;}

    void Update_tangent_submatrix5()
    {
      // if (subK[5] != NULL) MatDestroy(&subK[5]); // Prevent duplicate allocation leakage
      // MatDuplicate(subK[3], MAT_COPY_VALUES, &subK[5]);
      // MatAXPY(subK[5], tangent_alpha_RK, subK[4], DIFFERENT_NONZERO_PATTERN);

      if (subK[5] == NULL) // If subK [5] has not yet allocated memory, proceed with memory allocation
        MatDuplicate(subK[3], MAT_COPY_VALUES, &subK[5]);
      else // If subK [5] already has memory, update its contents directly
        MatCopy(subK[3], subK[5], DIFFERENT_NONZERO_PATTERN);

      MatAXPY(subK[5], tangent_alpha_RK, subK[4], DIFFERENT_NONZERO_PATTERN);
    }

  private:
    // Private data
    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const std::unique_ptr<const ALocal_NBC> nbc;
    const std::unique_ptr<const ALocal_EBC> ebc;
    const std::unique_ptr<PLocAssem_Block_VMS_NS_HERK> locassem;

    const int nLocBas, snLocBas, dof_sol, dof_mat_v, dof_mat_p, num_ebc, nlgn;

    // Runge Kuta butcher table coefficients
    double tangent_alpha_RK;

    // Essential boundary condition
    void EssBC_KG();

    void EssBC_K();

    void EssBC_G();

    // Natural boundary condition for substep of the HERK
    void NatBC_G_HERK_Sub( const double &curr_time, const double &dt,
        const int &substep_index,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr );

    // Natural boundary condition for laststep of the HERK
    void NatBC_G_HERK_Final( const double &curr_time, const double &dt,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr );

    // Natural boundary condition for finalstep of the HERK
    void NatBC_G_HERK_Pressure( const double &curr_time, const double &dt );

    void GetLocal(const double * const &array, const int * const &IEN,
        double * const &local_array) const
    {
      for(int ii=0; ii<nLocBas; ++ii)
      {
        const int offset1 = ii * dof_sol;
        const int offset2 = IEN[ii] * dof_sol;
        for(int jj=0; jj<dof_sol; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, double * const &local_array) const
    {
      for(int ii=0; ii<in_locbas; ++ii)
      {
        const int offset1 = ii * dof_sol;
        const int offset2 = IEN[ii] * dof_sol;
        for(int jj=0; jj<dof_sol; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    std::vector<double> GetLocal( const std::vector<double> &array,
      const std::vector<int> &IEN, const int &in_locbas, const int &in_dof ) const
    {
      std::vector<double> out( in_locbas * in_dof, 0.0 );
      for(int ii=0; ii<in_locbas; ++ii)
        for(int jj=0; jj<in_dof; ++jj)
          out[ii * in_dof + jj] = array[IEN[ii] * in_dof + jj];
  
      return out;
    }
};

#endif
