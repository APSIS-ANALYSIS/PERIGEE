#ifndef PGASSEM_HERK_BLOCK_NS_FEM_HPP
#define PGASSEM_HERK_BLOCK_NS_FEM_HPP
// ==================================================================
// PGAssem_HERK_Block_NS_FEM.hpp
//
// Parallel golbal assembly based on PETSc, using AIJ matrix format.
// The assembly routine is designed for classical C0 FEM method, which
// means we do not need extraction operators and local mesh sizes.
//
// The assembly is for the NS equations written in VMS and HERK formualtion
// of NS equations. The input solution vectors contains
//  [ pressure; velocity/dot velocity ].
//
// Author: Yujie Sun 
// Date Created: Mar. 4 2025
// ==================================================================
#include "PETSC_Tools.hpp"
#include "PDNSolution_NS.hpp"

class PGAssem_HERK_Block_NS_FEM
{
  public:
    Mat subK[5];

    Vec G;
    Vec subG[2];

    // Constructor
    PGAssem_HERK_Block_NS_FEM(
        std::unique_ptr<ALocal_IEN> in_locien,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode,
        std::unique_ptr<ALocal_NBC> in_nbc,
        std::unique_ptr<ALocal_EBC> in_ebc,
        std::unique_ptr<ALocal_WeakBC> in_wbc,
        std::unique_ptr<IPLocAssem_2x2Block> in_locassem,
        const int &in_nz_estimate=60 );

    // Destructor
    ~PGAssem_HERK_Block_NS_FEM();

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

    // ------------------------------------------------------------------------
    // ! Clear K and G to be zero
    // ------------------------------------------------------------------------
    void Clear_KG()
    {
      MatZeroEntries(K);
      VecSet(G, 0.0);
    }

    // ------------------------------------------------------------------------
    // ! Clear G to be zero
    // ------------------------------------------------------------------------
    void Clear_G() {VecSet(G, 0.0);}

    // Nonzero pattern estimate for the NS equations
    void Assem_nonzero_estimate();
    
    // Assembly mass matrix and residual vector
    void Assem_mass_residual(
        const PDNSolution * const &sol_a );

    // Assembly the residual vector and tangent matrix for the sub-step of HERK
    void Assem_tangent_residual_substep(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt );

    // Assembly the residual vector and tangent matrix for the final step of HERK
    void Assem_tangent_residual_finalstep(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt );

    // Assembly the residual vector and tangent matrix for the pres stage of HERK 
    void Assem_tangent_residual_presstage(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt );

  private:
    // Private data
    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const std::unique_ptr<const ALocal_NBC> nbc;
    const std::unique_ptr<const ALocal_EBC> ebc;
    const std::unique_ptr<IPLocAssem_2x2Block> locassem;

    const int nLocBas, snLocBas, dof_sol, dof_mat_v, dof_mat_p, num_ebc, nlgn;

    // Natural boundary condition for substep of the HERK
    void NatBC_G_HERK_Sub( const double &curr_time, const double &dt,
        const int &substep_index,
        const Runge_Kutta_Butcher * const &tm_RK_ptr );

    // Natural boundary condition for laststep of the HERK
    void NatBC_G_HERK_Final( const double &curr_time, const double &dt,
        const Runge_Kutta_Butcher * const &tm_RK_ptr );

    // Natural boundary condition for finalstep of the HERK
    void NatBC_G_HERK_Pres( const double &curr_time, const double &dt );

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
};

#endif
