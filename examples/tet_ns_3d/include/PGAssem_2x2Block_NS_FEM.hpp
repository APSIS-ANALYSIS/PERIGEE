#ifndef PGASSEM_2X2BLOCK_NS_FEM_HPP
#define PGASSEM_2X2BLOCK_NS_FEM_HPP
// ==================================================================
// PGAssem_2x2Block_NS_FEM.hpp
//
// Parallel global assembly using nest matrix.
// The assembly is designed for C0 Lagrangian elements, which does
// not need extraction operator or local mesh sizes.
//
// Author: Ju Liu
// Date Created: Aug. 11 2020
// ==================================================================
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"
#include "IAGlobal_Mesh_Info.hpp"
#include "IPLocAssem_2x2Block.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"
#include "ALocal_NodalBC.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_EBC.hpp"
#include "IGenBC.hpp"

class PGAssem_2x2Block_NS_FEM
{
  public:
    Mat K;
    Mat subK[4];
    
    Vec G;
    Vec subG[2];

    PGAssem_2x2Block_NS_FEM(
        IPLocAssem_2x2Block * const &locassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quads,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &part_nbc,
        const ALocal_EBC * const &part_ebc,
        const IGenBC * const &gbc,
        const int &in_nz_estimate=60 );

    virtual ~PGAssem_2x2Block_NS_FEM();

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

    
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

  private:
    // --------------------------------------------------------------
    // Data:
    const int nLocBas, dof_sol, dof_mat_v, dof_mat_p, num_ebc;

    int snLocBas;

    PetscInt * row_index_v, * row_index_p, * srow_index_v, * srow_index_p;

    double * array_v, * array_p, * array_dot_v, * array_dot_p;

    double * local_v, * local_p, * local_dot_v, * local_dot_p;

    double * local_vs, * local_ps, * local_dot_vs, * local_dot_ps;

    int * IEN_e, * LSIEN;

    double * ectrl_x, * ectrl_y, * ectrl_z;

    double * sctrl_x, * sctrl_y, * sctrl_z;

    // --------------------------------------------------------------
    // Functions:
    void EssBC_KG_v( const ALocal_NodalBC * const &nbc_part );

    void EssBC_G_v( const ALocal_NodalBC * const &nbc_part );

};

#endif
