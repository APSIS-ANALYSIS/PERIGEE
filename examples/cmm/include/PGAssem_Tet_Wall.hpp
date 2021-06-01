#ifndef PGASSEM_TET_WALL_HPP
#define PGASSEM_TET_WALL_HPP
// ============================================================================
// PGAssem_Tet_Wall.hpp
//
// Parallel global assembly of wall mechanics. The input solution vector
// contains [disp], [dot_disp], [velo], and [dot_velo].
//
// ============================================================================
#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"
#include "PDNSolution_NS.hpp"
#include "PDNSolution_Wall_Disp.hpp"

class PGAssem_Tet_Wall : public IPGAssem
{
  public:
    PGAssem_Tet_Wall( 
        const IPLocAssem * const &locassem_ptr,
        const int &in_nlocbas,
        const int &in_snlocbas,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &part_nbc,
        const int &in_nz_estimate );

    virtual ~PGAssem_Tet_Wall();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NodalBC * const &nbc_part );

  private:
    const int nLocBas, dof_sol, dof_mat, nlgn, snLocBas;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const ALocal_NodalBC * const &nbc_part, const int &field );

    void EssBC_G( const ALocal_NodalBC * const &nbc_part, const int &field );

};

#endif
