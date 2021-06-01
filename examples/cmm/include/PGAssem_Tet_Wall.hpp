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
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &part_nbc,
        const ALocal_EBC * const &part_ebc,
        const int &in_nz_estimate );

    virtual ~PGAssem_Tet_Wall();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NodalBC * const &nbc_part );

    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &sol_wall_disp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementw,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_Ring_NodalBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_wall_part,
        const IGenBC * const &gbc );

  private:
    const int nLocBas, dof_sol, dof_mat, nlgn, snLocBas;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const ALocal_NodalBC * const &nbc_part, const int &field );

    void EssBC_G( const ALocal_NodalBC * const &nbc_part, const int &field );

    void WallMembrane_KG( const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &sol_wall_disp,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_w,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_wall_part );

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

    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, const int &in_dof,
        double * const &local_array) const
    {
      for(int ii=0; ii<in_locbas; ++ii)
      {
        const int offset1 = ii * in_dof;
        const int offset2 = IEN[ii] * in_dof;
        for(int jj=0; jj<in_dof; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

};

#endif
