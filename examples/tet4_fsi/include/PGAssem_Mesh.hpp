#ifndef PGASSEM_MESH_HPP
#define PGASSEM_MESH_HPP
// ============================================================================
// PGAssem_Mesh.hpp
//
// Parallel Global Assembly for Mesh motion equations.
//
// Author: Ju Liu
// Date: Dec. 30 2021
// ============================================================================
#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_Mesh : public IPGAssem
{
  public:
    PGAssem_Mesh( IPLocAssem * const &locassem_ptr,
        ALocal_Elem const * const &alelem_ptr,
        ALocal_IEN const * const &aien_ptr,
        APart_Node const * const &pnode_ptr,
        ALocal_NBC const * const &part_nbc,
        ALocal_EBC const * const &part_ebc,
        const int &in_nz_estimate );

    virtual ~PGAssem_Mesh();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NBC * const &nbc_part );

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

  private:
    const int nLocBas, snLocBas, dof, num_ebc, nlgn;
    
    void EssBC_KG( const ALocal_NBC * const &nbc_part );

    void EssBC_G(  const ALocal_NBC * const &nbc_part );

    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    std::vector<double> GetLocal( const std::vector<double> &array,
        const std::vector<int> &IEN, const int &in_locbas ) const
    {
      std::vector<double> out( in_locbas * dof, 0.0 );
      for(int ii=0; ii<in_locbas; ++ii)
        for(int jj=0; jj<dof; ++jj)
          out[ii*dof + jj] = array[ IEN[ii] * dof + jj ];
      
      return out;
    }

};

#endif
