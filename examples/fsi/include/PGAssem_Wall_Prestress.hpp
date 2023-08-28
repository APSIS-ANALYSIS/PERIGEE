#ifndef PGASSEM_WALL_PRESTRESS_HPP
#define PGASSEM_WALL_PRESTRESS_HPP
// ============================================================================
// PGAssem_Wall_Prestress.hpp
// 
// Parallel Global Assembly for wall prestress generation.
//
// Date: Jan 30 2022
// ============================================================================
#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_Wall_Prestress : public IPGAssem
{
  public:
    PGAssem_Wall_Prestress( 
        IPLocAssem_2x2Block * const &locassem_s_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_v,
        const ALocal_IEN * const &aien_p,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        const ALocal_NBC * const &part_nbc_v,
        const ALocal_NBC * const &part_nbc_p,
        const ALocal_EBC * const &part_ebc,
        const int &in_nz_estimate = 60 );

    virtual ~PGAssem_Wall_Prestress();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p );

    virtual void Assem_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const Tissue_prestress * const &ps_ptr );

    virtual void Assem_Tangent_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const Tissue_prestress * const &ps_ptr );

    virtual void Update_Wall_Prestress(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        const IQuadPts * const &quadv,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        Tissue_prestress * const &ps_ptr ) const;

  private:
    const int nLocBas, snLocBas, num_ebc, nlgn_v, nlgn_p;

    void EssBC_KG( const ALocal_NBC * const &nbc_v, const ALocal_NBC * const &nbc_p );

    void EssBC_G( const ALocal_NBC * const &nbc_v, const ALocal_NBC * const &nbc_p );

    void NatBC_G( const double &curr_time,
        const PDNSolution * const &pres,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_v,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p );

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
