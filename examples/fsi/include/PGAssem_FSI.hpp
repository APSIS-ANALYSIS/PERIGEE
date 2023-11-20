#ifndef PGASSEM_FSI_HPP
#define PGASSEM_FSI_HPP
// ============================================================================
// PGAssem_FSI.hpp
// 
// Parallel global assembly for FSI problems using the unified continuum
// formulation and segregated algorithm.
// 
// Author: Ju Liu
// Date: Jan 2 2022
// ============================================================================
#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"
#include "PDNSolution_V.hpp"

class PGAssem_FSI : public IPGAssem
{
  public:
    PGAssem_FSI( 
        IPLocAssem_2x2Block * const &locassem_f_ptr,
        IPLocAssem_2x2Block * const &locassem_s_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quads,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_v,
        const ALocal_IEN * const &aien_p,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        const ALocal_NBC * const &part_nbc_v,
        const ALocal_NBC * const &part_nbc_p,
        const ALocal_EBC * const &part_ebc,
        const IGenBC * const &gbc,
        const int &in_nz_estimate = 60 );

    virtual ~PGAssem_FSI();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const APart_Node * const &pnode_v,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    virtual void Assem_mass_residual(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
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
        const ALocal_EBC * const &ebc_part,
        const Tissue_prestress * const &ps_ptr,
        const Tissue_property * const &tp_ptr );

    virtual void Assem_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const PDNSolution * const &dot_velo_np1,
        const PDNSolution * const &velo_np1,
        const PDNSolution * const &disp_np1,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
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
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const Tissue_prestress * const &ps_ptr,
        const Tissue_property * const &tp_ptr );

    virtual void Assem_Tangent_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const PDNSolution * const &dot_velo_np1,
        const PDNSolution * const &velo_np1,
        const PDNSolution * const &disp_np1,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
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
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const Tissue_prestress * const &ps_ptr,
        const Tissue_property * const &tp_ptr );

    // Assembly routine for the surface integrals for flow rates
    // and averaged pressure
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const int &ebc_id );

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id );

  private:
    const int nLocBas, snLocBas, num_ebc, nlgn_v, nlgn_p;

    void EssBC_KG( const ALocal_NBC * const &nbc_v, const ALocal_NBC * const &nbc_p );

    void EssBC_G( const ALocal_NBC * const &nbc_v, const ALocal_NBC * const &nbc_p );

    void NatBC_G( const double &curr_time, const double &dt,
        const PDNSolution * const &disp,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_v,
        const ALocal_EBC * const &ebc_part );

    // Resistance BC for fluid
    void NatBC_Resis_G( const double &curr_time, const double &dt,
        const PDNSolution * const &disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_v,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    void NatBC_Resis_KG( const double &curr_time, const double &dt,
        const PDNSolution * const &disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_v,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    // Backflow stabilization
    void BackFlow_G( 
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_v,
        const ALocal_EBC * const &ebc_part );

    void BackFlow_KG( const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_v,
        const ALocal_EBC * const &ebc_part );

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
