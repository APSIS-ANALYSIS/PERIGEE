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
        const IGenBC * const &gbc,
        std::unique_ptr<ALocal_IEN> in_locien_v,
        std::unique_ptr<ALocal_IEN> in_locien_p,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode_v,
        std::unique_ptr<APart_Node> in_pnode_p,
        std::unique_ptr<ALocal_NBC> in_nbc_v,
        std::unique_ptr<ALocal_NBC> in_nbc_p,
        std::unique_ptr<ALocal_EBC> in_ebc_v,
        std::unique_ptr<ALocal_EBC> in_ebc_p,
        std::unique_ptr<IPLocAssem_2x2Block> in_locassem_f,
        std::unique_ptr<IPLocAssem_2x2Block> in_locassem_s,
        std::unique_ptr<Tissue_prestress> in_ps,  
        const int &in_nz_estimate=60 );

    virtual ~PGAssem_FSI();

    virtual void Assem_nonzero_estimate(
        const IGenBC * const &gbc );

    virtual void Assem_mass_residual(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres );

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
        const IGenBC * const &gbc );

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
        const IGenBC * const &gbc );

    // Assembly routine for the surface integrals for flow rates
    // and averaged pressure
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const int &ebc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        const int &ebc_id );

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id );

  private:
    const std::unique_ptr<const ALocal_IEN> locien_v;
    const std::unique_ptr<const ALocal_IEN> locien_p;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode_v;
    const std::unique_ptr<const APart_Node> pnode_p;
    const std::unique_ptr<const ALocal_NBC> nbc_v;
    const std::unique_ptr<const ALocal_NBC> nbc_p;
    const std::unique_ptr<const ALocal_EBC> ebc_v;
    const std::unique_ptr<const ALocal_EBC> ebc_p;
    const std::unique_ptr<IPLocAssem_2x2Block> locassem_f;
    const std::unique_ptr<IPLocAssem_2x2Block> locassem_s;
    const std::unique_ptr<const Tissue_prestress> ps;

    const int nLocBas, snLocBas, num_ebc, nlgn_v, nlgn_p;

    void EssBC_KG();

    void EssBC_G();

    void NatBC_G( const double &curr_time, const double &dt,
        const PDNSolution * const &disp );

    // Resistance BC for fluid
    void NatBC_Resis_G( const double &curr_time, const double &dt,
        const PDNSolution * const &disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &velo,
        const IGenBC * const &gbc );

    void NatBC_Resis_KG( const double &curr_time, const double &dt,
        const PDNSolution * const &disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &velo,
        const IGenBC * const &gbc );

    // Backflow stabilization
    void BackFlow_G( 
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo );

    void BackFlow_KG( const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo );

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
