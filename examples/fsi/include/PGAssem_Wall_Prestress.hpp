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
#include "ALocal_Elem.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_EBC.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_Wall_Prestress : public IPGAssem
{
  public:
    PGAssem_Wall_Prestress(
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
        std::unique_ptr<IPLocAssem_2x2Block> in_locassem_s,
        std::unique_ptr<Tissue_prestress> in_ps,
        const int &in_nz_estimate=60 );

    virtual ~PGAssem_Wall_Prestress();

    virtual void Assem_nonzero_estimate();

    virtual void Assem_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres );

    virtual void Assem_Tangent_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres );

    virtual void Update_Wall_Prestress(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres ) const;

    virtual void write_prestress_hdf5() const {ps->write_prestress_hdf5();}

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
    const std::unique_ptr<IPLocAssem_2x2Block> locassem_s;
    const std::unique_ptr<Tissue_prestress> ps;

    const int nLocBas, snLocBas, num_ebc, nlgn_v, nlgn_p;

    void EssBC_KG();

    void EssBC_G();

    void NatBC_G( const double &curr_time,
        const PDNSolution * const &pres );

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
