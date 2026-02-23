#ifndef PGASSEM_SOLID_FEM_HPP
#define PGASSEM_SOLID_FEM_HPP
// ============================================================================
// PGAssem_Solid_FEM.hpp
//
// Parallel global assembly for hyperelastic solid with mixed u-p formulation.
//
// Date: Jan 29 2026
// ==========================================================================
#include "IPGAssem.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_NBC_Solid.hpp"
#include "ALocal_EBC.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_Solid_FEM : public IPGAssem
{
  public:
    PGAssem_Solid_FEM(
        std::unique_ptr<ALocal_IEN> in_locien,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode,
        std::unique_ptr<ALocal_NBC_Solid> in_nbc,
        std::unique_ptr<ALocal_EBC> in_ebc,
        std::unique_ptr<IPLocAssem_2x2Block> in_locassem,
        const int &in_nz_estimate = 60 );

    virtual ~PGAssem_Solid_FEM();

    virtual void Assem_nonzero_estimate() override;

    virtual void Assem_mass_residual(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres ) override;

    virtual void Assem_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres ) override;

    virtual void Assem_Tangent_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres ) override;

    void Apply_Dirichlet_BC(
        const double &time,
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &disp,
        PDNSolution * const &velo ) const;

    const ALocal_NBC_Solid * get_nbc() const { return nbc.get(); }

  private:
    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const std::unique_ptr<const ALocal_NBC_Solid> nbc;
    const std::unique_ptr<const ALocal_EBC> ebc;
    const std::unique_ptr<IPLocAssem_2x2Block> locassem;

    const int num_ebc;
    const int nLocBas, snLocBas, dof_mat, nlgn, nqpv;

    std::vector<double> zero_prestress;

    void EssBC_KG( const int &field );
    void EssBC_G( const int &field );

    void NatBC_G( const double &curr_time, const double &dt );

    std::vector<double> GetLocal( const std::vector<double> &array,
        const std::vector<int> &IEN, const int &in_locbas, const int &in_dof ) const
    {
      std::vector<double> out( in_locbas * in_dof, 0.0 );
      for(int ii=0; ii<in_locbas; ++ii)
        for(int jj=0; jj<in_dof; ++jj)
          out[ii * in_dof + jj] = array[ IEN[ii] * in_dof + jj ];

      return out;
    }

    PGAssem_Solid_FEM() = delete;
};

#endif
