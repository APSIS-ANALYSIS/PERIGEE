#ifndef PGASSEM_TRANSPORT_NS_FEM_HPP
#define PGASSEM_TRANSPORT_NS_FEM_HPP
// ============================================================================
// PGAssem_Transport_NS_FEM.hpp
//
// Global assembly for a weakly coupled scalar transport equation in the NS
// examples.
//
// Date: Jul. 5 2026
// ============================================================================
#include "IPGAssem.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_InflowBC.hpp"
#include "PETSc_Tools.hpp"
#include "PLocAssem_Transport_VMS_NS_GenAlpha.hpp"

class PGAssem_Transport_NS_FEM : public IPGAssem
{
  public:
    PGAssem_Transport_NS_FEM(
        const ALocal_InflowBC * const in_inflow,
        std::unique_ptr<ALocal_EBC> in_ebc,
        std::unique_ptr<ALocal_IEN> in_locien,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode,
        const FEType &in_elemType,
        const int &in_nqp_s,
        std::unique_ptr<PLocAssem_Transport_VMS_NS_GenAlpha> in_locassem,
        const int &in_nz_estimate = 60 );

    virtual ~PGAssem_Transport_NS_FEM();

    using IPGAssem::Assem_residual;
    using IPGAssem::Assem_tangent_residual;

    virtual void Assem_nonzero_estimate();

    void Assem_mass_residual(
        const PDNSolution * const &phi,
        const PDNSolution * const &flow_sol );

    void Assem_residual(
        const PDNSolution * const &dot_phi,
        const PDNSolution * const &phi,
        const PDNSolution * const &flow_sol,
        const double &curr_time,
        const double &dt );

    void Assem_tangent_residual(
        const PDNSolution * const &dot_phi,
        const PDNSolution * const &phi,
        const PDNSolution * const &flow_sol,
        const double &curr_time,
        const double &dt );

    bool has_hi_model() const { return locassem->has_hi_model(); }

    double get_hi_beta() const { return locassem->get_hi_beta(); }

    double Assem_surface_HI_flowrate(
        const PDNSolution * const &transport_sol,
        const PDNSolution * const &flow_sol,
        const int &ebc_id ) const;

    double Assem_surface_HI_integral(
        const PDNSolution * const &transport_sol,
        const int &ebc_id ) const;

  private:
    const ALocal_InflowBC * const inflow;
    const std::unique_ptr<const ALocal_EBC> ebc;
    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const FEType elemType;
    const int nqps;
    const std::unique_ptr<FEAElement> elements;
    const std::unique_ptr<IQuadPts> quads;
    const std::unique_ptr<PLocAssem_Transport_VMS_NS_GenAlpha> locassem;
    const int nLocBas;
    const int snLocBas;
    const int dof_mat;
    const int nlgn;

    void Apply_inflow_BC_G();

    void Apply_inflow_BC_KG();

    void GetLocal(
        const double * const &array,
        const int * const &IEN,
        const int &in_dof,
        double * const &local_array ) const
    {
      for(int ii=0; ii<nLocBas; ++ii)
      {
        const int offset1 = ii * in_dof;
        const int offset2 = IEN[ii] * in_dof;
        for(int jj=0; jj<in_dof; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }
};

#endif
