#ifndef PLOCASSEM_TRANSPORT_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_TRANSPORT_VMS_NS_GENALPHA_HPP
// ============================================================================
// PLocAssem_Transport_VMS_NS_GenAlpha.hpp
//
// Local assembly for a weakly coupled scalar transport equation using
// generalized-alpha in time and a SUPG/RBVMS-style stabilization in space.
//
// Date: Jul. 5 2026
// ============================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "FEAElementFactory.hpp"
#include "IHIModel.hpp"
#include "QuadPtsFactory.hpp"
#include "Math_Tools.hpp"
#include "SymmTensor2_3D.hpp"

class PLocAssem_Transport_VMS_NS_GenAlpha
{
  public:
    PLocAssem_Transport_VMS_NS_GenAlpha(
        const FEType &in_type, const int &in_nqp_v, const double &in_ct,
        const double &in_vis_mu,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        std::unique_ptr<IHIModel> in_hi_model );

    ~PLocAssem_Transport_VMS_NS_GenAlpha();

    int get_dof() const { return 1; }

    int get_nLocBas() const { return nLocBas; }

    void Zero_Tangent_Residual();

    void Zero_Residual();

    void Assem_Estimate();

    void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_phi,
        const double * const &phi,
        const double * const &flow_sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_phi,
        const double * const &phi,
        const double * const &flow_sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Mass_Residual(
        const double * const &phi,
        const double * const &flow_sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    bool has_hi_model() const { return hi_model != nullptr; }

    double get_hi_beta() const
    {
      return hi_model == nullptr ? 1.0 : hi_model->get_beta();
    }

    PetscScalar * Tangent;
    PetscScalar * Residual;

  private:
    const FEType elemType;
    const int nqpv;
    const std::unique_ptr<FEAElement> elementv;
    const std::unique_ptr<IQuadPts> quadv;
    const std::unique_ptr<IHIModel> hi_model;
    const double alpha_f, alpha_m, gamma;
    const double Ct;
    const double vis_mu;
    const int nLocBas;
    const int vec_size;

    void print_info() const;

    double get_tau(
        const double &dt,
        const std::array<double, 9> &dxi_dx,
        const double &ux, const double &uy, const double &uz ) const;

    double get_scalar_shear_stress(
        const double &u_x, const double &u_y, const double &u_z,
        const double &v_x, const double &v_y, const double &v_z,
        const double &w_x, const double &w_y, const double &w_z ) const;

    double get_source( const double &tau_s ) const
    {
      return hi_model == nullptr ? 0.0 : hi_model->get_source(tau_s);
    }
};

#endif
