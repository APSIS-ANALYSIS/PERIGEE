#ifndef PLOCASSEM_2x2BLOCK_TET4_ALE_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_2x2BLOCK_TET4_ALE_VMS_NS_GENALPHA_HPP
// ============================================================================
// PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha.hpp
//
// This is a local assembly routine for ALE-VMS formulation of the 3D
// Navier-Stokes equations with Generalized-alpha for time stepping.
//
// Date: Jan 1 2022
// Author: Ju Liu
// ============================================================================
#include "IPLocAssem_2x2Block.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha : public IPLocAssem_2x2Block
{
  public:
    PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nqp,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta );

    virtual ~PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha();

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
      for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;
    
      for(int ii=0; ii<vec_size_0 * vec_size_0; ++ii) Tangent00[ii] = 0.0;
      for(int ii=0; ii<vec_size_0 * vec_size_1; ++ii) Tangent01[ii] = 0.0;
      for(int ii=0; ii<vec_size_1 * vec_size_0; ++ii) Tangent10[ii] = 0.0;
      for(int ii=0; ii<vec_size_1 * vec_size_1; ++ii) Tangent11[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size_0; ++ii) sur_Residual0[ii] = 0.0;
      for(int ii=0; ii<sur_size_1; ++ii) sur_Residual1[ii] = 0.0;

      for(int ii=0; ii<sur_size_0 * sur_size_0; ++ii) sur_Tangent00[ii] = 0.0;
      for(int ii=0; ii<sur_size_0 * sur_size_1; ++ii) sur_Tangent01[ii] = 0.0;
      for(int ii=0; ii<sur_size_1 * sur_size_0; ++ii) sur_Tangent10[ii] = 0.0;
      for(int ii=0; ii<sur_size_1 * sur_size_1; ++ii) sur_Tangent11[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
      for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;
    }

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size_0; ++ii) sur_Residual0[ii] = 0.0;
      for(int ii=0; ii<sur_size_1; ++ii) sur_Residual1[ii] = 0.0;
    }


  private:
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta, CI, CT;

    const int nLocBas_v, nLocBas_p, snLocBas_v, snLocBas_p;
    const int vec_size_0, vec_size_1, sur_size_0, sur_size_1, nqp;

    void print_info() const;
};

#endif
