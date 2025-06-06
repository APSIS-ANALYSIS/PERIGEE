#ifndef PLOCASSEM_VMS_NS_GENALPHA_INTERFACE_HPP
#define PLOCASSEM_VMS_NS_GENALPHA_INTERFACE_HPP
// ==================================================================
// PLocAssem_VMS_NS_GenAlpha_Interface.hpp
// 
// Parallel Local Assembly routine for VMS, Nitsche's method on 
// interfaces and Gen-alpha based NS solver.
//
// Author: Xuanming Huang
// Date: Jun 24 2024
// ==================================================================
#include "PLocAssem_VMS_NS_GenAlpha_WeakBC.hpp"
#include "Math_Tools.hpp"

class PLocAssem_VMS_NS_GenAlpha_Interface : public PLocAssem_VMS_NS_GenAlpha_WeakBC
{
  public:
    PLocAssem_VMS_NS_GenAlpha_Interface(
        const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s,
        const TimeMethod_GenAlpha * const &tm_gAlpha, const double &in_rho,
        const double &in_vis_mu, const double &in_beta,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0,
        const double &in_C_bI = 4.0 );

    virtual ~PLocAssem_VMS_NS_GenAlpha_Interface();

    virtual void print_info() const;

    virtual void Zero_Tangent_ss()
    {
      for(int ii{0}; ii<vec_size*vec_size; ++ii)
        Tangent_ss[ii] = 0.0;
    }

    virtual void Zero_Tangent_sr()
    {
      for(int ii{0}; ii<vec_size*vec_size; ++ii)
        Tangent_sr[ii] = 0.0;
    }

    virtual void Zero_Tangent_rr()
    {
      for(int ii{0}; ii<vec_size*vec_size; ++ii)
        Tangent_rr[ii] = 0.0;
    }

    virtual void Zero_Tangent_rs()
    {
      for(int ii{0}; ii<vec_size*vec_size; ++ii)
        Tangent_rs[ii] = 0.0;
    }

    virtual void Zero_Residual_s()
    {
      for(int ii{0}; ii<vec_size; ++ii)
        Residual_s[ii] = 0.0;
    }

    virtual void Zero_Residual_r()
    {
      for(int ii{0}; ii<vec_size; ++ii)
        Residual_r[ii] = 0.0;
    }

    virtual void Assem_Residual_itf_fixed(
      const int &fixed_qua,
      const double &fixed_qw,
      const double &dt,
      const FEAElement * const &fixed_elementv,
      const FEAElement * const &rotated_elementv,
      const double * const &fixed_local_sol, 
      const double * const &rotated_local_sol );

    virtual void Assem_Residual_itf_rotated(
      const int &rotated_qua,
      const double &rotated_qw,
      const double &dt,
      const FEAElement * const &rotated_elementv,
      const FEAElement * const &fixed_elementv,
      const double * const &rotated_local_sol,
      const double * const &rotated_local_mvelo, 
      const double * const &fixed_local_sol );

    virtual void Assem_Diag_Tangent_Residual_itf_fixed(
      const int &fixed_qua,
      const double &fixed_qw,
      const double &dt,
      const FEAElement * const &fixed_elementv,
      const FEAElement * const &rotated_elementv,
      const double * const &fixed_local_sol, 
      const double * const &rotated_local_sol );

    virtual void Assem_Diag_Tangent_Residual_itf_rotated(
      const int &rotated_qua,
      const double &rotated_qw,
      const double &dt,
      const FEAElement * const &rotated_elementv,
      const FEAElement * const &fixed_elementv,
      const double * const &rotated_local_sol,
      const double * const &rotated_local_mvelo, 
      const double * const &fixed_local_sol );

    virtual void Assem_Tangent_itf_MF_fixed(
      const int &fixed_qua,
      const double &fixed_qw,
      const double &dt,
      const FEAElement * const &fixed_elementv,
      const FEAElement * const &rotated_elementv,
      const double * const &fixed_local_sol,
      const double * const &rotated_local_sol );

    virtual void Assem_Tangent_itf_MF_rotated(
      const int &rotated_qua,
      const double &rotated_qw,
      const double &dt,
      const FEAElement * const &rotated_elementv,
      const FEAElement * const &fixed_elementv,
      const double * const &rotated_local_sol,
      const double * const &fixed_local_sol,
      const double * const &rotated_local_mvelo );
};

#endif
