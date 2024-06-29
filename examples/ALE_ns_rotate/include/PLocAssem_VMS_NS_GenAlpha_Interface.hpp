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
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const int &elemtype,
        const Vector_3 &point_xyz, const Vector_3 &angular,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0,
        const double &in_C_bI = 4.0 );

    virtual ~PLocAssem_VMS_NS_GenAlpha_Interface();

    PetscScalar * Tangent_ss;
    PetscScalar * Tangent_sr;
    PetscScalar * Tangent_rs;
    PetscScalar * Tangent_rr;

    PetscScalar * Residual_s;
    PetscScalar * Residual_r;

    virtual void print_info() const;

    virtual void Zero_Tangent_Residual_itf()
    {
      for(int ii{0}; ii<vec_size; ++ii)
      {
        Residual_s[ii] = 0.0;
        Residual_r[ii] = 0.0;
      }
      for(int ii{0}; ii<vec_size*vec_size; ++ii)
      {
        Tangent_ss[ii] = 0.0;
        Tangent_sr[ii] = 0.0;
        Tangent_rs[ii] = 0.0;
        Tangent_rr[ii] = 0.0;
      }
    }

    virtual void Zero_Residual_itf()
    {
      for(int ii{0}; ii<vec_size; ++ii)
      {
        Residual_s[ii] = 0.0;
        Residual_r[ii] = 0.0;
      }
    }

    virtual void Assem_Residual_itf(); // 积分点loop在此函数之外，传入具体某个积分点的、在两个区域中的basis，解和法向量等

    virtual void Assem_Tangent_Residual_itf(const double &dt); // 积分点loop在此函数之外，传入具体某个积分点的、在两个区域中的basis，解和法向量等
};

#endif