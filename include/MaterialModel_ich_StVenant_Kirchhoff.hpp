#ifndef MATERIALMODEL_ICH_STVENANT_KIRCHHOF_HPP
#define MATERIALMODEL_ICH_STVENANT_KIRCHHOF_HPP
// ==================================================================
// MaterialModel_ich_StVenant_Kirchhoff.hpp
//
// Date: Aug. 14 2024
// Author: Ju Liu
// ==================================================================
#include "IMaterialModel_ich.hpp"

class MaterialModel_ich_StVenant_Kirchhoff : public IMaterialModel_ich
{
  public:
    MaterialModel_ich_StVenant_Kirchhoff( const double &in_mu ) : mu( in_mu ) {};

    virtual ~MaterialModel_ich_StVenant_Kirchhoff() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_ich_StVenant_Kirchhoff: \n");
      SYS_T::commPrint("\t  Shear modulus mu = %e \n", mu);
    }

    virtual std::string get_model_name() const {return std::string("StVenant_Kirchhoff");}

    virtual double get_elastic_mu() const {return mu;}

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {
      const double detFm0d67 = std::pow( F.det(), -2.0 / 3.0 );

      const auto C = STen2::gen_right_Cauchy_Green( F );

    	const auto S_tilde = mu * ( detFm0d67 * C - STen2::gen_id() );

      return detFm0d67 * STen2::gen_DEV_part( S_tilde, C );
    }

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_ich ) const
    {
      const auto S_ich = get_PK_2nd( F );

      // First PK stress
      P_ich = F * S_ich;

      constexpr double pt67 = 2.0 / 3.0;
      const double detFm0d67 = std::pow( F.det(), -pt67 );

      const auto C = STen2::gen_right_Cauchy_Green( F );

      auto CC_tilde = 2.0 * std::pow( F.det(), -2.0*pt67 ) * mu * STen4::gen_symm_id();

      Tensor4_3D PP = Ten4::gen_P( C );

      CC_tilde.TenPMult( PP );

      // Elasticity tensor
      auto CC_ich = CC_tilde;

      const auto PP_tilde = STen4::gen_Ptilde( STen2::inverse(C) );

      const auto S_tilde = mu * ( detFm0d67 * C - STen2::gen_id() );

      CC_ich += pt67 * detFm0d67 * S_tilde.MatContraction( C ) * PP_tilde;

      CC_ich.add_SymmOutProduct( -pt67, STen2::inverse(C), S_ich );

      return CC_ich;
    }

    virtual double get_energy( const Tensor2_3D &F ) const
    {
      const auto C_tilde = std::pow( F.det(), -2.0 / 3.0 ) * STen2::gen_right_Cauchy_Green( F );

      const auto E_tilde = 0.5 * ( C_tilde - STen2::gen_id() );

      return mu * E_tilde.MatContraction( E_tilde );
    }

    virtual Vector_3 get_fibre_dir( const int &dir ) const
    {
      return Vector_3();
    }

  private:
    const double mu;
};

#endif
