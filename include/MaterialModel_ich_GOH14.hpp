#ifndef MATERIALMODEL_ICH_GOH14_HPP
#define MATERIALMODEL_ICH_GOH14_HPP
// ============================================================================
// MaterialModel_ich_GOH14.hpp
//
// The isochoric isotropic part is the same as the GOH06_Incompressible. The 
// anisotropic part is modified by replacing the modified right Cauchy-Green 
// tensor by the right Cauchy-Green tensor.
//
// Ref: Nolan et al. J Mech Behav Biomed Mater 39 (2014) 48-60
//
// Date: Aug. 21 2024
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ============================================================================
#include "IMaterialModel_ich.hpp"
#include "Math_Tools.hpp"

class MaterialModel_ich_GOH14 : public IMaterialModel_ich
{
  public:
    MaterialModel_ich_GOH14( const double &in_mu,
    const double &in_f1the, const double &in_f1phi,
    const double &in_f2the, const double &in_f2phi,
    const double &in_fk1, const double &in_fk2,
    const double &in_fkd)
    : mu( in_mu ), f1_the(in_f1the* MATH_T::PI / 180.0),
    f1_phi(in_f1phi*MATH_T::PI / 180.0), f2_the(in_f2the * MATH_T::PI / 180.0),
    f2_phi(in_f2phi*MATH_T::PI / 180.0), fk1(in_fk1), fk2(in_fk2), fkd(in_fkd),
    a1( sin(f1_the) * cos(f1_phi), sin(f1_the) * sin(f1_phi), cos(f1_the) ),
    a2( sin(f2_the) * cos(f2_phi), sin(f2_the) * sin(f2_phi), cos(f2_the) ) {}

    virtual ~MaterialModel_ich_GOH14() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_ich_GOH14: \n");
      SYS_T::commPrint("\t  Shear modulus mu           = %e \n", mu);
      SYS_T::commPrint("\t  Dispersion parameter kappa = %e \n", fkd);
      SYS_T::commPrint("\t  Stress-like parameter k1   = %e \n", fk1);
      SYS_T::commPrint("\t  Dimensionless parameter k2 = %e \n", fk2);
      SYS_T::commPrint("\t  Mean angle theta_1 (deg)   = %e \n", f1_the*180/static_cast<double>(MATH_T::PI));
      SYS_T::commPrint("\t  Mean angle phi_1   (deg)   = %e \n", f1_phi*180/static_cast<double>(MATH_T::PI));
      SYS_T::commPrint("\t  Mean angle theta_1 (rad)   = %e \n", f1_the);
      SYS_T::commPrint("\t  Mean angle phi_1   (rad)   = %e \n", f1_phi);
      SYS_T::commPrint("\t  Mean angle theta_2 (deg)   = %e \n", f2_the*180/static_cast<double>(MATH_T::PI));
      SYS_T::commPrint("\t  Mean angle phi_2   (deg)   = %e \n", f2_phi*180/static_cast<double>(MATH_T::PI));
      SYS_T::commPrint("\t  Mean angle theta_2 (rad)   = %e \n", f2_the);
      SYS_T::commPrint("\t  Mean angle phi_2   (rad)   = %e \n", f2_phi);
      SYS_T::commPrint("\t  Mean direction a1 of the 1st family of fibres = [%e, %e, %e] \n", a1(0), a1(1), a1(2));
      SYS_T::commPrint("\t  Mean direction a2 of the 2nd family of fibres = [%e, %e, %e] \n", a2(0), a2(1), a2(2));
    }

    virtual std::string get_model_name() const {return std::string("GOH14");}

    virtual double get_elastic_mu() const {return mu;}

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {
      const auto C = STen2::gen_right_Cauchy_Green( F );

      const double I1 = C.tr();
      const double I4 = C.VecMatVec( a1, a1 );
      const double I6 = C.VecMatVec( a2, a2 );

      const double fE1 = fkd * I1 + ( 1.0 - 3.0 * fkd ) * I4 - 1.0;
      const double fE2 = fkd * I1 + ( 1.0 - 3.0 * fkd ) * I6 - 1.0;

      const double dfpsi1_dfE1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
      const double dfpsi2_dfE2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

      const auto H_f1 = fkd * STen2::gen_id() + ( 1.0 - 3.0 * fkd ) * STen2::gen_dyad( a1 );
      const auto H_f2 = fkd * STen2::gen_id() + ( 1.0 - 3.0 * fkd ) * STen2::gen_dyad( a2 );

      return mu * std::pow( F.det(), - 2.0/3.0 ) * STen2::gen_DEV_part( STen2::gen_id(), C ) +
             2.0  * dfpsi1_dfE1 * H_f1 + 2.0  * dfpsi2_dfE2 * H_f2;
    }

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_ich ) const
    {
      constexpr double pt67 = 2.0 / 3.0;

      const auto C = STen2::gen_right_Cauchy_Green( F );
      
      const double val = mu * std::pow( F.det(), -pt67 );

      const double I1 = C.tr();
      const double I4 = C.VecMatVec( a1, a1 );
      const double I6 = C.VecMatVec( a2, a2 );

      const double fE1 = fkd * I1 + ( 1.0 - 3.0 * fkd ) * I4 - 1.0;
      const double fE2 = fkd * I1 + ( 1.0 - 3.0 * fkd ) * I6 - 1.0;

      const double dfpsi1_dfE1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
      const double dfpsi2_dfE2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

      const auto H_f1 = fkd * STen2::gen_id() + ( 1.0 - 3.0 * fkd ) * STen2::gen_dyad( a1 );
      const auto H_f2 = fkd * STen2::gen_id() + ( 1.0 - 3.0 * fkd ) * STen2::gen_dyad( a2 );

      const auto S_ich =  val * STen2::gen_DEV_part( STen2::gen_id(), C ) +
             2.0  * dfpsi1_dfE1 * H_f1 + 2.0  * dfpsi2_dfE2 * H_f2;

      // First PK stress
      P_ich = F * S_ich;

      const double d2fpsi1_dfE1 = fk1 * ( 1.0 + 2.0 * fk2 * fE1 * fE1 ) * std::exp( fk2 * fE1 * fE1 );
      const double d2fpsi2_dfE2 = fk1 * ( 1.0 + 2.0 * fk2 * fE2 * fE2 ) * std::exp( fk2 * fE2 * fE2 );

      // Elasticity tensor
      auto CC_ich = pt67 * val * I1 * STen4::gen_Ptilde( STen2::inverse(C) );

      const auto S_iso = val * STen2::gen_DEV_part(STen2::gen_id(), C );

      CC_ich.add_SymmOutProduct( -pt67, STen2::inverse(C), S_iso );

      CC_ich.add_OutProduct( 4.0 * d2fpsi1_dfE1, H_f1 );
      CC_ich.add_OutProduct( 4.0 * d2fpsi2_dfE2, H_f2 );

      return CC_ich;
    }

    virtual double get_energy( const Tensor2_3D &F ) const
    {
      const auto C = STen2::gen_right_Cauchy_Green( F );

      const double I1 = C.tr();
      const double detFm0d67 = std::pow( F.det(), -2.0/3.0 );

      const double I4 = C.VecMatVec( a1, a1 );
      const double I6 = C.VecMatVec( a2, a2 );

      const double fE1 = fkd * I1 + ( 1.0 - 3.0 * fkd ) * I4 - 1.0;
      const double fE2 = fkd * I1 + ( 1.0 - 3.0 * fkd ) * I6 - 1.0;

      return 0.5 * mu *( detFm0d67 * I1 - 3.0) + 
             0.5 * ( fk1 / fk2 ) * ( std::exp( fk2 * fE1 * fE1 ) - 1.0 ) + 
             0.5 * ( fk1 / fk2 ) * ( std::exp( fk2 * fE2 * fE2 ) - 1.0 );
    }

    virtual Vector_3 get_fibre_dir (const int &dir) const
    {
      if (dir == 0) return a1;
      else if (dir == 1) return a2;
      else
      {
        SYS_T::print_fatal("Error:MaterialModel_ich_GOH14, wrong fibre direction.");
        return Vector_3();
      }
    }

  private:
    const double mu, f1_the, f1_phi, f2_the, f2_phi;
    const double fk1, fk2, fkd;
    const Vector_3 a1, a2;
};

#endif
