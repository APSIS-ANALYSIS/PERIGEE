#ifndef MATERIALMODEL_ICH_GOH06_HPP
#define MATERIALMODEL_ICH_GOH06_HPP
// ============================================================================
// MaterialModel_ich_GOH06.hpp
// 
// Anisotropic model from T.C.Gasser & R.W.Ogden & G.A.Holzapfel. Journal of the Royal Society Interface 3:15-35, 2006. 
//
// Date: Aug. 19 2024
// Author: Ju Liu, Jiawei Luo
// Contact: liujuy@gmail.com
// ============================================================================

#include "IMaterialModel_ich.hpp"
#include "Math_Tools.hpp"

class MaterialModel_ich_GOH06 : public IMaterialModel_ich
{
  public:
    MaterialModel_ich_GOH06( const double &in_mu, 
    const double &in_f1the, const double &in_f1phi,
    const double &in_f2the, const double &in_f2phi,
    const double &in_fk1, const double &in_fk2,
    const double &in_fkd) : mu( in_mu ), f1_the(in_f1the* MATH_T::PI / 180.0), f1_phi(in_f1phi*MATH_T::PI / 180.0), f2_the(in_f2the * MATH_T::PI / 180.0), f2_phi(in_f2phi*MATH_T::PI / 180.0) , fk1(in_fk1), fk2(in_fk2), fkd(in_fkd) 
    {
      a1(0) = sin(f1_the) * cos(f1_phi);
      a1(1) = sin(f1_the) * sin(f1_phi);
      a1(2) = cos(f1_the);

      a2(0) = sin(f2_the) * cos(f2_phi);
      a2(1) = sin(f2_the) * sin(f2_phi);
      a2(2) = cos(f2_the);
    }

    virtual ~MaterialModel_ich_GOH06() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_ich_GOH06: \n");
      SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
      SYS_T::commPrint("\t  Dispersion parameter kappa   = %e \n", fkd);
      SYS_T::commPrint("\t  Dimensionless parameter k1   = %e \n", fk1);
      SYS_T::commPrint("\t  Stress-like parameter k2   = %e \n", fk2);
      SYS_T::commPrint("\t  Mean angle Theta_1 (deg) of the 1st family of fibres = %e \n", f1_the*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Phi_1 (deg) of the 1st family of fibres = %e \n", f1_phi*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Theta_1 (rad) of the 1st family of fibres = %e \n", f1_the);
      SYS_T::commPrint("\t  Mean angle Phi_1 (rad) of the 1st family of fibres = %e \n", f1_phi);
      SYS_T::commPrint("\t  Mean angle Theta_2 (deg) of 2nd family of fibres = %e \n", f2_the*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Phi_2 (deg) of the 2nd family of fibres = %e \n", f2_phi*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Theta_2 (rad) of 2nd family of fibres = %e \n", f2_the);
      SYS_T::commPrint("\t  Mean angle Phi_2 (rad) of the 2nd family of fibres = %e \n", f2_phi);
      SYS_T::commPrint("\t  Mean direction a1 of the 1st family of fibres = [%e, %e, %e] \n", a1(0), a1(1), a1(2));
      SYS_T::commPrint("\t  Mean direction a2 of the 2nd family of fibres = [%e, %e, %e] \n", a2(0), a2(1), a2(2));

    }

    virtual std::string get_model_name() const {return std::string("GOH06");}

    virtual double get_elastic_mu() const {return mu;}

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {
      const auto CC = STen2::gen_right_Cauchy_Green(F);
      const double Inva1 = CC.tr();
      const double detFm0d67 = std::pow(F.det(), -2.0/3.0);

      // dPsi_iso/dC_tilde = dPsi_tilde_iso/dInva1_tilde * dInva1_tilde / dC_tilde = (0.5 * mu) * I
      // S_sio = 2*dPsi_iso/dC = 2*dPsi_iso/dC_tilde : dC_tilde/dC= mu * J^(-2/3) * P:I
      const auto S_iso = mu * detFm0d67 * STen2::gen_DEV_part(STen2::gen_id(), CC );

      const double Inva4_1 = CC.VecMatVec(a1, a1);
      const double Inva4_2 = CC.VecMatVec(a2, a2);

      const double fE1 = detFm0d67 * (fkd*Inva1 + (1.0 -3.0*fkd)*Inva4_1) - 1.0;
      const double fE2 = detFm0d67 * (fkd*Inva1 + (1.0 -3.0*fkd)*Inva4_2) - 1.0;

      const double dfpsi1_dfE1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1);
      const double dfpsi2_dfE2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2);

      Tensor2_3D a1xa1, a2xa2;
      a1xa1.gen_outprod(a1);
      a2xa2.gen_outprod(a2);

      auto H_f1 = fkd * Ten2::gen_id() + (1 - 3*fkd) * a1xa1;
      auto H_f2 = fkd * Ten2::gen_id() + (1 - 3*fkd) * a2xa2;

      // dPsi_fi/dC_tilde = dPsi_fi/dEi_tilde * H_fi
      // S_fi = 2 * dPsi_fi/dC = 2 * dPsi_fi/dC_tilde : dC_tilde/dC = 2 * J^(-2/3) *dPsi_fi/dEi_tilde *  P: H_fi
      const auto S_fi1 = 2.0 * detFm0d67 * dfpsi1_dfE1 * STen2::gen_DEV_part(STen2::gen_symm_part(H_f1), CC );
      const auto S_fi2 = 2.0 * detFm0d67 * dfpsi2_dfE2 * STen2::gen_DEV_part(STen2::gen_symm_part(H_f2), CC );

      return S_iso + S_fi1 + S_fi2;









      auto out = STen2::gen_DEV_part( STen2::gen_id(), CC );
      return mu * std::pow(CC.det(), -1.0/3.0) * out;
    }

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_iso ) const
    {
      constexpr double pt67 = 2.0 / 3.0;
      const auto CC = STen2::gen_right_Cauchy_Green(F);
      const double val = mu * std::pow(CC.det(), -pt67 * 0.5);
      
      const auto S_iso = val * STen2::gen_DEV_part( STen2::gen_id(), CC );
     
      // First PK stress 
      P_iso = F * S_iso;
     
      // Elasticity tensor 
      const auto invCC = STen2::inverse(CC);
      auto out = pt67 * val * CC.tr() * STen4::gen_Ptilde( invCC );
      
      out.add_SymmOutProduct(-pt67, invCC, S_iso);

      return out;
    }

    virtual double get_energy( const Tensor2_3D &F ) const
    {
      const auto CC = STen2::gen_right_Cauchy_Green(F);
      const double Inva1 = CC.tr();
      const double detFm0d67 = std::pow(F.det(), -2.0/3.0);

      const double Inva4_1 = CC.VecMatVec(a1, a1);
      const double Inva4_2 = CC.VecMatVec(a2, a2);

      const double fE1 = detFm0d67 * (fkd*Inva1 + (1.0 -3.0*fkd)*Inva4_1) - 1.0;
      const double fE2 = detFm0d67 * (fkd*Inva1 + (1.0 -3.0*fkd)*Inva4_2) - 1.0;

      const double Psi_iso = 0.5 * mu *( detFm0d67 * Inva1 - 3.0);
      const double Psi_fi1 = 0.5 * (fk1 / fk2) * (std::exp(fk2*fE1*fE1) - 1.0);
      const double Psi_fi2 = 0.5 * (fk1 / fk2) * (std::exp(fk2*fE2*fE2) - 1.0);

      return Psi_iso + Psi_fi1 + Psi_fi2;
    }

  private:
    const double mu;
    double f1_the, f1_phi, f2_the, f2_phi;
    double fk1, fk2, fkd;
    Vector_3 a1, a2;
};

#endif
