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
        const double &in_fkd) : mu( in_mu ), 
    f1_the(in_f1the * MATH_T::PI / 180.0), f1_phi(in_f1phi * MATH_T::PI / 180.0), 
    f2_the(in_f2the * MATH_T::PI / 180.0), f2_phi(in_f2phi * MATH_T::PI / 180.0) , 
    fk1(in_fk1), fk2(in_fk2), fkd(in_fkd),
    a1( sin(f1_the) * cos(f1_phi), sin(f1_the) * sin(f1_phi), cos(f1_the) ),
    a2( sin(f2_the) * cos(f2_phi), sin(f2_the) * sin(f2_phi), cos(f2_the) ) {}

    virtual ~MaterialModel_ich_GOH06() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_ich_GOH06: \n");
      SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
      SYS_T::commPrint("\t  Dispersion parameter kappa   = %e \n", fkd);
      SYS_T::commPrint("\t  Stress-like parameter k1   = %e \n", fk1);
      SYS_T::commPrint("\t  Dimensionless parameter k2   = %e \n", fk2);
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

      const auto H_f1 = fkd * STen2::gen_id() + (1 - 3*fkd) * STen2::gen_dyad(a1);
      const auto H_f2 = fkd * STen2::gen_id() + (1 - 3*fkd) * STen2::gen_dyad(a2);

      // dPsi_fi/dC_tilde = dPsi_fi/dEi_tilde * H_fi
      // S_fi = 2 * dPsi_fi/dC = 2 * dPsi_fi/dC_tilde : dC_tilde/dC = 2 * J^(-2/3) *dPsi_fi/dEi_tilde *  P: H_fi
      const auto S_fi1 = 2.0 * detFm0d67 * dfpsi1_dfE1 * STen2::gen_DEV_part(H_f1, CC );
      const auto S_fi2 = 2.0 * detFm0d67 * dfpsi2_dfE2 * STen2::gen_DEV_part(H_f2, CC );

      return S_iso + S_fi1 + S_fi2;
    }

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
        Tensor2_3D &P_iso ) const
    {
      constexpr pt67 = 2.0 / 3.0;
      const auto CC = STen2::gen_right_Cauchy_Green(F);
      const double Inva1 = CC.tr();
      const auto invCC = STen2::inverse(CC);
      const double detFm0d67 = std::pow(F.det(), -pt67);

      const auto S_iso = mu * detFm0d67 * STen2::gen_DEV_part(STen2::gen_id(), CC );

      P_iso = F * S_iso;  // Tensor2_3D

      // Holzapfel p255. eqn(6.168)
      // PKstiff_tilde = 4J^(-4/3) * d^2Psi_sio/dC_tilde^2 = 4J^(-4/3) * d ((0.5 * mu) * I)/ dC_tilde = 0
      // P_tilde = C^(-1)odot C^(-1) - 1/3 * C^(-1) otimes C^(-1)
      // 2/3*J^(-2/3)*S_tilde:C = 2/3*J^(-2/3)*0.5*mu*I:C=1/3*J^(-2/3)*mu*trC
      auto PKstiff_iso = pt67 * detFm0d67 * mu * CC.tr() * STen4::gen_Ptilde( invCC );

      // -2/3 * (C^(-1) otimes S_iso + S_iso otimes C^(-1) )
      PKstiff_iso.add_SymmOutProduct( -pt67, invCC, S_iso );

      const double Inva4_1 = CC.VecMatVec(a1, a1);
      const double Inva4_2 = CC.VecMatVec(a2, a2);

      const double fE1 = detFm0d67 * (fkd*Inva1 + (1.0 -3.0*fkd)*Inva4_1) - 1.0;
      const double fE2 = detFm0d67 * (fkd*Inva1 + (1.0 -3.0*fkd)*Inva4_2) - 1.0;

      const double dfpsi1_dfE1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1);
      const double dfpsi2_dfE2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2);

      const auto H_f1 = fkd * STen2::gen_id() + (1 - 3*fkd) * STen2::gen_dyad(a1);
      const auto H_f2 = fkd * STen2::gen_id() + (1 - 3*fkd) * STen2::gen_dyad(a2);

      const auto S_fi1 = 2.0 * detFm0d67 * dfpsi1_dfE1 * STen2::gen_DEV_part(H_f1, CC );
      const auto S_fi2 = 2.0 * detFm0d67 * dfpsi2_dfE2 * STen2::gen_DEV_part(H_f2, CC );

      const auto S = S_iso + S_fi1 + S_fi2;

      //  dPsi_fi/dC_tilde = dPsi_fi/dEi_tilde * H_fi
      // d^2Psi_fi/dC_tilde^2 = d(dPsi_fi/dEi_tilde * H_fi) / dC_tilde = d^2Psi_fi / dEi_tilde^2 * H_fi otimes H_fi + 0
      const double d2fpsi1_dfE1 = fk1 * (1.0 + 2.0*fk2*fE1*fE1) * std::exp(fk2*fE1*fE1);
      const double d2fpsi2_dfE2 = fk1 * (1.0 + 2.0*fk2*fE2*fE2) * std::exp(fk2*fE2*fE2);

      // PKstiff_fi = 2 dS_fi / dC = 2d(2 * J^(-2/3) *dPsi_fi/dEi_tilde *  P: H_fi) / dC = 2(P:H_fi) otimes d(2J^(-2/3)dPsi_fi/dEi_tilde) /dC + 4J^(-2/3) dPsi_fi/dEi_tilde *d(P:H_fi) /dC
      // 2(P:H_fi) otimes d(2J^(-2/3)dPsi_fi/dEi_tilde) /dC = -4/3*J^(-2/3) *dPsi_fi/dEi_tilde * (P:H_fi) otimes C^(-1) + 4J^(-4/3)d^2Psi_fi/dEi_tilde^2* (P:H_fi) otimes (P:H_fi)
      // 4J^(-2/3) dPsi_fi/dEi_tilde * d(P:H_fi)/dC = -4/3*J^(-2/3)* dPsi_fi/dEi_tilde * C^(-1)otimes H_fi:P^T + 4/3*J^(-2/3)* dPsi_fi/dEi_tilde* (H_fi:C)C^(-1) odot C^(-1):P^T

      const double val1 = 2.0 * pt67 * detFm0d67 * dfpsi1_dfE1 * (fkd * Inva1 + (1.0-3.0*fkd)*Inva4_1);
      const double val2 = 2.0 * pt67 * detFm0d67 * dfpsi2_dfE2 * (fkd * Inva1 + (1.0-3.0*fkd)*Inva4_2);

      auto PKstiff_fi1 = val1 * STen4::gen_Ptilde( invCC );
      auto PKstiff_fi2 = val2 * STen4::gen_Ptilde( invCC );

      const double val = 4.0 * detFm0d67 * detFm0d67;

      PKstiff_fi1.add_OutProduct(val * d2fpsi1_dfE1, STen2::gen_DEV_part(H_f1, CC));
      PKstiff_fi2.add_OutProduct(val * d2fpsi2_dfE2, STen2::gen_DEV_part(H_f2, CC));


      PKstiff_fi1.add_SymmOutProduct(-pt67, invCC, S_fi1);
      PKstiff_fi2.add_SymmOutProduct(-pt67, invCC, S_fi2);

      return PKstiff_iso + PKstiff_fi1 + PKstiff_fi2;
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

    Vector_3 get_fibre_dir(const int &dir) const
    {
      if (dir == 0) return a1;
      else if (dir == 1) return a2;
      else
      {
        SYS_T::print_fatal("Error:MaterialModel_ich_GOH06, wrong fibre direction.");
        return Vector_3();
      }
    }

  private:
    const double mu;
    const double f1_the, f1_phi, f2_the, f2_phi;
    const double fk1, fk2, fkd;
    Vector_3 a1, a2;
};

#endif
