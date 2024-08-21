#ifndef MATERIALMODEL_ICH_GOH14_HPP
#define MATERIALMODEL_ICH_GOH14_HPP
// ============================================================================
// MaterialModel_ich_GOH14.hpp
//
// MaterialModel_GOH14_ST91_Mixed.hpp
// Quasi-incompressible Gasser-Odgen-Holzapfel model. The model's volumetric
// part is based on the Simo-Taylor91 penalization. The isochoric isotropic
// part is the same as the GOH06_Incompressible. The anisotropic part is
// modified by replacing the modified right Cauchy-Green tensor by the right
// Cauchy-Green tensor.
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
    f2_phi(in_f2phi*MATH_T::PI / 180.0), fk1(in_fk1), fk2(in_fk2), fkd(in_fkd)
    {
      a1(0) = sin(f1_the) * cos(f1_phi);
      a1(1) = sin(f1_the) * sin(f1_phi);
      a1(2) = cos(f1_the);

      a2(0) = sin(f2_the) * cos(f2_phi);
      a2(1) = sin(f2_the) * sin(f2_phi);
      a2(2) = cos(f2_the);
    }

    virtual ~MaterialModel_ich_GOH14() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_ich_GOH14: \n");
      SYS_T::commPrint("\t  Shear modulus mu = %e \n", mu);
      SYS_T::commPrint("\t  Dispersion parameter kappa = %e \n", fkd);
      SYS_T::commPrint("\t  Stress-like parameter k1 = %e \n",   fk1);
      SYS_T::commPrint("\t  Dimensionless parameter k2 = %e \n", fk2);
      SYS_T::commPrint("\t  Mean angle Theta_1 (deg) of the 1st family of fibres = %e \n", f1_the*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Phi_1 (deg) of the 1st family of fibres = %e \n",   f1_phi*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Theta_1 (rad) of the 1st family of fibres = %e \n", f1_the);
      SYS_T::commPrint("\t  Mean angle Phi_1 (rad) of the 1st family of fibres = %e \n",   f1_phi);
      SYS_T::commPrint("\t  Mean angle Theta_2 (deg) of 2nd family of fibres = %e \n",     f2_the*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Phi_2 (deg) of the 2nd family of fibres = %e \n",   f2_phi*180/MATH_T::PI);
      SYS_T::commPrint("\t  Mean angle Theta_2 (rad) of 2nd family of fibres = %e \n",     f2_the);
      SYS_T::commPrint("\t  Mean angle Phi_2 (rad) of the 2nd family of fibres = %e \n",   f2_phi);
      SYS_T::commPrint("\t  Mean direction a1 of the 1st family of fibres = [%e, %e, %e] \n", a1(0), a1(1), a1(2));
      SYS_T::commPrint("\t  Mean direction a2 of the 2nd family of fibres = [%e, %e, %e] \n", a2(0), a2(1), a2(2));
    }

    virtual std::string get_model_name() const {return std::string("GOH14");}

    virtual double get_elastic_mu() const {return mu;}

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {

    }

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_iso ) const
    {

    }

    virtual double get_energy( const Tensor2_3D &F ) const
    {

    }

    virtual Vector_3 get_fibre_dir (const int &dir) const
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
    const double mu, f1_the, f1_phi, f2_the, f2_phi;
    const double fk1, fk2, fkd;

    Vector_3 a1, a2;
};

#endif
