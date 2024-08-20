#ifndef MATERIALMODEL_ICH_NEOHOOKEAN_HPP
#define MATERIALMODEL_ICH_NEOHOOKEAN_HPP

#include "IMaterialModel_ich.hpp"

class MaterialModel_ich_NeoHookean : public IMaterialModel_ich
{
  public:
    MaterialModel_ich_NeoHookean( const double &in_mu ) : mu( in_mu ) {}

    virtual ~MaterialModel_ich_NeoHookean() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_ich_NeoHookean: \n");
      SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
    }

    virtual std::string get_model_name() const {return std::string("Neo-Hookean");}

    virtual double get_elastic_mu() const {return mu;}

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {
      const auto CC = STen2::gen_right_Cauchy_Green(F);
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
      return 0.5 * mu * ( std::pow(F.det(), -2.0/3.0) * CC.tr() - 3.0 );
    }

    virtual Vector_3 get_fibre_dir (const int &dir) const
    {
      return Vector_3();
    }

  private:
    const double mu;
};

#endif
