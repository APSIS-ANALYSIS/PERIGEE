#ifndef MATERIALMODEL_MIXED_NEOHOOKEAN_HPP
#define MATERIALMODEL_MIXED_NEOHOOKEAN_HPP

#include "IMaterialModel_mixed.hpp"

class MaterialModel_mixed_NeoHookean : public IMaterialModel_mixed
{
  public:
    MaterialModel_mixed_NeoHookean( std::unique_ptr<IMaterialModel_vol> vmodel,
        const double &in_mu ) : IMaterialModel_mixed(std::move(vmodel)), mu( in_mu ) {}

    virtual ~MaterialModel_mixed_NeoHookean() = default;

    virtual void print_info() const
    {}

    virtual double get_elastic_mu() const {return mu;}

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {
      const double val = mu * std::pow(F.det(), -2.0/3.0);
      const SymmTensor2_3D CC = STen2::gen_left_Cauchy_Green(F);
      Tensor2_3D out = Ten4::gen_P( CC.convert_to_full() ) * Ten2::gen_id();
      return val * STen2::gen_symm_part(out);
    }

    virtual SymmTensor4_3D get_PK_Stiffness( const Tensor2_3D &F,
       Tensor2_3D &P_iso ) const
    {
      constexpr double pt67 = 2.0 / 3.0;
      const double val = mu * std::pow(F.det(), -pt67);
      const SymmTensor2_3D CC = STen2::gen_left_Cauchy_Green(F);
      const SymmTensor2_3D invCC = STen2::inverse(CC);
      const SymmTensor2_3D S_iso = val * STen2::gen_symm_part( Ten4::gen_P( CC.convert_to_full() ) * Ten2::gen_id() );
     
      // First PK stress 
      P_iso = F * S_iso;
     
      // Elasticity tensor 
      SymmTensor4_3D out = pt67 * val * CC.tr() * STen4::gen_Ptilde( invCC );
      
      out.add_SymmOutProduct(-pt67, invCC, S_iso);

      return out;
    }

  private:
    const double mu;
};

#endif
