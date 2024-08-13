#ifndef MATERIALMODEL_NEOHOOKEAN_HPP
#define MATERIALMODEL_NEOHOOKEAN_HPP

#include "IMaterialModel_mixed.hpp"

class MaterialModel_NeoHookean : public IMaterialModel_mixed
{
  public:
    MaterialModel_NeoHookean( std::unique_ptr<IMaterialModel_vol> vmodel,
        const double &in_mu ) : IMaterialModel_mixed(std::move(vmodel)), 
    elastic_mu( in_mu ) {}

    virtual ~MaterialModel_NeoHookean() = default;

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const
    {
      const double val = std::pow(F.det(), -2.0/3.0);
      const SymmTensor2_3D CC = STen2::gen_left_Cauchy_Green(F);
      Tensor2_3D out = Ten4::gen_P( CC.convert_to_full() ) * Ten2::gen_id();
      return elastic_mu * val * STen2::gen_symm_part(out);
    }

    virtual SymmTensor4_3D get_Stiffness( const Tensor2_3D &F ) const
    {
      constexpr pt67 = 2.0 / 3.0;
      const double val = elastic_mu * std::pow(F.det(), -pt67);
      const SymmTensor2_3D CC = STen2::gen_left_Cauchy_Green(F);
      const SymmTensor2_3D invCC = STen2::inverse(CC);
      
      SymmTensor4_3D out = pt67 * val * CC.tr() * STen4::gen_Ptilde( invCC );
      
      const SymmTensor2_3D S_iso = val * STen2::gen_symm_part( Ten4::gen_P( CC.convert_to_full() ) * Ten2::gen_id() );
      out.add_SymmOutProduct(-pt67, invCC, S_iso);

      return out;
    }

  private:
    double elastic_mu;
};

#endif
