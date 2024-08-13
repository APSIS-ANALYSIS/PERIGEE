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

  private:
    double elastic_mu;

};

#endif
