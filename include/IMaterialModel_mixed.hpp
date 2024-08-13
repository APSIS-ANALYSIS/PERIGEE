#ifndef IMATERIALMODEL_MIXED_HPP
#define IMATERIALMODEL_MIXED_HPP
// ============================================================================
// IMaterialModel_mixed.hpp
// ============================================================================

#include "IMaterialModel_vol.hpp"
#include "SymmTensor4_3D.hpp"

class IMaterialModel_mixed
{
  public:
    IMaterialModel_mixed( std::unique_ptr<IMaterialModel_vol> in_vmodel ) 
      : vmodel(std::move(in_vmodel)) {};

    virtual ~IMaterialModel_mixed() = default;

    virtual Tensor2_3D get_PK_1st( const Tensor2_3D &F ) const
    {
      return F * get_PK_2nd(F);
    }

    virtual SymmTensor2_3D get_PK_2nd( const Tensor2_3D &F ) const = 0;

    double get_rho_0() const {return vmodel->get_rho_0();}
    
    double get_rho( const double &pp ) const {return vmodel->get_rho(pp);}

    double get_drho_dp( const double &pp ) const { return vmodel->get_drho_dp(pp);}

    double get_beta( const double &pp ) const {return vmodel->get_beta(pp);}

    double get_dbeta_dp( const double &pp ) const {return vmodel->get_dbeta_dp(pp);}

    double get_elastic_kappa() const {return vmodel->get_elastic_kappa();}

  protected:
    std::unique_ptr<IMaterialModel_vol> vmodel;
};

#endif
