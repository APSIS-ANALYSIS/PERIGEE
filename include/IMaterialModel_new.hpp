#ifndef IMATERIALMODEL_NEW_HPP
#define IMATERIALMODEL_NEW_HPP
// ============================================================================
// IMaterialModel_new.hpp
// ============================================================================

#include "IMaterialModel_vol.hpp"

class IMaterialModel_new
{
  public:
    IMaterialModel_new( std::unique_ptr<IMaterialModel_vol> in_vmodel ) 
      : vmodel(std::move(in_vmodel)) {};

    virtual ~IMaterialModel_new() = default;

    double get_rho_0() const {return vmodel->get_rho_0();}
    
    double get_rho( const double &pp ) const {return vmodel->get_rho(pp);}

    double get_drho_dp( const double &pp ) const { return vmodel->get_drho_dp(pp);}

    double get_beta( const double &pp ) const {return vmodel->get_beta(pp);}

    double get_dbeta_dp( const double &pp ) const {return vmodel->get_dbeta_dp(pp);}

  protected:
    std::unique_ptr<IMaterialModel_vol> vmodel;
};

#endif
