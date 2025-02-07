#ifndef VISCOSITYMODEL_CARREAU_HPP
#define VISCOSITYMODEL_CARREAU_HPP
// ============================================================================
// ViscosityModel_Carreau.hpp
//
// Interface for Carreau non-Newtonian model
// ============================================================================
#include "IViscosityModel.hpp"

class ViscosityModel_Carreau : public IViscosityModel
{
  public:
    ViscosityModel_Carreau() = delete;

    ViscosityModel_Carreau( const double &in_mu_inf, const double &in_mu_0,
                            const double &in_lambda, const double &in_n_pli );

    ViscosityModel_Carreau( const char * const &fname = "viscosity_model.h5");

    virtual ~ViscosityModel_Carreau();
    
    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Carreau";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "viscosity_model.h5" ) const;
    
    virtual double get_mu( const double &D_xx, const double &D_yy,
                           const double &D_zz, const double &D_yz,
                           const double &D_xz, const double &D_xy ) const;

    virtual double get_mu( const Matrix_3x3 &grad_velo ) const;

    virtual double get_dmu_dI1( const double &D_xx, const double &D_yy,
                                const double &D_zz) const;

    virtual double get_dmu_dI1( const Matrix_3x3 &grad_velo ) const;

    virtual double get_dmu_dI2( const double &D_xx, const double &D_yy,
                                const double &D_zz, const double &D_yz,
                                const double &D_xz, const double &D_xy ) const;

    virtual double get_dmu_dI2( const Matrix_3x3 &grad_velo ) const;

    virtual double get_dmu_dI3( const double &D_xx, const double &D_yy,
                                const double &D_zz, const double &D_yz,
                                const double &D_xz, const double &D_xy ) const;

    virtual double get_dmu_dI3( const Matrix_3x3 &grad_velo) const;

		private:
    // ----------------------------------------------------------------------------
    // mu_inf : viscosity as shear rate tends to infinity
    // mu_0   : viscosity as shear rate tends to 0
    // lambda : a time constant
    // n_pli  : Power law index, an index shows the dependency of viscosity on shear rate
    // ----------------------------------------------------------------------------
    double mu_inf, mu_0, lambda, n_pli;
};

#endif
