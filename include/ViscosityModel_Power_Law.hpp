#ifndef VISCOSITYMODEL_POWER_LAW_HPP
#define VISCOSITYMODEL_POWER_LAW_HPP
// ============================================================================
// ViscosityModel_Power_Law.hpp
//
// Interface for Power Law non-Newtonian model
// ============================================================================
#include "IViscosityModel.hpp"

class ViscosityModel_Power_Law : public IViscosityModel
{
  public:
    ViscosityModel_Power_Law() = delete;
    
    ViscosityModel_Power_Law( const double &in_m_cons, const double &in_n_pli );

    ViscosityModel_Power_Law( const char * const &fname = "viscosity_model.h5");

    virtual ~ViscosityModel_Power_Law();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Power Law";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "viscosity_model.h5" ) const;

    // virtual double get_mu( const double &D_xx, const double &D_yy,
    //                        const double &D_zz, const double &D_yz,
    //                        const double &D_xz, const double &D_xy ) const;

    // virtual double get_mu( const Matrix_3x3 &grad_velo ) const;

    // virtual double get_dmu_dI1( const double &D_xx, const double &D_yy,
    //                             const double &D_zz) const;

    // virtual double get_dmu_dI1( const Matrix_3x3 &grad_velo ) const;

    // virtual double get_dmu_dI2( const double &D_xx, const double &D_yy,
    //                             const double &D_zz, const double &D_yz,
    //                             const double &D_xz, const double &D_xy ) const;

    // virtual double get_dmu_dI2( const Matrix_3x3 &grad_velo ) const;

    // virtual double get_dmu_dI3( const double &D_xx, const double &D_yy,
    //                             const double &D_zz, const double &D_yz,
    //                             const double &D_xz, const double &D_xy ) const;

    // virtual double get_dmu_dI3( const Matrix_3x3 &grad_velo) const;
    
    virtual double get_mu( const Tensor2_3D &grad_velo ) const;

    virtual double get_dmu_dI1( const Tensor2_3D &grad_velo ) const;

    virtual double get_dmu_dI2( const Tensor2_3D &grad_velo ) const;

    virtual double get_dmu_dI3( const Tensor2_3D &grad_velo ) const;

    private:
    // ------------------------------------------------------------------------
    // m_cons : consistency ( Pa * s ^ n ) when n = 1, it is the same as viscosity
    // n_pli  : an index shows the dependency of viscosity on shear rate
    // ------------------------------------------------------------------------
    double m_cons, n_pli;
};

#endif
