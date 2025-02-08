#ifndef VISCOSITYMODEL_NEWTONIAN_HPP
#define VISCOSITYMODEL_NEWTONIAN_HPP
// ============================================================================
// ViscosityModel_Newtonian.hpp
//
// Interface for Newtonian model
// ============================================================================
#include "IViscosityModel.hpp"

class ViscosityModel_Newtonian : public IViscosityModel
{
  public:
    ViscosityModel_Newtonian() = delete;

    ViscosityModel_Newtonian( const double &in_mu );

    ViscosityModel_Newtonian( const char * const &fname = "viscosity_model.h5");

    virtual ~ViscosityModel_Newtonian();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Newtonian";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "viscosity_model.h5" ) const;

    // virtual double get_mu( const double &D_xx, const double &D_yy,
    //                        const double &D_zz, const double &D_yz,
    //                        const double &D_xz, const double &D_xy ) const
    // {
    //   return mu;
    // }

    // virtual double get_mu( const Matrix_3x3 &grad_velo ) const
    // {
    //   return mu;
    // }

    // virtual double get_dmu_dI1( const double &D_xx, const double &D_yy,
    //                             const double &D_zz) const
    // {
    //   return 0.0;
    // }

    // virtual double get_dmu_dI1( const Matrix_3x3 &grad_velo ) const
    // {
    //   return 0.0;
    // }

    // virtual double get_dmu_dI2( const double &D_xx, const double &D_yy,
    //                             const double &D_zz, const double &D_yz,
    //                             const double &D_xz, const double &D_xy ) const
    // {
    //   return 0.0;
    // }

    // virtual double get_dmu_dI2( const Matrix_3x3 &grad_velo ) const
    // {
    //   return 0.0;
    // }

    // virtual double get_dmu_dI3( const double &D_xx, const double &D_yy,
    //                             const double &D_zz, const double &D_yz,
    //                             const double &D_xz, const double &D_xy ) const
    // {
    //   return 0.0;
    // }

    // virtual double get_dmu_dI3( const Matrix_3x3 &grad_velo) const
    // {
    //   return 0.0;
    // }

    virtual double get_mu( const Tensor2_3D &grad_velo ) const
    {
      return mu;
    }

    virtual double get_dmu_dI1( const Tensor2_3D &grad_velo ) const
    {
      return 0.0;
    }

    virtual double get_dmu_dI2( const Tensor2_3D &grad_velo) const
    {
      return 0.0;
    }

    virtual double get_dmu_dI3( const Tensor2_3D &grad_velo ) const
    {
      return 0.0;
    }
         
    private:
    // ------------------------------------------------------------------------
    // mu : viscosity
    // ------------------------------------------------------------------------
    double mu;
};

#endif
