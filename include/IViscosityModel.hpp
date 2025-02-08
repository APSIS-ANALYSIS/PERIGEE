#ifndef IVISCOSITYMODEL_HPP
#define IVISCOSITYMODEL_HPP
// ============================================================================
// IViscosityModel.hpp
//
// Interface for fluid viscosity model (including Newtonian and Non-Newtonian
// fluid models)
// ============================================================================
#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"
#include "SymmTensor2_3D.hpp"

class IViscosityModel
{
  public:
    IViscosityModel() = default;

    virtual ~IViscosityModel() = default;

    virtual void print_info() const = 0;

    virtual std::string get_model_name() const
    {
      const std::string output = "unknown";
      SYS_T::commPrint("Warning: IViscosityModel::get_model_name() is not implemented. \n");
      return output;
    }

    virtual void write_hdf5( const char * const &fname = "viscosity_model.h5" ) const
    {
      SYS_T::commPrint("Warning: IViscosityModel::write_hdf5() is not implemented. \n");
    }

    // virtual double get_mu( const double &D_xx, const double &D_yy, const double &D_zz,
    //    const double &D_yz, const double &D_xz, const double &D_xy ) const = 0;

		// virtual double get_dmu_dI1( const double &D_xx, const double &D_yy,
    //                             const double &D_zz ) const = 0;

    // virtual double get_dmu_dI2( const double &D_xx, const double &D_yy,
    //                             const double &D_zz, const double &D_yz,
    //                             const double &D_xz, const double &D_xy ) const = 0;

		// virtual double get_dmu_dI3( const double &D_xx, const double &D_yy,
    //                             const double &D_zz, const double &D_yz,
    //                             const double &D_xz, const double &D_xy ) const = 0;

    virtual double get_mu( const Tensor2_3D &grad_velo ) const = 0;

    virtual double get_dmu_dI1( const Tensor2_3D &grad_velo ) const = 0;

    virtual double get_dmu_dI2( const Tensor2_3D &grad_velo ) const = 0;

    virtual double get_dmu_dI3( const Tensor2_3D &grad_velo ) const = 0;
};

#endif
