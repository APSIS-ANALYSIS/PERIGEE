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

    virtual double get_mu( const SymmTensor2_3D &strain_rate ) const = 0;

    virtual double get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const = 0;

    virtual double get_dmu_dI2( const SymmTensor2_3D &strain_rate ) const = 0;

    virtual double get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const = 0;
};

#endif
