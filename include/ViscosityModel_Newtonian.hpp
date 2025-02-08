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

    virtual ~ViscosityModel_Newtonian() = default;

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Newtonian";
      return mname;
    }

    virtual double get_mu( const SymmTensor2_3D &strain_rate ) const
    {
      return mu;
    }

    virtual double get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const
    {
      return 0.0;
    }

    virtual double get_dmu_dI2( const SymmTensor2_3D &strain_rate) const
    {
      return 0.0;
    }

    virtual double get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const
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
