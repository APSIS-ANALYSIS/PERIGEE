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

    ViscosityModel_Newtonian( const double &in_mu ) : mu(in_mu) {};

    ~ViscosityModel_Newtonian() override = default;

    void print_info() const override
    {
      SYS_T::commPrint("\t  ViscosityModel_Newtonian:: \n");
      SYS_T::commPrint("\t  Viscosity mu = %e \n", mu);
    }

    std::string get_model_name() const override
    {
      return "Newtonian";
    }

    double get_mu( const SymmTensor2_3D &strain_rate ) const override
    {
      return mu;
    }

    double get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const override
    {
      return 0.0;
    }

    double get_dmu_dI2( const SymmTensor2_3D &strain_rate) const override
    {
      return 0.0;
    }

    double get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const override
    {
      return 0.0;
    }

  private:
    // ------------------------------------------------------------------------
    // mu : viscosity
    // ------------------------------------------------------------------------
    const double mu;
};

#endif
