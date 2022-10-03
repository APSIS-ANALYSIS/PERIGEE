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

    ViscosityModel_Power_Law( const double &in_m, const double &in_n );

		ViscosityModel_Power_Law( const char * const &fname = "viscosity_model.h5");

		virtual ~ViscosityModel_Power_Law();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Power Law";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "viscosity_model.h5" ) const;

		private:
// ----------------------------------------------------------------------------
// m      : consistency ( Pa * s ^ n ) when n = 1, it is the same as viscosity
// n      : an index shows the dependency of viscosity on shear rate
// ----------------------------------------------------------------------------
			double m, n;
};

#endif
