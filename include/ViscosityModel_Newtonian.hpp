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

		private:
// ----------------------------------------------------------------------------
// mu : viscosity
// ----------------------------------------------------------------------------
			double mu;
};

#endif
