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
														const double &in_lambda, const double &in_n );

		ViscosityModel_Carreau( const char * const &fname = "viscosity_model.h5");

		virtual ~ViscosityModel_Carreau();

		private:
// ----------------------------------------------------------------------------
// mu_inf : viscosity as shear rate tends to infinity
// mu_0   : viscosity as shear rate tends to 0
// lambda : a time constant
// n      : an index shows the dependency of viscosity on shear rate
// ----------------------------------------------------------------------------
			double mu_inf, mu_0, lambda, n;
};

#endif
