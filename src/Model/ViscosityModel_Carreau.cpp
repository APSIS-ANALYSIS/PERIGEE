#include "ViscosityModel_Carreau.hpp"

ViscosityModel_Carreau::ViscosityModel_Carreau( const double &in_mu_inf,
const double &in_mu_0, const double &in_lambda, const double &in_n )
: mu_inf( in_mu_inf ), mu_0( in_mu_0 ), lambda( in_lambda ), n( in_n )
{}

ViscosityModel_Carreau::ViscosityModel_Carreau( const char * const &fname )
{
	hid_t h5file = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	HDF5_Reader * h5r = new HDF5_Reader( h5file );

	SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
			"Error: ViscosityModel_Carreau constructor does not match h5 file.\n");

	mu_inf = h5r -> read_doubleScalar("/", "mu_inf");
	mu_0   = h5r -> read_doubleScalar("/", "mu_0");
	lambda = h5r -> read_doubleScalar("/", "lambda");
	n      = h5r -> read_doubleScalar("/", "n");

	delete h5r;
	H5Fclose(h5file);
}

ViscosityModel_Carreau::~ViscosityModel_Carreau()
{}

