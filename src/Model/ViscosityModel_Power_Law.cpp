#include "ViscosityModel_Power_Law.hpp"

ViscosityModel_Power_Law::ViscosityModel_Power_Law( const double &in_m,
                                                    const double &in_n )
: m( in_m ), n( in_n )
{}

ViscosityModel_Power_Law::ViscosityModel_Power_Law( const char * const &fname )
{
	hid_t h5file = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	HDF5_Reader * h5r = new HDF5_Reader( h5file );

	SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
			"Error: ViscosityModel_Power_Law constructor does not match h5 file.\n");

	m      = h5r -> read_doubleScalar("/", "m");
	n      = h5r -> read_doubleScalar("/", "n");

	delete h5r;
	H5Fclose(h5file);
}

ViscosityModel_Power_Law::~ViscosityModel_Power_Law()
{}

