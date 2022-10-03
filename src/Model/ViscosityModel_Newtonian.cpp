#include "ViscosityModel_Newtonian.hpp"

ViscosityModel_Newtonian::ViscosityModel_Newtonian( const double &in_mu )
: mu( in_mu )
{}

ViscosityModel_Newtonian::ViscosityModel_Newtonian( const char * const &fname ) 
{
  hid_t h5file = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
      "Error: ViscosityModel_Newtonian constructor does not match h5 file.\n");

  mu = h5r -> read_doubleScalar("/", "mu");

  delete h5r;
  H5Fclose(h5file);
}

ViscosityModel_Newtonian::~ViscosityModel_Newtonian()
{}

