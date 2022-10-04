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

void ViscosityModel_Newtonian::print_info() const
{
  SYS_T::commPrint("\t  ViscosityModel_Newtonian:: \n");
  SYS_T::commPrint("\t  Viscosity mu = %e \n", mu);
}

void ViscosityModel_Newtonian::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    HDF5_Writer * h5w = new HDF5_Writer( file_id );

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "mu", mu );

    delete h5w;
    H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

