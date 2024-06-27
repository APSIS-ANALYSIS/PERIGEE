#include "MaterialModel_Linear_Elasticity.hpp"

MaterialModel_Linear_Elasticity::MaterialModel_Linear_Elasticity(
    const char * const &fname)
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_Linear_Elasticity constructor does not match h5 file.\n" );

  modulus_E = h5r -> read_doubleScalar("/", "modulus_E");
  nu        = h5r -> read_doubleScalar("/", "nu");
  
  delete h5r; H5Fclose(h5file);
}

void MaterialModel_Linear_Elasticity::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_Linear_Elasticity: \n");
  SYS_T::commPrint("\t  Young's Modulus E     = %e \n", modulus_E);
  SYS_T::commPrint("\t  Possion's ratio nu    = %e \n", nu);
  SYS_T::commPrint("\t  Lame parameter lambda = %e \n", get_elastic_lambda());
  SYS_T::commPrint("\t  Lame parameter mu     = %e \n", get_elastic_mu());
}

void MaterialModel_Linear_Elasticity::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar("modulus_E", modulus_E);
    h5w -> write_doubleScalar("nu", nu);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

Tensor2_3D MaterialModel_Linear_Elasticity::get_Cauchy_stress( const Tensor2_3D &F ) const
{
  const double lambda = get_elastic_lambda();
  const double mu = get_elastic_mu();
  const double l2mu = lambda + 2.0 * mu;

  return Tensor2_3D(l2mu * F(0) + lambda * (F(4) + F(8)), mu * ( F(1) + F(3) ), mu * ( F(2) + F(6) ),
                    mu * ( F(1) + F(3) ), l2mu * F(4) + lambda * (F(0) + F(8)), mu * ( F(5) + F(7) ),
                    mu * ( F(2) + F(6) ), mu * ( F(5) + F(7) ), l2mu * F(8) + lambda * (F(0) + F(4)));
}

// EOF
