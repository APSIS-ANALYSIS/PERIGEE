#include "MaterialModel_StVenant_Kirchhoff_Simo85.hpp"

MaterialModel_StVenant_Kirchhoff_Simo85::MaterialModel_StVenant_Kirchhoff_Simo85( 
    const double &in_E, const double &in_nu )
: E( in_E ), nu( in_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{
}

MaterialModel_StVenant_Kirchhoff_Simo85::MaterialModel_StVenant_Kirchhoff_Simo85(
    const char * const &fname):
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{ 
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_StVenant_Kirchhoff_Simo85 constructor does not match h5 file.\n" );

  E  = h5r -> read_doubleScalar("/", "E");
  nu = h5r -> read_doubleScalar("/", "nu");

}

MaterialModel_StVenant_Kirchhoff_Simo85::~MaterialModel_StVenant_Kirchhoff_Simo85()
{}


void MaterialModel_StVenant_Kirchhoff_Simo85::print_info() const
{
  SYS_T::commPrint( "\t  MaterialModel_StVenant_Kirchhoff_Simo85: \n");
  SYS_T::commPrint( "\t  Young's Modulus E  = %e \n", E);
  SYS_T::commPrint( "\t  Possion's ratio nu = %e \n", nu);
  SYS_T::commPrint( "\t  Lame coeff lambda  = %e \n", lambda);
  SYS_T::commPrint( "\t  Shear modulus mu   = %e \n", mu);
  SYS_T::commPrint( "\t  Bulk modulus kappa = %e \n", kappa);
}

void MaterialModel_StVenant_Kirchhoff_Simo85::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar("E", E);
    h5w -> write_doubleScalar("nu", nu);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}


void MaterialModel_StVenant_Kirchhoff_Simo85::get_PK( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S )
{
  C.MatMultTransposeLeft(F);
  detF = F.det();
  S.copy(C);
  S.inverse();
  S.scale(kappa * std::log(detF));
  S.AXPY(mu, C);
  S.AXPY(-1.0 * mu, I);
  P.MatMult(F,S);
}


void MaterialModel_StVenant_Kirchhoff_Simo85::get_PK_Stiffness( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC)
{
  CC.gen_zero();
  const double val = -2.0 * kappa * std::log(detF);
  
  C.MatMultTransposeLeft(F);
  detF = F.det();
  S.copy(C);
  S.inverse();
  
  // Use S now because it is still C^-1 now. 
  CC.add_SymmProduct(val, S, S);
  CC.add_SymmProduct(2.0 * mu, I, I);
  CC.add_OutProduct( kappa, S, S);

  S.scale(kappa * std::log(detF));
  S.AXPY(mu, C);
  S.AXPY(-1.0 * mu, I);
  P.MatMult(F,S);
}


double MaterialModel_StVenant_Kirchhoff_Simo85::get_strain_energy( 
    const Matrix_3x3 &F )
{
  C.MatMultTransposeLeft(F);
  C.AXPY(-1.0, I);
  C.scale(0.5);
  C.MatMult(C,C);
  const double trE2 = C.tr();
  
  detF = F.det();

  return 0.5 * kappa * std::log(detF) * std::log(detF) + mu * trE2;
}

// EOF
