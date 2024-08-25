#include "MaterialModel_Vorp03_Incompressible_Mixed.hpp"

MaterialModel_Vorp03_Incompressible_Mixed::MaterialModel_Vorp03_Incompressible_Mixed( const double &in_rho,
        const double &in_c1, const double &in_c2 )
: pt33( 1.0/3.0 ), mpt67( -2.0 * pt33 ), rho0(in_rho),
  c1(in_c1), c2(in_c2)
{}

MaterialModel_Vorp03_Incompressible_Mixed::MaterialModel_Vorp03_Incompressible_Mixed(
    const char * const &fname )
: pt33( 1.0/3.0 ), mpt67( -2.0 * pt33 )
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_Vorp03_Incompressible_Mixed constructor does not match h5 file.\n" );

  rho0 = h5r -> read_doubleScalar("/", "rho0");
  c1   = h5r -> read_doubleScalar("/", "c1");
  c2   = h5r -> read_doubleScalar("/", "c2");

  delete h5r; H5Fclose(h5file);
}

MaterialModel_Vorp03_Incompressible_Mixed::~MaterialModel_Vorp03_Incompressible_Mixed()
{}

void MaterialModel_Vorp03_Incompressible_Mixed::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_Vorp03_Incompressible_Mixed: \n");
  SYS_T::commPrint("\t  Density rho  = %e \n", rho0);
  SYS_T::commPrint("\t  Material constant 1 = %e \n", c1);
  SYS_T::commPrint("\t  Material constant 2 = %e \n", c2);
}

void MaterialModel_Vorp03_Incompressible_Mixed::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "rho0", rho0 );
    h5w -> write_doubleScalar("c1", c1);
    h5w -> write_doubleScalar("c2", c2);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

void MaterialModel_Vorp03_Incompressible_Mixed::get_PK(
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  const double detCm0d67 = std::pow(F.det(), 2.0 * mpt67);
  const double I2 = C.I2();
  const double coe = 2.0 * c1 + 4.0 * c2 * (detCm0d67 * I2 - 3);
  S.copy(Cinv);
  S.scale( mpt67 * I2 );
  S.AXPI( C.tr() );
  S.AXPY( -1.0, C );
  S.scale( coe * detCm0d67 );
  P.MatMult(F,S);
}

void MaterialModel_Vorp03_Incompressible_Mixed::get_PK_Stiffness( 
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S, Tensor4_3D &CC ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  Tensor2_3D C1, C2;
  Tensor2_3D I; I.gen_id();
  Tensor4_3D II; II.gen_symm_id();
  const double detCm0d67 = std::pow(F.det(), 2.0 * mpt67);
  const double I2 = C.I2();
  const double coe = 2.0 * c1 + 4.0 * c2 * (detCm0d67 * I2 - 3);
  S.copy(Cinv);
  S.scale( mpt67 * I2 );
  S.AXPI( C.tr() );
  S.AXPY( -1.0, C );
  C1.copy(S);
  S.scale( coe * detCm0d67 );
  P.MatMult(F,S);
  
  C2.copy(C);
  C2.scale(-1.0);
  C2.AXPI( C.tr() );

  CC.gen_zero();
  const double val1 = 8.0 * c2 * detCm0d67 * detCm0d67;
  const double val2 = 2.0 * mpt67 * detCm0d67 * coe;
  const double val3 = 2.0 * coe * detCm0d67;
  CC.add_OutProduct(val1, C1, C1);
  CC.add_OutProduct(val2, C2, Cinv);
  CC.add_OutProduct(val2, Cinv, C2);
  CC.add_OutProduct(val3 * 4.0 * I2 / 9.0, Cinv, Cinv);
  CC.add_SymmProduct(2.0 * val3 * pt33 * I2, Cinv, Cinv);
  CC.add_OutProduct(val3, I, I);
  CC.AXPY(-val3, II);
}

double MaterialModel_Vorp03_Incompressible_Mixed::get_strain_energy( 
    const Tensor2_3D &F ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  const double detCm0d67 = std::pow(F.det(), 2.0 * mpt67);
  const double I2_iso = detCm0d67 * C.I2();
  return c1 * (I2_iso - 3.0) + c2 * (I2_iso - 3) * (I2_iso - 3);
}

// EOF
