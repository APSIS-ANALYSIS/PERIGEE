#include "MaterialModel_StVenant_Kirchhoff_M94_Mixed.hpp"

MaterialModel_StVenant_Kirchhoff_M94_Mixed::MaterialModel_StVenant_Kirchhoff_M94_Mixed( const double &in_E, const double &in_nu )
: rho0(1.0), E( in_E ), nu( in_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}

MaterialModel_StVenant_Kirchhoff_M94_Mixed::MaterialModel_StVenant_Kirchhoff_M94_Mixed( const double &in_rho0, const double &in_E, const double &in_nu )
: rho0( in_rho0 ), E( in_E ), nu( in_nu ), 
  lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}

MaterialModel_StVenant_Kirchhoff_M94_Mixed::MaterialModel_StVenant_Kirchhoff_M94_Mixed(
	const char * const &fname)
: pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ), I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_StVenant_Kirchhoff_M94_Mixed constructor does not match h5 file.\n" );

  rho0   =  h5r -> read_doubleScalar("/", "rho0");
  E      =  h5r -> read_doubleScalar("/", "E");
  nu     =  h5r -> read_doubleScalar("/", "nu");
  lambda =  h5r -> read_doubleScalar("/", "lambda");
  mu     =  h5r -> read_doubleScalar("/", "mu");
  kappa  =  h5r -> read_doubleScalar("/", "kappa");

  delete h5r; H5Fclose(h5file);
}

MaterialModel_StVenant_Kirchhoff_M94_Mixed::~MaterialModel_StVenant_Kirchhoff_M94_Mixed()
{}

void MaterialModel_StVenant_Kirchhoff_M94_Mixed::print_info() const
{
  SYS_T::commPrint( "\t  MaterialModel_StVenant_Kirchhoff_M94_Mixed: \n");
  SYS_T::commPrint( "\t  Young's Modulus E  = %e \n", E);
  SYS_T::commPrint( "\t  Possion's ratio nu = %e \n", nu);
  SYS_T::commPrint( "\t  Lame coeff lambda  = %e \n", lambda);
  SYS_T::commPrint( "\t  Shear modulus mu   = %e \n", mu);
  SYS_T::commPrint( "\t  Bulk modulus kappa = %e \n", kappa);
  SYS_T::commPrint( "\t  rho_0              = %e \n", rho0);
}

void MaterialModel_StVenant_Kirchhoff_M94_Mixed::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "rho0", rho0 );
    h5w -> write_doubleScalar("E", E);
    h5w -> write_doubleScalar("nu", nu);
    h5w -> write_doubleScalar("lambda", lambda);
    h5w -> write_doubleScalar("mu", mu);
    h5w -> write_doubleScalar("kappa", kappa);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

void MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_PK( 
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  
  const double detFm0d67 = std::pow( F.det(), mpt67 );

  const double mpt33CdC = (-1.0) * pt33 * C.MatContraction(C);
  Tensor2_3D PxC(Cinv); PxC.scale(mpt33CdC); PxC += C;
  
  Tensor2_3D PxI(Cinv); PxI.scale( (-1.0) * pt33 * C.tr() ); PxI.AXPI( 1.0 );

  S.copy(PxC); S.scale(detFm0d67);  S.AXPY( -1.0, PxI );
  S.scale( mu*detFm0d67 );

  P.MatMult(F,S);
}

void MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_PK_Stiffness( 
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S, Tensor4_3D &CC ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  
  const double detFm0d67 = std::pow(F.det(), mpt67);
  const double CdC = C.MatContraction(C);
  const double mpt33CdC = (-1.0) * pt33 * CdC;
  
  Tensor2_3D PxC(Cinv); PxC.scale(mpt33CdC); PxC += C;
  Tensor2_3D PxI(Cinv); PxI.scale( (-1.0) * pt33 * C.tr() ); PxI.AXPI( 1.0 );

  S.copy(PxC); S.scale(detFm0d67);  S.AXPY( -1.0, PxI );
  S.scale( mu*detFm0d67 );

  P.MatMult(F,S);

  CC.gen_symm_id();

  CC.add_OutProduct( (-1.0) * pt33, Cinv, C );
  CC.add_OutProduct( (-1.0) * pt33, C, Cinv );

  const double val1 = pt33 * pt33 * CdC;

  CC.add_OutProduct( val1, Cinv, Cinv );

  CC.scale( 2.0 * mu * detFm0d67 * detFm0d67 );

  const double val2 = 2.0 * pt33 * detFm0d67 * mu * ( detFm0d67 * CdC - C.tr() );

  CC.add_SymmProduct(val2, Cinv, Cinv);
  CC.add_OutProduct((-1.0)*pt33*val2, Cinv, Cinv);

  CC.add_OutProduct(mpt67, Cinv, S);
  CC.add_OutProduct(mpt67, S, Cinv);
}

double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_strain_energy( 
    const Tensor2_3D &F ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  const double detFm0d67 = std::pow(F.det(), mpt67);
  
  Tensor2_3D Es( C ); Es.scale(0.5 * detFm0d67); Es.AXPY(-0.5, I);

  Es.MatMult(Es,Es);

  return mu * Es.tr();
}

double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_rho( 
    const double &p ) const
{
  const double pk = p / kappa;

  return rho0 * (1.0 + pk);
}

double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_drho_dp( 
    const double &p ) const
{
  return rho0 / kappa;
}

double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_beta( 
    const double &p ) const
{
  return 1.0 / (p + kappa);
}

double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_dbeta_dp( 
    const double &p ) const
{
  return (-1.0) / ( (p+kappa) * (p+kappa) );
}

// EOF
